/* -------------------------------------------------------------------------
 * Implements schrodinger::rdkit_extensions:: text block <-> rdkit
 mol conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#include "schrodinger/rdkit_extensions/convert.h"

#include <functional>

#include <boost/algorithm/string.hpp>
#include <boost/beast/core/detail/base64.hpp>

#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/ChemReactions/ReactionParser.h>
#include <rdkit/GraphMol/ChemReactions/ReactionPickler.h>
#include <rdkit/GraphMol/Chirality.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <rdkit/GraphMol/DetermineBonds/DetermineBonds.h>
#include <rdkit/GraphMol/DistGeomHelpers/Embedder.h>
#include <rdkit/GraphMol/FileParsers/CDXMLParser.h>
#include <rdkit/GraphMol/FileParsers/FileParsers.h>
#include <rdkit/GraphMol/FileParsers/MolFileStereochem.h>
#include <rdkit/GraphMol/FileParsers/MolSupplier.h>
#include <rdkit/GraphMol/FileParsers/MolWriters.h>
#include <rdkit/GraphMol/MarvinParse/MarvinParser.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/MolPickler.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/SmilesParse/SmartsWrite.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <rdkit/GraphMol/StereoGroup.h>
#include <rdkit/GraphMol/SubstanceGroup.h>
#include <rdkit/GraphMol/inchi.h>

#include <zstd.h>

#include "schrodinger/rdkit_extensions/capture_rdkit_log.h"
#include "schrodinger/rdkit_extensions/atomistic_conversions.h"
#include "schrodinger/rdkit_extensions/constants.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/rdkit_extensions/fasta/to_rdkit.h"
#include "schrodinger/rdkit_extensions/fasta/to_string.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"
#include "schrodinger/rdkit_extensions/helm/to_string.h"
#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/rgroup.h"
#include "schrodinger/rdkit_extensions/stereochemistry.h"

namespace schrodinger
{
namespace rdkit_extensions
{

namespace
{

template <auto CreateFn, auto ResetFn, auto FreeFn, typename Ptr_t>
class ZstdContextMgr : private boost::noncopyable
{
  public:
    ZstdContextMgr() = default;
    ~ZstdContextMgr()
    {
        FreeFn(m_ctx);
    }

    Ptr_t get()
    {
        return m_ctx;
    }
    void reset()
    {
        ResetFn(m_ctx, ZSTD_reset_session_only);
    }

  private:
    Ptr_t m_ctx = CreateFn();
};

using ZstdCompressionContext =
    ZstdContextMgr<ZSTD_createCCtx, ZSTD_CCtx_reset, ZSTD_freeCCtx, ZSTD_CCtx*>;

using ZstdDecompressionContext =
    ZstdContextMgr<ZSTD_createDCtx, ZSTD_DCtx_reset, ZSTD_freeDCtx, ZSTD_DCtx*>;

/**
 * @internal
 * Guess by attempting conversions, and returning the first valid parse
 */
template <typename T> boost::shared_ptr<T> auto_detect(
    const std::string& text, const std::vector<Format>& formats,
    std::function<boost::shared_ptr<T>(const std::string&, const Format)>
        text_to_)
{
    for (const auto& fmt : formats) {
        try {
            return text_to_(text, fmt);
        } catch (const std::exception&) {
            continue;
        }
    }
    throw std::invalid_argument("Unable to determine format");
}

void attachment_point_dummies_to_molattachpt_property(RDKit::RWMol& rdk_mol)
{
    rdk_mol.updatePropertyCache(false);
    std::vector<RDKit::Atom*> dummies_to_remove;
    for (auto atom : rdk_mol.atoms()) {
        if (is_attachment_point_dummy(*atom)) {
            dummies_to_remove.push_back(atom);

            // This should find only one neighbor
            for (auto parent : rdk_mol.atomNeighbors(atom)) {
                int attachment_point_type{1};
                if (parent->hasProp(RDKit::common_properties::molAttachPoint)) {
                    attachment_point_type = -1;
                }
                parent->setNoImplicit(false);
                parent->setProp(RDKit::common_properties::molAttachPoint,
                                attachment_point_type);
            }
        }
    }
    rdk_mol.beginBatchEdit();
    for (auto atom : dummies_to_remove) {
        rdk_mol.removeAtom(atom);
    }
    rdk_mol.commitBatchEdit();

    rdk_mol.updatePropertyCache(false);
}

bool molattachpt_property_to_attachment_point_dummies(RDKit::RWMol& rdk_mol)
{
    std::unordered_set<int> new_attachment_dummies;
    for (auto& atom : rdk_mol.atoms()) {
        int attach_point_type{0};
        if (atom->getPropIfPresent(RDKit::common_properties::molAttachPoint,
                                   attach_point_type)) {
            auto num_dummies = (attach_point_type == -1 ? 2 : 1);
            auto explicit_h_count = atom->getNumExplicitHs();
            for (auto i = 1; i <= num_dummies; ++i) {

                auto dummy_atom = new RDKit::QueryAtom(0);
                dummy_atom->setQuery(RDKit::makeAtomNullQuery());

                dummy_atom->setProp(RDKit::common_properties::atomLabel,
                                    ATTACHMENT_POINT_LABEL_PREFIX +
                                        std::to_string(i));

                bool update_label = false;
                bool take_ownership = true;
                auto dummy_idx =
                    rdk_mol.addAtom(dummy_atom, update_label, take_ownership);
                rdk_mol.addBond(atom, dummy_atom, RDKit::Bond::SINGLE);

                new_attachment_dummies.insert(dummy_idx);

                RDKit::MolOps::setTerminalAtomCoords(rdk_mol, dummy_idx,
                                                     atom->getIdx());
            }
            if (explicit_h_count != 0) {
                atom->setNumExplicitHs(explicit_h_count - num_dummies);
            }
            atom->clearProp(RDKit::common_properties::molAttachPoint);
        }
    }

    if (!new_attachment_dummies.empty()) {
        RDKit::Chirality::reapplyMolBlockWedging(rdk_mol);
        RDKit::MolOps::assignChiralTypesFromBondDirs(rdk_mol);
        rdk_mol.updatePropertyCache(false);
        return true;
    }
    return false;
}

void preserve_wiggly_bonds(RDKit::RWMol& mol)
{
    for (auto& bond : mol.bonds()) {
        int wiggly_bond_v2000{0};
        int wiggly_bond_v3000{0};
        bond->getPropIfPresent(RDKit::common_properties::_MolFileBondStereo,
                               wiggly_bond_v2000);
        bond->getPropIfPresent(RDKit::common_properties::_MolFileBondCfg,
                               wiggly_bond_v3000);
        if (wiggly_bond_v2000 == 4 || wiggly_bond_v3000 == 2) {
            bond->setBondDir(RDKit::Bond::BondDir::UNKNOWN);
        }
    }
}

void fix_cxsmiles_rgroups(RDKit::ROMol& mol)
{
    for (auto atom : mol.atoms()) {
        if (atom->getAtomicNum() == 0) {

            // In SMILES RGroups have no query, but SMARTS & MDL
            // they should have an empty query
            if (atom->hasQuery()) {
                auto queryAtom = static_cast<RDKit::QueryAtom*>(atom);
                auto query_label = queryAtom->getQuery()->getTypeLabel();
                if (!query_label.empty()) {
                    continue;
                }
            }

            // If this is a null query, check for a proper RGroup atom label
            std::string atom_label;
            if (atom->getPropIfPresent(RDKit::common_properties::atomLabel,
                                       atom_label) &&
                atom_label.find("_R") == 0) {
                int idx = -1;
                size_t end = 0;
                try {
                    idx = std::stoi(atom_label.substr(2), &end);
                } catch (const std::invalid_argument&) {
                    continue;
                }
                if (idx > -1 && end == atom_label.size() - 2) {
                    // RDKit::setAtomRLabel() wouldn't work for R0
                    // (this will be corrected lated on)
                    atom->setProp(RDKit::common_properties::_MolFileRLabel,
                                  static_cast<unsigned>(idx));
                }
            }
        }
    }
}

/**
 * @internal Update the mol for extended SMILES/SMARTS output
 */
void adjust_for_extended_smiles_smarts_format(RDKit::ROMol& mol)
{
    // Don't pollute extension with coordinates
    mol.clearConformers();

    // Take any R-groups that are denoted using _MolFileRLabel and make sure
    // that the atom labels are set correctly. This ensures that R-group numbers
    // can be correctly exported to extended SMILES. This function also clears
    // any isotope values that are actually R-group numbers, as well as any
    // dummy label properties, so that those values don't appear in the SMILES.
    for (auto* atom : mol.atoms()) {
        auto r_group_num = get_r_group_number(atom);
        if (r_group_num.has_value()) {
            auto rlabel = R_GROUP_LABEL_PREFIX + std::to_string(*r_group_num);
            atom->setProp(RDKit::common_properties::atomLabel, rlabel);
            atom->clearProp(RDKit::common_properties::dummyLabel);
            if (atom->getIsotope() == *r_group_num) {
                atom->setIsotope(0);
            }
        }
    }
}

void adjust_for_extended_smiles_smarts_format(
    const RDKit::ChemicalReaction& rxn)
{
    for (auto reactant : rxn.getReactants()) {
        adjust_for_extended_smiles_smarts_format(*reactant);
    }
    for (auto product : rxn.getProducts()) {
        adjust_for_extended_smiles_smarts_format(*product);
    }
}

void fix_r0_rgroup(RDKit::ROMol& mol)
{
    int highest_rlabel{0};
    std::vector<RDKit::Atom*> atoms_with_r0;
    for (auto atom : mol.atoms()) {
        int rlabel{-1};
        if (atom->getPropIfPresent(RDKit::common_properties::_MolFileRLabel,
                                   rlabel)) {
            highest_rlabel = std::max(highest_rlabel, rlabel);
            if (rlabel == 0) {
                atoms_with_r0.push_back(atom);
            }
        }
    }
    if (!atoms_with_r0.empty()) {
        ++highest_rlabel;
        for (auto atom : atoms_with_r0) {
            RDKit::setAtomRLabel(atom, highest_rlabel);
        }
    }
}

void throw_parse_error(const CaptureRDErrorLog& rd_error_log,
                       const std::string& text)
{
    throw std::invalid_argument(rd_error_log.messages() +
                                "Failed to parse text:\n" + text);
}

template <typename T>
RDKit::RWMol* read_mol(T& reader, const CaptureRDErrorLog& rd_error_log)
{

    boost::shared_ptr<RDKit::ROMol> mol = nullptr;
    try {
        mol.reset(reader.next());
    } catch (std::runtime_error&) { // EOF hit
        return nullptr;
    }
    if (mol == nullptr) {
        return nullptr;
    }
    if (!reader.atEnd()) {
        throw std::invalid_argument(
            "Single structure required; multiple present" +
            rd_error_log.messages());
    }
    return new RDKit::RWMol(*mol);
}

void adjust_for_mdl_v2k_format(RDKit::RWMol& mol)
{
    auto nsgroups = getSubstanceGroups(mol).size();
    if (mol.getNumAtoms() > 999u || mol.getNumBonds() > 999u ||
        nsgroups > 999u) {
        throw std::runtime_error("Mol incompatible with MDL V2000 format");
    }
    // MDL V2000 format does not support enhanced stereo groups; if RDKit
    // detects then on write, it will force V3000 format to prevent
    // information loss; here, we explicitly clear enhanced stereo so that
    // V2000 can be written as requested.
    mol.setStereoGroups({});
}

// This will allow us to combine absolute stereo groups into one. e.g.
// |a:1,a:3| -> |a:1,3|
void combine_absolute_stereo_groups(RDKit::RWMol& mol)
{
    using namespace RDKit;

    // get copy of stereo groups
    auto stereo_groups = mol.getStereoGroups();

    // clear stereo groups
    mol.setStereoGroups({});

    // we only need to merge the atoms; other info can be accessed via the
    // original stereo group
    std::vector<Atom*>* combined_absolute_stereo_atoms = nullptr;
    std::vector<std::pair<StereoGroup*, std::vector<Atom*>>> combined;
    for (auto& sg : stereo_groups) {
        // always add OR and AND stereo groups
        if (sg.getGroupType() != StereoGroupType::STEREO_ABSOLUTE) {
            combined.push_back({&sg, sg.getAtoms()});
        }
        // this is the first ABS stereo group
        else if (combined_absolute_stereo_atoms == nullptr) {
            combined.push_back({&sg, sg.getAtoms()});
            combined_absolute_stereo_atoms = &combined.back().second;
        }
        // this is an ABS stereo group, but not the first one
        else {
            std::ranges::transform(
                sg.getAtoms(),
                std::back_inserter(*combined_absolute_stereo_atoms),
                [](auto& entry) { return entry; });
        }
    }

    // now construct the new stereo groups
    std::vector<StereoGroup> new_stereogroups;
    for (const auto& [sg, atoms] : combined) {
        new_stereogroups.push_back(
            {sg->getGroupType(), atoms, {}, sg->getReadId()});
        new_stereogroups.back().setWriteId(sg->getWriteId());
    }

    mol.setStereoGroups(std::move(new_stereogroups));
}

void set_xyz_title(RDKit::RWMol& mol)
{
    std::string title;
    mol.getPropIfPresent(RDKit::common_properties::_Name, title);

    int charge = 0;
    for (auto atom : mol.atoms()) {
        charge += atom->getFormalCharge();
    }
    if (charge != 0) {
        title += title.empty() ? "" : "; ";
        title += "charge=" + std::to_string(charge);
    }

    int multiplicity = 1;
    if (mol.getPropIfPresent("i_m_Spin_multiplicity", multiplicity) &&
        multiplicity != 1) {
        title += title.empty() ? "" : "; ";
        title += "multiplicity=" + std::to_string(multiplicity);
    }

    mol.setProp(RDKit::common_properties::_Name, title);
}

/**
 * Attempts to read the XYZ title for an encoded charge, on common convention
 * in Schrodinger usage of .xyz files. If no "charge=" string is present,
 * assumes total charge is 0
 */
[[maybe_unused]] int get_xyz_charge(const std::string& xyz_block)
{
    std::vector<std::string> lines;
    boost::split(lines, xyz_block, boost::is_any_of("\n"));
    if (lines.size() < 2) {
        return 0;
    }
    // title should be the second line of xyz input
    std::vector<std::string> title_tokens;
    boost::split(title_tokens, lines[1], boost::is_any_of("; "));
    for (const auto& token : title_tokens) {
        if (boost::starts_with(token, "charge=")) {
            try {
                return std::stoi(token.substr(7));
            } catch (const std::invalid_argument&) {
                return 0;
            }
        }
    }
    return 0;
}

std::string zstd_compress(std::string&& byte_array)
{
    static ZstdCompressionContext cctx;

    constexpr int zstd_compression_level = 1;

    const auto cBuffSize = ZSTD_compressBound(byte_array.size());
    if (ZSTD_isError(cBuffSize)) {
        // failed to calculate the compression bounds
        return byte_array;
    }

    std::string compressed_byte_array(cBuffSize, '\0');
    const auto cSize = ZSTD_compressCCtx(
        cctx.get(), compressed_byte_array.data(), cBuffSize, byte_array.data(),
        byte_array.size(), zstd_compression_level);

    if (ZSTD_isError(cSize)) {
        // Compression failed; reset the context and just return the
        // uncompressed array, in case it works?
        cctx.reset();
        return byte_array;
    } else if (cSize >= byte_array.size()) {
        // The compression worked, but not enough to offset the extra size due
        // to the header. If this is the case, prefer the uncompressed data,
        // which will be smaller and we won't have to decompress.
        return byte_array;
    }

    // trim the output to the actual compressed size
    compressed_byte_array.resize(cSize);

    return compressed_byte_array;
}

std::string zstd_decompress(std::string&& byte_array)
{
    static ZstdDecompressionContext dctx;

    const auto rSize =
        ZSTD_getFrameContentSize(byte_array.data(), byte_array.size());

    if (rSize == ZSTD_CONTENTSIZE_ERROR) {
        // Probably not compressed by zstd -- just return the byte array
        return byte_array;
    } else if (ZSTD_isError(rSize)) {
        // failed to read the uncompressed size. Maybe not compressed?
        return byte_array;
    }

    std::string decompressed_byte_array(rSize, '\0');
    const auto dSize =
        ZSTD_decompressDCtx(dctx.get(), decompressed_byte_array.data(), rSize,
                            byte_array.data(), byte_array.size());

    if (ZSTD_isError(dSize)) {
        // Decompression failed; reset the context and return the uncompressed
        // array, in case it works?
        dctx.reset();
        return byte_array;
    }

    return decompressed_byte_array;
}

template <typename T, typename U> boost::shared_ptr<T>
base64_to_mol(const std::string& text,
              std::function<void(const std::string&, U*)> from_pickle_func)
{
    // text is about 1/3 longer than the decoded data
    std::string byte_array(text.size(), '\0');
    const auto& [sz, read_sz] = boost::beast::detail::base64::decode(
        byte_array.data(), text.data(), text.size());
    byte_array.resize(sz);

    if (byte_array.empty()) {
        return nullptr; // failed decoding
    }

    const auto decompressed_byte_array = zstd_decompress(std::move(byte_array));
    boost::shared_ptr<T> mol_or_rxn(new T());
    try {
        from_pickle_func(decompressed_byte_array, mol_or_rxn.get());
    } catch (const std::exception&) {
        return nullptr; // PicklerException base class
    }
    return mol_or_rxn;
}

template <typename T> std::string mol_to_base64(
    const T* mol_or_rxn,
    std::function<void(const T*, std::string&, unsigned int)> pickle_func)
{
    std::string byte_array;
    pickle_func(mol_or_rxn, byte_array, RDKit::PicklerOps::AllProps);
    auto compressed_byte_array = zstd_compress(std::move(byte_array));

    // base64 has a 1/3 size overhead in respect to binary, but we give
    // it some extra space to be safe. Base64 is also padded, which makes
    // the minimum output size 4 bytes.
    constexpr size_t min_padded_size = 4;
    std::string text(
        std::max(2 * compressed_byte_array.size(), min_padded_size), '\0');
    const auto sz = boost::beast::detail::base64::encode(
        text.data(), compressed_byte_array.data(),
        compressed_byte_array.size());
    text.resize(sz);
    return text;
}

// A very simple check to see if a string can be parsed as SMILES without
// actually trying to parse it. We want this function because RDKit can
// parse things containing (e.g.) "[#6]" or "~" as a valid SMILES (!!!)
bool can_be_smiles(const std::string& text)
{
    for (auto c = text.begin(); c != text.end(); ++c) {
        if (*c == ' ' || *c == '\t') {
            // If the first thing we see is the separator between the
            // SMILES/SMARTS and the extension, it means any non-SMILES symbol
            // would be in the extension, so this can be a SMILES.
            return true;
        } else if (*c == '~') { // "Any" bond
            return false;
        } else if (*c == '[') { // e.g. [#6] or [13#6]
            for (++c; c != text.end() && *c >= '0' && *c <= '9'; ++c)
                ; // Skip isotope if present
            if (*c == '#') {
                // Not a SMILES either
                return false;
            }
        }
    }
    // We didn't see any of the forbidden symbols
    return true;
}

} // unnamed namespace

boost::shared_ptr<RDKit::RWMol> to_rdkit(const std::string& text,
                                         const Format format)
{
    if (format == Format::AUTO_DETECT) {
        return auto_detect<RDKit::RWMol>(
            text,
            {
                Format::RDMOL_BINARY_BASE64,
                Format::MDL_MOLV3000,
                Format::MAESTRO,
                Format::INCHI,
                Format::PDB,
                Format::MOL2,
                Format::XYZ,
                Format::MRV,
#ifndef __EMSCRIPTEN__
                // These formats don't parse correctly in WASM builds and may
                // crash the Sketcher.  This #ifndef should be removed as part
                // of SKETCH-2357.
                Format::CDXML,
#endif
                // Attempt SMILES before SMARTS, given not all SMARTS are SMILES
                Format::SMILES,
                Format::SMARTS,
                // Guess at HELM after guessing atomistic formats
                Format::HELM,
            },
            &to_rdkit);
    }

    CaptureRDErrorLog rd_error_log;

    boost::shared_ptr<RDKit::RWMol> mol;

    bool sanitize = false;
    bool removeHs = false;
    switch (format) {
        case Format::RDMOL_BINARY_BASE64:
            mol = base64_to_mol<RDKit::RWMol, RDKit::ROMol>(
                text, (void (*)(const std::string&, RDKit::ROMol*)) &
                          RDKit::MolPickler::molFromPickle);
            break;
        case Format::SMILES:
        case Format::EXTENDED_SMILES: {
            if (!can_be_smiles(text)) {
                throw std::invalid_argument(text + " is not a valid SMILES");
            }

            int debugParse = 0;
            mol.reset(RDKit::SmilesToMol(text, debugParse, sanitize));
            if (mol != nullptr) {
                fix_cxsmiles_rgroups(*mol);
            }

            break;
        }
        case Format::SMARTS:
        case Format::EXTENDED_SMARTS:
            mol.reset(RDKit::SmartsToMol(text));
            if (mol != nullptr) {
                fix_cxsmiles_rgroups(*mol);
            }
            break;
        case Format::MDL_MOLV2000:
        case Format::MDL_MOLV3000: {
            // SDMolSupplier will preserve structure level properties
            RDKit::SDMolSupplier reader;
            bool strict_parsing = false;
            reader.setData(text, sanitize, removeHs, strict_parsing);
            mol.reset(read_mol(reader, rd_error_log));
            break;
        }
        case Format::MAESTRO: {
            RDKit::MaeMolSupplier reader;
            try { // parse error on init()
                reader.setData(text, sanitize, removeHs);
            } catch (std::runtime_error&) {
                mol.reset();
                break;
            }
            mol.reset(read_mol(reader, rd_error_log));
            // workaround for SHARED-10515
            if (mol != nullptr && mol->hasProp("i_m_ct_stereo_status")) {
                mol->clearProp("i_m_ct_stereo_status");
            }
            break;
        }
        case Format::INCHI: {
            RDKit::ExtraInchiReturnValues aux;
            mol.reset(RDKit::InchiToMol(text, aux, sanitize, removeHs));
            break;
        }
        case Format::INCHI_KEY:
            throw std::invalid_argument("Cannot read from INCHI_KEY");
        case Format::PDB:
            mol.reset(RDKit::PDBBlockToMol(text, sanitize, removeHs));
            break;
        case Format::MOL2:
            mol.reset(RDKit::Mol2BlockToMol(text, sanitize, removeHs));
            break;
        case Format::XYZ:
            try {
                mol.reset(RDKit::XYZBlockToMol(text));
            } catch (const RDKit::FileParseException&) {
                break; // leave mol a nullptr to trigger throw_parse_error below
            }
#ifndef __EMSCRIPTEN__ // SKETCH-2080: RDKit::determineBonds malforms the WASM
            if (mol != nullptr) {
                // determineBonds requires all explicit hydrogens present; if
                // they are not, it may throw. We assume here that if a user
                // tries to read XYZ without explicit hydrogens, shame on them.
                // We also kekulize if possible given most users expect that.
                auto charge = get_xyz_charge(text);
                RDKit::determineBonds(*mol, /*use_huckel*/ false, charge);
                RDKit::MolOps::KekulizeIfPossible(*mol);
                rdkit_extensions::removeHs(*mol);
            }
#endif
            break;
        case Format::MRV:
            try {
                mol.reset(RDKit::MrvBlockToMol(text, sanitize, removeHs));
            } catch (const std::exception&) {
                mol.reset();
            }
            break;
        case Format::CDXML:
#ifndef __EMSCRIPTEN__
            // This format doesn't parse correctly in WASM builds and may crash
            // the Sketcher. See SKETCH-2357.
            mol.reset(new RDKit::RWMol());
            try { // combine multiple molecules into one like SmilesToMol
                for (auto& m : RDKit::CDXMLToMols(text, sanitize, removeHs)) {
                    mol->insertMol(*m);
                }
            } catch (const RDKit::FileParseException&) {
                mol.reset();
            }
            break;
#endif
        case Format::HELM:
            mol = helm::helm_to_rdkit(text);
            break;
        case Format::FASTA_PEPTIDE:
            mol = fasta::peptide_fasta_to_rdkit(text);
            break;
        case Format::FASTA_DNA:
            mol = fasta::dna_fasta_to_rdkit(text);
            break;
        case Format::FASTA_RNA:
            mol = fasta::rna_fasta_to_rdkit(text);
            break;
        case Format::FASTA:
            throw std::invalid_argument(
                "Can't determine FASTA format. Please use any of the "
                "FASTA_PEPTIDE, FASTA_DNA or FASTA_RNA formats to convert "
                "input.");
        default:
            throw std::invalid_argument("Unsupported import format");
    }

    if (mol == nullptr) {
        throw_parse_error(rd_error_log, text);
    }

    // SKETCH-2190: we don't want to sanitize SMARTS because they are not
    // complete molecules, and sanitization may, e.g. create radicals on
    // (query) atoms that do not have their valence completely satisfied.
    if (format != Format::SMARTS && format != Format::EXTENDED_SMARTS) {
        apply_sanitization(*mol, Sanitization::PARTIAL);
    } else {
        mol->updatePropertyCache(false);
    }

    bool has_attchpt = molattachpt_property_to_attachment_point_dummies(*mol);
    if (!has_attchpt) {
        preserve_wiggly_bonds(*mol);
        fix_r0_rgroup(*mol);
    }

    assign_stereochemistry(*mol);
    combine_absolute_stereo_groups(*mol);

    return mol;
}

boost::shared_ptr<RDKit::ChemicalReaction>
to_rdkit_reaction(const std::string& text, const Format format)
{
    if (format == Format::AUTO_DETECT) {
        return auto_detect<RDKit::ChemicalReaction>(
            text,
            {
                Format::RDMOL_BINARY_BASE64,
                Format::MDL_MOLV2000,
                // Attempt SMILES before SMARTS, given not all SMARTS are SMILES
                Format::SMILES,
                Format::SMARTS,
            },
            &to_rdkit_reaction);
    }

    CaptureRDErrorLog rd_error_log;

    boost::shared_ptr<RDKit::ChemicalReaction> rxn;

    bool sanitize = false;
    bool removeHs = false;
    switch (format) {
        case Format::RDMOL_BINARY_BASE64:
            rxn =
                base64_to_mol<RDKit::ChemicalReaction, RDKit::ChemicalReaction>(
                    text,
                    (void (*)(const std::string&, RDKit::ChemicalReaction*)) &
                        RDKit::ReactionPickler::reactionFromPickle);
            break;
        case Format::SMILES:
        case Format::EXTENDED_SMILES:
        case Format::SMARTS:
        case Format::EXTENDED_SMARTS: {
            auto useSMILES = (format == Format::SMILES ||
                              format == Format::EXTENDED_SMILES) &&
                             can_be_smiles(text);
            rxn.reset(
                RDKit::RxnSmartsToChemicalReaction(text, nullptr, useSMILES));

            // Bootstrap coords to preserve double bond stereo (SKETCH-1254)
            double spacing = 2.0;
            bool updateProps = true;
            RDDepict::compute2DCoordsForReaction(*rxn, spacing, updateProps);

            break;
        }
        case Format::MDL_MOLV2000:
        case Format::MDL_MOLV3000:
            rxn.reset(
                RDKit::RxnBlockToChemicalReaction(text, sanitize, removeHs));
            break;
        default:
            throw std::invalid_argument("Unsupported reaction import format");
    }
    if (rxn == nullptr) {
        throw_parse_error(rd_error_log, text);
    }

    for (auto reactant : rxn->getReactants()) {
        fix_r0_rgroup(*reactant);
        fix_cxsmiles_rgroups(*reactant);
    }

    for (auto product : rxn->getProducts()) {
        fix_r0_rgroup(*product);
        fix_cxsmiles_rgroups(*product);
    }

    return rxn;
}

std::string to_string(const RDKit::ROMol& input_mol, const Format format)
{
    CaptureRDErrorLog rd_error_log;

    auto mol = [&]() -> boost::shared_ptr<RDKit::RWMol> {
        // Since SHARED-11054, we want to start allowing conversion between
        // atomistic and monomeric mols
        auto is_monomeric = isMonomeric(input_mol);
        auto is_seq_format =
            std::ranges::find(SEQ_FORMATS, format) != SEQ_FORMATS.end();
        if (is_monomeric && !is_seq_format) {
            auto atomistic_mol = toAtomistic(input_mol);
            // NOTE: MaeWriter will attempt to generate 2D coordinates for this
            // molecule without the dummy conformer. The coordinate generation
            // procedure can end up taking a very long time for large
            // biomolecules
            atomistic_mol->addConformer(
                new RDKit::Conformer(atomistic_mol->getNumAtoms()));
            return atomistic_mol;
        } else if (!is_monomeric &&
                   (is_seq_format && format != Format::RDMOL_BINARY_BASE64)) {
            return toMonomeric(input_mol);
        } else {
            return boost::make_shared<RDKit::RWMol>(input_mol);
        }
    }();

    bool include_stereo = true;
    int confId = -1;
    bool kekulize = false;

    combine_absolute_stereo_groups(*mol);

    switch (format) {
        case Format::RDMOL_BINARY_BASE64:
            return mol_to_base64<RDKit::ROMol>(
                mol.get(),
                (void (*)(const RDKit::ROMol*, std::string&, unsigned int)) &
                    RDKit::MolPickler::pickleMol);
        case Format::SMILES:
            return RDKit::MolToSmiles(*mol, include_stereo, kekulize);
        case Format::EXTENDED_SMILES:
            adjust_for_extended_smiles_smarts_format(*mol);
            return RDKit::MolToCXSmiles(*mol, include_stereo, kekulize);
        case Format::SMARTS:
            return RDKit::MolToSmarts(*mol);
        case Format::EXTENDED_SMARTS:
            adjust_for_extended_smiles_smarts_format(*mol);
            return RDKit::MolToCXSmarts(*mol);
        case Format::MDL_MOLV2000:
        case Format::MDL_MOLV3000: {
            attachment_point_dummies_to_molattachpt_property(*mol);
            if (format == Format::MDL_MOLV2000) {
                adjust_for_mdl_v2k_format(*mol);
            }
            bool forceV3000 = format == Format::MDL_MOLV3000;
            return RDKit::SDWriter::getText(*mol, confId, kekulize, forceV3000);
        }
        case Format::MAESTRO:
            return RDKit::MaeWriter::getText(*mol);
        case Format::INCHI: {
            RDKit::ExtraInchiReturnValues aux;
            return RDKit::MolToInchi(*mol, aux);
        }
        case Format::INCHI_KEY:
            return RDKit::MolToInchiKey(*mol);
        case Format::PDB:
            return RDKit::MolToPDBBlock(*mol);
        case Format::XYZ:
            // Require explicit hydrogens and a complete 3D conformer on
            // write; see GraphMol/DetermineBonds/DetermineBonds.h in the
            // RDKit
            RDKit::MolOps::addHs(*mol);
            if (RDKit::DGeomHelpers::EmbedMolecule(*mol, /*maxIterations*/ 0,
                                                   /*seed*/ 1234) == -1) {
                throw std::runtime_error("Unable to export 3D coordinates");
            }
            set_xyz_title(*mol);
            return RDKit::MolToXYZBlock(*mol);
        case Format::MRV:
            return RDKit::MolToMrvBlock(*mol, include_stereo, confId, kekulize);
        case Format::HELM:
            return helm::rdkit_to_helm(*mol);
        case Format::FASTA:
            return fasta::rdkit_to_fasta(*mol);
        case Format::FASTA_PEPTIDE:
        case Format::FASTA_DNA:
        case Format::FASTA_RNA:
            throw std::invalid_argument("FASTA formats are auto-detected on "
                                        "write, so please use the FASTA format"
                                        "to write output.");
        default:
            throw std::invalid_argument("Unsupported export format");
    }

    throw std::invalid_argument("Invalid format specified");
}

std::string to_string(const RDKit::ChemicalReaction& rxn, const Format format)
{
    switch (format) {
        case Format::RDMOL_BINARY_BASE64:
            return mol_to_base64<RDKit::ChemicalReaction>(
                &rxn, (void (*)(const RDKit::ChemicalReaction*, std::string&,
                                unsigned int)) &
                          RDKit::ReactionPickler::pickleReaction);
        case Format::SMILES:
            return RDKit::ChemicalReactionToRxnSmiles(rxn);
        case Format::EXTENDED_SMILES:
            adjust_for_extended_smiles_smarts_format(rxn);
            return RDKit::ChemicalReactionToCXRxnSmiles(rxn);
        case Format::SMARTS:
            return RDKit::ChemicalReactionToRxnSmarts(rxn);
        case Format::EXTENDED_SMARTS:
            adjust_for_extended_smiles_smarts_format(rxn);
            return RDKit::ChemicalReactionToCXRxnSmarts(rxn);
        case Format::MDL_MOLV2000:
        case Format::MDL_MOLV3000: {
            bool separateAgents = false;
            bool forceV3000 = format == Format::MDL_MOLV3000;
            return RDKit::ChemicalReactionToRxnBlock(rxn, separateAgents,
                                                     forceV3000);
        }
        default:
            throw std::invalid_argument("Unsupported reaction export format");
    }
    throw std::invalid_argument("Invalid format specified");
}

} // namespace rdkit_extensions
} // namespace schrodinger
