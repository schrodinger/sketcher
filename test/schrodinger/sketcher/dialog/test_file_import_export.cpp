
#define BOOST_TEST_MODULE test_file_import_export

#include <algorithm>

#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/test/checkexceptionmsg.h"

using namespace schrodinger::sketcher;
using schrodinger::rdkit_extensions::Format;

BOOST_TEST_DONT_PRINT_LOG_VALUE(Format);

namespace
{

/**
 * @return the import formats offered for an interface that allows everything
 * and a Sketcher that an import would replace
 */
FormatList<Format> all_import_formats()
{
    return get_import_formats(InterfaceType::ATOMISTIC_OR_MONOMERIC,
                              MoleculeType::EMPTY, true);
}

bool contains(const std::vector<std::string>& vec, const std::string& str)
{
    return std::find(vec.begin(), vec.end(), str) != vec.end();
}

/**
 * @return the labels of every entry in the given format list, in order
 */
std::vector<std::string> labels_of(const FormatList<Format>& formats)
{
    std::vector<std::string> labels;
    for (const auto& [format, label, extensions] : formats) {
        labels.push_back(label);
    }
    return labels;
}

/**
 * @return the extensions of the entry with the given label
 */
std::vector<std::string> extensions_for(const FormatList<Format>& formats,
                                        const std::string& label)
{
    for (const auto& [format, cur_label, extensions] : formats) {
        if (cur_label == label) {
            return extensions;
        }
    }
    return {};
}

} // namespace

BOOST_AUTO_TEST_CASE(test_get_import_export_formats)
{
    // 13 atomistic molecule entries + 2 reaction entries + HELM + FASTA
    BOOST_TEST(all_import_formats().size() == 17);
    // Formats with both compressed and uncompressed extensions are listed
    // twice; formats with no extensions at all are omitted
    BOOST_TEST(get_standard_export_formats().size() == 14);
    BOOST_TEST(get_reaction_export_formats().size() == 2);
    BOOST_TEST(get_image_export_formats().size() == 2);
}

BOOST_AUTO_TEST_CASE(test_get_standard_export_formats_includes_sequence_formats)
{
    auto formats = get_standard_export_formats();
    auto has_format = [&formats](Format target) {
        return std::ranges::any_of(formats, [target](const auto& e) {
            return std::get<0>(e) == target;
        });
    };
    BOOST_TEST(has_format(Format::HELM));
    BOOST_TEST(has_format(Format::FASTA));
}

/**
 * On input, SMILES and CXSMILES are read by the same parser, so they share a
 * single "SMILES" entry that claims both sets of extensions. SMARTS likewise
 * covers CXSMARTS, but has no extensions of its own and so is not offered in
 * a file dialog.
 */
BOOST_AUTO_TEST_CASE(test_smiles_entry_covers_cxsmiles)
{
    auto formats = all_import_formats();
    auto labels = labels_of(formats);
    BOOST_TEST(contains(labels, "SMILES"));
    BOOST_TEST(!contains(labels, "Extended SMILES"));
    BOOST_TEST(!contains(labels, "Extended SMARTS"));

    auto smiles_exts = extensions_for(formats, "SMILES");
    BOOST_TEST(contains(smiles_exts, ".smi"));
    BOOST_TEST(contains(smiles_exts, ".smiles"));
    BOOST_TEST(contains(smiles_exts, ".cxsmi"));
    BOOST_TEST(contains(smiles_exts, ".cxsmiles"));

    // SMARTS has no extensions, so it can't appear in a file dialog, but it
    // must still be offered where the user supplies the format themselves
    BOOST_TEST(!contains(labels, "SMARTS"));
    auto mol_labels = std::vector<std::string>();
    for (const auto& [format, label] : get_mol_import_formats()) {
        mol_labels.push_back(label);
    }
    BOOST_TEST(contains(mol_labels, "SMARTS"));
    BOOST_TEST(contains(mol_labels, "SMILES"));
}

/**
 * Compressed extensions get their own "[compressed]" entry so that no single
 * row in the dialog lists both compressed and uncompressed extensions.
 */
BOOST_AUTO_TEST_CASE(test_compressed_extensions_are_split_out)
{
    auto formats = all_import_formats();

    auto smiles_exts = extensions_for(formats, "SMILES");
    BOOST_TEST(!contains(smiles_exts, ".smigz"));
    BOOST_TEST(!contains(smiles_exts, ".smi.gz"));
    auto compressed_smiles_exts =
        extensions_for(formats, "SMILES [compressed]");
    BOOST_TEST(contains(compressed_smiles_exts, ".smigz"));
    BOOST_TEST(contains(compressed_smiles_exts, ".smi.gz"));

    // ".maezst" ends in "zst" without a preceding dot, so it is easy to
    // misclassify as uncompressed
    auto maestro_exts = extensions_for(formats, "Maestro");
    BOOST_TEST(maestro_exts == std::vector<std::string>{".mae"});
    auto compressed_maestro_exts =
        extensions_for(formats, "Maestro [compressed]");
    BOOST_TEST(contains(compressed_maestro_exts, ".maegz"));
    BOOST_TEST(contains(compressed_maestro_exts, ".maezst"));
    BOOST_TEST(contains(compressed_maestro_exts, ".mae.zst"));

    // Formats with no compressed extensions get no second entry
    BOOST_TEST(!contains(labels_of(formats), "MOL2 [compressed]"));

    // Within the atomistic group, the compressed entries are collected at the
    // end rather than interleaved with the formats they duplicate
    auto atomistic = labels_of(get_import_formats(InterfaceType::ATOMISTIC,
                                                  MoleculeType::EMPTY, true));
    BOOST_TEST(atomistic == (std::vector<std::string>{
                                "MDL SD", "Maestro", "SMILES", "InChI", "MOL2",
                                "PDB", "XYZ", "Marvin Document", "ChemDraw XML",
                                "MDL RXN", "Reaction SMILES",
                                "MDL SD [compressed]", "Maestro [compressed]",
                                "SMILES [compressed]", "PDB [compressed]"}));
}

/**
 * Sequence formats lead on a monomeric interface and trail elsewhere, and are
 * dropped entirely when an import can't produce a monomeric structure.
 */
BOOST_AUTO_TEST_CASE(test_sequence_formats_are_ordered_by_interface)
{
    auto is_seq = [](const std::string& label) {
        return label == "HELM" || label == "FASTA";
    };

    // Monomer tab: sequence formats first
    auto monomeric = labels_of(get_import_formats(InterfaceType::MONOMERIC,
                                                  MoleculeType::EMPTY, true));
    BOOST_TEST(monomeric == (std::vector<std::string>{"HELM", "FASTA"}));

    // Atomistic tab: no sequence formats at all
    auto atomistic = labels_of(get_import_formats(InterfaceType::ATOMISTIC,
                                                  MoleculeType::EMPTY, true));
    BOOST_TEST(atomistic.size() == 15);
    BOOST_TEST(std::none_of(atomistic.begin(), atomistic.end(), is_seq));

    // Both: sequence formats last, since the atomistic formats dominate
    auto both = labels_of(all_import_formats());
    BOOST_TEST(both.size() == 17);
    BOOST_TEST(is_seq(both[both.size() - 2]));
    BOOST_TEST(is_seq(both.back()));

    // An import that adds to an existing atomistic structure can't bring in a
    // sequence, so the sequence formats drop out
    auto adding_to_atomistic = labels_of(get_import_formats(
        InterfaceType::ATOMISTIC_OR_MONOMERIC, MoleculeType::ATOMISTIC, false));
    BOOST_TEST(std::none_of(adding_to_atomistic.begin(),
                            adding_to_atomistic.end(), is_seq));
    // ...but replacing it can
    auto replacing_atomistic = labels_of(get_import_formats(
        InterfaceType::ATOMISTIC_OR_MONOMERIC, MoleculeType::ATOMISTIC, true));
    BOOST_TEST(std::any_of(replacing_atomistic.begin(),
                           replacing_atomistic.end(), is_seq));
}

BOOST_AUTO_TEST_CASE(test_get_importable_mol_types)
{
    // A single-type interface pins the answer regardless of anything else
    BOOST_TEST(get_importable_mol_types(InterfaceType::ATOMISTIC,
                                        MoleculeType::MONOMERIC,
                                        true) == InterfaceType::ATOMISTIC);
    BOOST_TEST(get_importable_mol_types(InterfaceType::MONOMERIC,
                                        MoleculeType::ATOMISTIC,
                                        true) == InterfaceType::MONOMERIC);

    // Replacing the contents frees us up to import anything
    BOOST_TEST(get_importable_mol_types(InterfaceType::ATOMISTIC_OR_MONOMERIC,
                                        MoleculeType::ATOMISTIC, true) ==
               InterfaceType::ATOMISTIC_OR_MONOMERIC);
    // Adding to the contents limits us to what's already there
    BOOST_TEST(get_importable_mol_types(InterfaceType::ATOMISTIC_OR_MONOMERIC,
                                        MoleculeType::MONOMERIC,
                                        false) == InterfaceType::MONOMERIC);
    // An empty Sketcher imposes no limit
    BOOST_TEST(get_importable_mol_types(InterfaceType::ATOMISTIC_OR_MONOMERIC,
                                        MoleculeType::EMPTY, false) ==
               InterfaceType::ATOMISTIC_OR_MONOMERIC);
}

/**
 * Every FASTA extension maps to Format::FASTA, which can't be read, so the
 * ".fasta" entry would be useless without this substitution.
 */
BOOST_AUTO_TEST_CASE(test_resolve_ambiguous_import_format)
{
    using schrodinger::rdkit_extensions::to_rdkit;

    BOOST_TEST(resolve_ambiguous_import_format(Format::FASTA) ==
               Format::FASTA_PEPTIDE);
    // Anything that isn't ambiguous is passed straight through
    BOOST_TEST(resolve_ambiguous_import_format(Format::SMILES) ==
               Format::SMILES);
    BOOST_TEST(resolve_ambiguous_import_format(Format::HELM) == Format::HELM);

    // The substitution is what makes a FASTA file readable at all
    std::string fasta = ">seq\nAGL\n";
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(to_rdkit(fasta, Format::FASTA),
                                    std::invalid_argument, "FASTA_PEPTIDE");
    BOOST_TEST(to_rdkit(fasta, resolve_ambiguous_import_format(Format::FASTA))
                   ->getNumAtoms() > 0);
}
