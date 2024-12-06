// @copyright Schrodinger, LLC - All Rights Reserved

#include "schrodinger/rdkit_extensions/example.h"

#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/assign.hpp>
#include <boost/beast/core/detail/base64.hpp>
#include <boost/bimap.hpp>
#include <boost/filesystem.hpp>
#include <boost/functional/hash.hpp>
#include <boost/json.hpp>
#include <boost/noncopyable.hpp>
#include <boost/operators.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/irange.hpp>
#include <boost/serialization/strong_typedef.hpp>
#include <boost/shared_ptr.hpp>
#include <fmt/format.h>
#include <rdkit/GraphMol/CIPLabeler/CIPLabeler.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/ChemReactions/ReactionParser.h>
#include <rdkit/GraphMol/DetermineBonds/DetermineBonds.h>
#include <rdkit/GraphMol/DistGeomHelpers/Embedder.h>
#include <rdkit/GraphMol/FileParsers/MolWriters.h>
#include <rdkit/GraphMol/inchi.h>
#include <zstd.h>

namespace schrodinger
{
namespace rdkit_extensions
{

static std::string check_zstd(std::string&& input)
{
    auto cctx = ZSTD_createCCtx();
    auto cBuffSize = ZSTD_compressBound(input.size());
    std::string output(cBuffSize, '\0');
    ZSTD_compressCCtx(cctx, output.data(), cBuffSize, input.data(),
                      input.size(), 1);
    ZSTD_freeCCtx(cctx);
    return output;
}

bool dependency_test(const std::string& smiles)
{
    // fmt
    auto test = fmt::format("foo{}.{}.{}", 10, 1, 1);

    // zstd
    check_zstd("foo");

    // boost
    class Tag : public boost::totally_ordered<Tag>,
                boost::totally_ordered2<Tag, int>,
                boost::noncopyable
    {
      public:
        explicit Tag(){};
    };

    static const boost::bimap<char, std::string> test_map =
        boost::assign::list_of<boost::bimap<char, std::string>::relation>(
            'F', "Foo")('B', "Bar");

    boost::filesystem::exists(test);
    boost::trim(test);
    boost::algorithm::to_lower(test);
    boost::starts_with(test, "x");
    std::vector<std::string> tokens;
    boost::split(tokens, "foo,bar", boost::is_any_of(","));
    boost::hash<std::string_view> hasher;
    hasher(test);
    boost::combine(tokens, boost::irange(2));
    std::string byte_array(test.size(), '\0');
    boost::beast::detail::base64::decode(byte_array.data(), test.data(),
                                         test.size());
    boost::json::stream_parser p;
    boost::system::error_code ec;
    p.write(byte_array.data(), byte_array.size(), ec);

    // RDKit
    boost::shared_ptr<RDKit::RWMol> mol;
    mol.reset(RDKit::SmilesToMol(smiles, 0, true));
    RDKit::CIPLabeler::assignCIPLabels(*mol, 1000000);
    RDKit::DGeomHelpers::EmbedMolecule(*mol, 0, 1234);
    RDKit::determineBonds(*mol, false, 0);
    RDKit::ExtraInchiReturnValues aux;
    RDKit::MolToInchi(*mol, aux);
    boost::shared_ptr<RDKit::ChemicalReaction> rxn;
    rxn.reset(RDKit::RxnSmartsToChemicalReaction("C.C>>CC", nullptr));
    // broken on windows -- maeparser.lib(Writer.obj) : error LNK2019:
    // unresolved external symbol "int const boost::iostreams::zlib::no_flush"
    RDKit::MaeWriter::getText(*mol);

    return true;
}

} // namespace rdkit_extensions
} // namespace schrodinger
