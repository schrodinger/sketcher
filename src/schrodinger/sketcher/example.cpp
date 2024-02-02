#include "schrodinger/sketcher/example.h"

#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/assign.hpp>
#include <boost/beast/core/detail/base64.hpp>
#include <boost/bimap.hpp>
#include <boost/filesystem.hpp>
#include <boost/functional/hash.hpp>
#include <boost/noncopyable.hpp>
#include <boost/operators.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/irange.hpp>
#include <boost/serialization/strong_typedef.hpp>
#include <boost/shared_ptr.hpp>
#include <fmt/format.h>
#include <zstd.h>

namespace schrodinger
{
namespace sketcher
{

bool boost_library_works()
{
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

    std::string test = "foo";
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
    return true;
}

bool fmt_library_works()
{
    return fmt::format("fmt {}.{}.{}", 10, 1, 1) == "fmt 10.1.1";
}

bool maeparser_library_works()
{
    return false;
}

bool qt6_library_works()
{
    return false;
}

bool rdkit_library_works()
{
    return false;
}

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

bool zstd_library_works()
{
    return check_zstd("foo").size() > 0;
}

} // namespace sketcher
} // namespace schrodinger
