#include "schrodinger/sketcher/example.h"

#include <string>

#include <fmt/format.h>
#include <zstd.h>

namespace schrodinger
{
namespace sketcher
{

bool boost_library_works()
{
    return false;
}

bool fmt_library_works()
{
    return fmt::format("fmt {}.{}.{}", 10, 1, 1).size() > 0;
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
