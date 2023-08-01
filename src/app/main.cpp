#include <cstdio>

#include <fmt/format.h>

int main(int argc, char** argv)
{
    printf(fmt::format("The answer to life, the universe, and everything is {}",
                       42)
               .c_str());
    return 0;
}
