#include "schrodinger/sketcher/example.h"

#include <cstdio>

int main(int argc, char** argv)
{
    printf("\nThe following dependencies work!!!\n");
    if (schrodinger::sketcher::boost_library_works()) {
        printf("- boost\n");
    }
    if (schrodinger::sketcher::fmt_library_works()) {
        printf("- fmt\n");
    }
    if (schrodinger::sketcher::maeparser_library_works()) {
        printf("- maeparser\n");
    }
    if (schrodinger::sketcher::qt6_library_works()) {
        printf("- qt6\n");
    }
    if (schrodinger::sketcher::rdkit_library_works()) {
        printf("- rdkit\n");
    }
    if (schrodinger::sketcher::zstd_library_works()) {
        printf("- zstd\n");
    }
    return 0;
}
