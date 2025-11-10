# ------------------------------------------------------------------------------
# This file and SQLite3Config.cmake.in are based on
# https://github.com/algoritnl/sqlite-amalgamation
#
# ------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.24...3.31)

# -- Project -------------------------------------------------------------------

# Extract the version info from the header file.
file(STRINGS sqlite3.h _ver_line
    REGEX "^#define SQLITE_VERSION +\"([0-9]+\\.[0-9]+\\.[0-9]+)\""
    LIMIT_COUNT 1)
string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+"
    SQLITE_VERSION ${_ver_line})
message(STATUS "SQLite library version: ${SQLITE_VERSION}")
unset(_ver_line)

project(SQLite3
    VERSION ${SQLITE_VERSION}
    DESCRIPTION "SQLite: an embedded SQL database engine"
    HOMEPAGE_URL https://sqlite.org
    LANGUAGES C)

if(PROJECT_IS_TOP_LEVEL)
    option(BUILD_SHARED_LIBS "Build using shared libraries." OFF)
    include(CTest)
endif()

# -- Configuration -------------------------------------------------------------

# Define a function for multiple choice configuration variables.
function(choice CHOICE DESCRIPTION)
    cmake_parse_arguments(CHOICE "FORCE_ON_ERROR;FATAL_ERROR" "DEFAULT" "VALID_VALUES" ${ARGN})
    set(${CHOICE} ${CHOICE_DEFAULT} CACHE STRING "${DESCRIPTION}")
    set_property(CACHE ${CHOICE} PROPERTY STRINGS ${CHOICE_VALID_VALUES})
    if(NOT ${CHOICE} IN_LIST CHOICE_VALID_VALUES)
        set(message "Configuration ${CHOICE}=${${CHOICE}} is not valid.")
        if(CHOICE_FATAL_ERROR)
            message(FATAL_ERROR "${message}")
        elseif(CHOICE_FORCE_ON_ERROR)
            string(APPEND message "\nContinuing with the default: ${CHOICE}=${CHOICE_DEFAULT}.")
            message(WARNING "${message}")
            set_property(CACHE ${CHOICE} PROPERTY VALUE ${CHOICE_DEFAULT})
        else()
            message(WARNING "${message}")
        endif()
    endif()
endfunction()

choice(SQLITE_DEFAULT_MEMSTATUS "Enable the SQLITE_CONFIG_MEMSTATUS features by default."
    VALID_VALUES 0 1 DEFAULT 1 FORCE_ON_ERROR)
choice(SQLITE_DEFAULT_SYNCHRONOUS "Determine the default value of PRAGMA synchronous."
    VALID_VALUES 0 1 2 3 DEFAULT 2 FORCE_ON_ERROR)
choice(SQLITE_DEFAULT_WAL_SYNCHRONOUS "Determine the default value of PRAGMA synchronous."
    VALID_VALUES 0 1 2 3 DEFAULT ${SQLITE_DEFAULT_SYNCHRONOUS} FORCE_ON_ERROR)
choice(SQLITE_DQS "Enable DQS: 0 = None (recommended); 1 = DML; 2 = DDL; 3 = Both (default)."
    VALID_VALUES 0 1 2 3 DEFAULT 3 FORCE_ON_ERROR)
choice(SQLITE_THREADSAFE "Select the thread model: 0 = single-threaded; 1 = serialized (default); 2 = multi-threaded."
    VALID_VALUES 0 1 2 DEFAULT 1 FATAL_ERROR)

set(SQLITE_MAX_ALLOCATION_SIZE -1 CACHE STRING "Set an upper bound on the size of memory allocations.")
set(SQLITE_MAX_EXPR_DEPTH -1 CACHE STRING "Set the maximum expression tree depth.")

option(SQLITE_ENABLE_API_ARMOR "This option activates extra code that attempts to detect misuse of the SQLite API." OFF)
option(SQLITE_ENABLE_COLUMN_METADATA "Enable some extra APIs that are required by some common systems, including Ruby-on-Rails." OFF)
option(SQLITE_ENABLE_DBSTAT_VTAB "This option enables the dbstat virtual table." OFF)
option(SQLITE_ENABLE_EXPLAIN_COMMENTS "This option adds extra logic to SQLite that inserts comment text into the output of EXPLAIN." OFF)
option(SQLITE_ENABLE_FTS3 "Enable version 3 of the full-text search engine." OFF)
option(SQLITE_ENABLE_FTS3_PARENTHESIS "Enable operators AND and NOT and allow nested parentheses." OFF)
option(SQLITE_ENABLE_FTS4 "Enable versions 3 and 4 of the full-text search engine." OFF)
option(SQLITE_ENABLE_FTS5 "Enable version 5 of the full-text search engine." OFF)
option(SQLITE_ENABLE_GEOPOLY "Enable the Geopoly extension." OFF)
option(SQLITE_ENABLE_ICU "Enable the ICU (International Components for Unicode) extension" OFF)
option(SQLITE_ENABLE_MATH_FUNCTIONS "Enable the built-in SQL math functions." ON)
option(SQLITE_ENABLE_NORMALIZE "Enable the sqlite3_normalized_sql() API." OFF)
option(SQLITE_ENABLE_PREUPDATE_HOOK "Enable the pre-update hook feature." OFF)
option(SQLITE_ENABLE_RBU "Enable the code that implements the RBU extension." OFF)
option(SQLITE_ENABLE_RTREE "This option causes SQLite to include support for the R*Tree index extension." OFF)
option(SQLITE_ENABLE_SESSION "Enable the session extension." OFF)
option(SQLITE_ENABLE_STMTVTAB "Enable the SQLITE_STMT virtual table logic." OFF)

option(SQLITE_LIKE_DOESNT_MATCH_BLOBS "Operators LIKE and GLOB always return FALSE if either operand is a BLOB." OFF)

option(SQLITE_OMIT_AUTOINIT "Omit calling sqlite3_initialize() automatically. The application is required to invoke it directly." OFF)
option(SQLITE_OMIT_COMPILEOPTION_DIAGS "Omit compile option interfaces." OFF)
option(SQLITE_OMIT_DECLTYPE "Omit the ability to return the declared type of columns from queries." OFF)
option(SQLITE_OMIT_DEPRECATED "Omit deprecated interfaces and features." OFF)
option(SQLITE_OMIT_DESERIALIZE "Omit sqlite3_serialize() and sqlite3_deserialize() interfaces." OFF)
option(SQLITE_OMIT_JSON "Omit the JSON SQL functions." OFF)
option(SQLITE_OMIT_LOAD_EXTENSION "Omit the extension loading mechanism." OFF)
option(SQLITE_OMIT_PROGRESS_CALLBACK "Omit the capability to issue progress callbacks." OFF)
option(SQLITE_OMIT_SHARED_CACHE "Omit support for shared-cache mode." OFF)

option(SQLITE_STRICT_SUBTYPE "This option helps to identify problems with application-defined SQL functions." ON)

option(SQLITE_USE_ALLOCA "Use alloca for memory allocations if possible." OFF)

if(SQLITE_ENABLE_COLUMN_METADATA AND SQLITE_OMIT_DECLTYPE)
    message(SEND_ERROR "Options SQLITE_ENABLE_COLUMN_META and SQLITE_OMIT_DECLTYPE are mutually exclusive.")
endif()

if(SQLITE_ENABLE_SESSION AND NOT SQLITE_ENABLE_PREUPDATE_HOOK)
    message(SEND_ERROR "Option SQLITE_ENABLE_SESSION requires option SQLITE_ENABLE_PREUPDATE_HOOK.")
endif()

if(NOT SQLITE_MAX_ALLOCATION_SIZE GREATER_EQUAL -1)
    message(SEND_ERROR "Configuration SQLITE_MAX_ALLOCATION_SIZE=${SQLITE_MAX_ALLOCATION_SIZE} is not valid.")
endif()

if(NOT SQLITE_MAX_EXPR_DEPTH GREATER_EQUAL -1)
    message(SEND_ERROR "Configuration SQLITE_MAX_EXPR_DEPTH=${SQLITE_MAX_EXPR_DEPTH} is not valid.")
endif()

# -- Introspection -------------------------------------------------------------

include(CMakePushCheckState)
cmake_push_check_state(RESET)

if(CMAKE_C_COMPILER_ID MATCHES "(GNU|Clang)")
    list(APPEND CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)
    list(APPEND CMAKE_REQUIRED_LIBRARIES m)
endif()

include(CheckIncludeFile)
include(CheckSymbolExists)
include(CheckTypeSize)

check_include_file(inttypes.h HAVE_INTTYPES_H)
check_include_file(stdint.h HAVE_STDINT_H)
check_type_size(int8_t INT8_T)
check_type_size(int16_t INT16_T)
check_type_size(uint8_t UINT8_T)
check_type_size(uint16_t UINT16_T)
check_type_size(uint32_t UINT32_T)

check_symbol_exists(posix_fallocate fcntl.h HAVE_POSIX_FALLOCATE)
check_symbol_exists(lstat "fcntl.h;sys/stat.h" HAVE_LSTAT)

check_include_file(malloc.h HAVE_MALLOC_H)
check_symbol_exists(malloc_usable_size malloc.h HAVE_MALLOC_USABLE_SIZE)

check_symbol_exists(isnan math.h HAVE_ISNAN)
check_symbol_exists(log2 math.h HAVE_LOG2)
check_symbol_exists(log10 math.h HAVE_LOG10)

check_symbol_exists(strchrnul string.h HAVE_STRCHRNUL)
check_symbol_exists(strerror_r string.h HAVE_STRERROR_R)

check_symbol_exists(gmtime_r time.h HAVE_GMTIME_R)
check_symbol_exists(localtime_r time.h HAVE_LOCALTIME_R)
check_symbol_exists(localtime_s time.h HAVE_LOCALTIME_S)
check_symbol_exists(nanosleep time.h HAVE_NANOSLEEP)

check_symbol_exists(fdatasync unistd.h HAVE_FDATASYNC)
check_symbol_exists(gethostuuid unistd.h HAVE_GETHOSTUUID)
check_symbol_exists(pread unistd.h HAVE_PREAD)
check_symbol_exists(pread64 unistd.h HAVE_PREAD64)
check_symbol_exists(pwrite unistd.h HAVE_PWRITE)
check_symbol_exists(pwrite64 unistd.h HAVE_PWRITE64)
check_symbol_exists(readlink unistd.h HAVE_READLINK)
check_symbol_exists(usleep unistd.h HAVE_USLEEP)

check_symbol_exists(utime utime.h HAVE_UTIME)

cmake_pop_check_state()

# -- External Dependencies -----------------------------------------------------

if(NOT SQLITE_THREADSAFE STREQUAL "0")
    set(THREADS_PREFER_PTHREAD_FLAG TRUE)
    find_package(Threads REQUIRED)
endif()

if(SQLITE_ENABLE_ICU)
    find_package(ICU REQUIRED COMPONENTS i18n uc)
endif()

if(UNIX)
    find_library(MATH_LIBRARY m)
endif()

# -- Library -------------------------------------------------------------------

add_library(SQLite3)
add_library(SQLite::SQLite3 ALIAS SQLite3)

if(SQLITE_THREADSAFE STREQUAL "0")
    message(STATUS "Building the SQLite3 library in single-threading mode")
    set(SQLITE_THREADMODEL_SUFFIX "-st")
elseif(SQLITE_THREADSAFE STREQUAL "1")
    message(STATUS "Building the SQLite3 library in serialized mode")
    set(SQLITE_THREADMODEL_SUFFIX "")
elseif(SQLITE_THREADSAFE STREQUAL "2")
    message(STATUS "Building the SQLite3 library in multi-threading mode")
    set(SQLITE_THREADMODEL_SUFFIX "-mt")
endif()

set_target_properties(SQLite3 PROPERTIES
    OUTPUT_NAME sqlite3${SQLITE_THREADMODEL_SUFFIX}
    DEBUG_POSTFIX "-dbg"
    VERSION ${SQLite3_VERSION}
    SOVERSION ${SQLite3_VERSION_MAJOR})

target_sources(SQLite3
    PRIVATE sqlite3.c
    INTERFACE FILE_SET HEADERS
        FILES sqlite3.h sqlite3ext.h)

target_include_directories(SQLite3
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>)

if(WIN32 AND BUILD_SHARED_LIBS)
    target_compile_definitions(SQLite3
        PRIVATE "SQLITE_API=__declspec(dllexport)"
        INTERFACE "SQLITE_API=__declspec(dllimport)")
endif()

if(SQLITE_ENABLE_NORMALIZE)
    target_compile_definitions(SQLite3 PUBLIC SQLITE_ENABLE_NORMALIZE)
endif()

if(SQLITE_ENABLE_PREUPDATE_HOOK)
    target_compile_definitions(SQLite3 PUBLIC SQLITE_ENABLE_PREUPDATE_HOOK)
endif()

if(SQLITE_ENABLE_SESSION)
    target_compile_definitions(SQLite3 PUBLIC SQLITE_ENABLE_SESSION)
endif()

if(SQLITE_OMIT_COMPILEOPTION_DIAGS)
    target_compile_definitions(SQLite3 PUBLIC SQLITE_OMIT_COMPILEOPTION_DIAGS)
endif()

if(SQLITE_OMIT_DEPRECATED)
    target_compile_definitions(SQLite3 PUBLIC SQLITE_OMIT_DEPRECATED)
endif()

if(SQLITE_OMIT_DESERIALIZE)
    target_compile_definitions(SQLite3 PUBLIC SQLITE_OMIT_DESERIALIZE)
endif()

if(SQLITE_OMIT_LOAD_EXTENSION)
    target_compile_definitions(SQLite3 PUBLIC SQLITE_OMIT_LOAD_EXTENSION)
endif()

target_compile_definitions(SQLite3
    PRIVATE
        SQLITE_DEFAULT_MEMSTATUS=${SQLITE_DEFAULT_MEMSTATUS}
        SQLITE_DEFAULT_SYNCHRONOUS=${SQLITE_DEFAULT_SYNCHRONOUS}
        SQLITE_DEFAULT_WAL_SYNCHRONOUS=${SQLITE_DEFAULT_WAL_SYNCHRONOUS}
        SQLITE_DQS=${SQLITE_DQS}
        SQLITE_THREADSAFE=${SQLITE_THREADSAFE}
        $<$<NOT:$<STREQUAL:${SQLITE_MAX_ALLOCATION_SIZE},-1>>:SQLITE_MAX_ALLOCATION_SIZE=${SQLITE_MAX_ALLOCATION_SIZE}>
        $<$<NOT:$<STREQUAL:${SQLITE_MAX_EXPR_DEPTH},-1>>:SQLITE_MAX_EXPR_DEPTH=${SQLITE_MAX_EXPR_DEPTH}>
        $<$<BOOL:${SQLITE_ENABLE_API_ARMOR}>:SQLITE_ENABLE_API_ARMOR>
        $<$<BOOL:${SQLITE_ENABLE_COLUMN_METADATA}>:SQLITE_ENABLE_COLUMN_METADATA>
        $<$<BOOL:${SQLITE_ENABLE_DBSTAT_VTAB}>:SQLITE_ENABLE_DBSTAT_VTAB>
        $<$<BOOL:${SQLITE_ENABLE_EXPLAIN_COMMENTS}>:SQLITE_ENABLE_EXPLAIN_COMMENTS>
        $<$<BOOL:${SQLITE_ENABLE_FTS3}>:SQLITE_ENABLE_FTS3>
        $<$<BOOL:${SQLITE_ENABLE_FTS3_PARENTHESIS}>:SQLITE_ENABLE_FTS3_PARENTHESIS>
        $<$<BOOL:${SQLITE_ENABLE_FTS4}>:SQLITE_ENABLE_FTS4>
        $<$<BOOL:${SQLITE_ENABLE_FTS5}>:SQLITE_ENABLE_FTS5>
        $<$<BOOL:${SQLITE_ENABLE_GEOPOLY}>:SQLITE_ENABLE_GEOPOLY>
        $<$<BOOL:${SQLITE_ENABLE_ICU}>:SQLITE_ENABLE_ICU>
        $<$<BOOL:${SQLITE_ENABLE_MATH_FUNCTIONS}>:SQLITE_ENABLE_MATH_FUNCTIONS>
        $<$<BOOL:${SQLITE_ENABLE_RBU}>:SQLITE_ENABLE_RBU>
        $<$<BOOL:${SQLITE_ENABLE_RTREE}>:SQLITE_ENABLE_RTREE>
        $<$<BOOL:${SQLITE_ENABLE_STMTVTAB}>:SQLITE_ENABLE_STMTVTAB>
        $<$<BOOL:${SQLITE_LIKE_DOESNT_MATCH_BLOBS}>:SQLITE_LIKE_DOESNT_MATCH_BLOBS>
        $<$<BOOL:${SQLITE_OMIT_AUTOINIT}>:SQLITE_OMIT_AUTOINIT>
        $<$<BOOL:${SQLITE_OMIT_DECLTYPE}>:SQLITE_OMIT_DECLTYPE>
        $<$<BOOL:${SQLITE_OMIT_JSON}>:SQLITE_OMIT_JSON>
        $<$<BOOL:${SQLITE_OMIT_PROGRESS_CALLBACK}>:SQLITE_OMIT_PROGRESS_CALLBACK>
        $<$<BOOL:${SQLITE_OMIT_SHARED_CACHE}>:SQLITE_OMIT_SHARED_CACHE>
        $<$<BOOL:${SQLITE_STRICT_SUBTYPE}>:SQLITE_STRICT_SUBTYPE=1>
        $<$<BOOL:${SQLITE_USE_ALLOCA}>:SQLITE_USE_ALLOCA>)

target_compile_definitions(SQLite3
    PRIVATE
        $<$<PLATFORM_ID:Windows>:SQLITE_USE_SEH>
        $<$<BOOL:${HAVE_FDATASYNC}>:HAVE_FDATASYNC>
        $<$<BOOL:${HAVE_GETHOSTUUID}>:HAVE_GETHOSTUUID>
        $<$<BOOL:${HAVE_GMTIME_R}>:HAVE_GMTIME_R>
        $<$<BOOL:${HAVE_INT8_T}>:HAVE_INT8_T>
        $<$<BOOL:${HAVE_INT16_T}>:HAVE_INT16_T>
        $<$<BOOL:${HAVE_INTTYPES_H}>:HAVE_INTTYPES_H>
        $<$<BOOL:${HAVE_ISNAN}>:HAVE_ISNAN>
        $<$<BOOL:${HAVE_LOCALTIME_R}>:HAVE_LOCALTIME_R>
        $<$<BOOL:${HAVE_LOCALTIME_S}>:HAVE_LOCALTIME_S>
        HAVE_LOG2=$<BOOL:${HAVE_LOG2}>
        HAVE_LOG10=$<BOOL:${HAVE_LOG10}>
        $<$<BOOL:${HAVE_LSTAT}>:HAVE_LSTAT>
        $<$<BOOL:${HAVE_MALLOC_H}>:HAVE_MALLOC_H>
        $<$<BOOL:${HAVE_MALLOC_USABLE_SIZE}>:HAVE_MALLOC_USABLE_SIZE>
        HAVE_NANOSLEEP=$<BOOL:${HAVE_NANOSLEEP}>
        $<$<BOOL:${HAVE_POSIX_FALLOCATE}>:HAVE_POSIX_FALLOCATE>
        $<$<BOOL:${HAVE_PREAD}>:HAVE_PREAD>
        $<$<BOOL:${HAVE_PREAD64}>:HAVE_PREAD64>
        $<$<BOOL:${HAVE_PWRITE}>:HAVE_PWRITE>
        $<$<BOOL:${HAVE_PWRITE64}>:HAVE_PWRITE64>
        $<$<BOOL:${HAVE_READLINK}>:HAVE_READLINK>
        $<$<BOOL:${HAVE_STDINT_H}>:HAVE_STDINT_H>
        $<$<BOOL:${HAVE_STRCHRNUL}>:HAVE_STRCHRNUL>
        $<$<BOOL:${HAVE_STRERROR_R}>:HAVE_STRERROR_R>
        $<$<BOOL:${HAVE_UINT8_T}>:HAVE_UINT8_T>
        $<$<BOOL:${HAVE_UINT16_T}>:HAVE_UINT16_T>
        $<$<BOOL:${HAVE_UINT32_T}>:HAVE_UINT32_T>
        $<$<BOOL:${HAVE_USLEEP}>:HAVE_USLEEP>
        $<$<BOOL:${HAVE_UTIME}>:HAVE_UTIME>)

if(MATH_LIBRARY)
    target_link_libraries(SQLite3 PUBLIC ${MATH_LIBRARY})
endif()

if(TARGET Threads::Threads)
    target_link_libraries(SQLite3 PUBLIC Threads::Threads)
endif()

if(NOT SQLITE_OMIT_LOAD_EXTENSION)
    target_link_libraries(SQLite3 PUBLIC ${CMAKE_DL_LIBS})
endif()

if(SQLITE_ENABLE_ICU)
    target_link_libraries(SQLite3 PUBLIC ICU::i18n ICU::uc)
endif()

# -- Interactive shell ---------------------------------------------------------

option(SQLITE_BUILD_SHELL "Build the SQLite3 interactive shell." ON)

if(SQLITE_BUILD_SHELL)
    add_executable(Shell)
    add_executable(SQLite::Shell ALIAS Shell)

    set_target_properties(Shell PROPERTIES
        OUTPUT_NAME sqlite3
        INSTALL_RPATH "$ORIGIN/../lib")

    target_sources(Shell
        PRIVATE shell.c)

    target_link_libraries(Shell PRIVATE SQLite::SQLite3)

    if(UNIX)
        find_package(PkgConfig)
        if(PKG_CONFIG_FOUND)
            pkg_search_module(Readline IMPORTED_TARGET readline)
        endif()
    endif()

    find_package(ZLIB)

    if(TARGET PkgConfig::Readline)
        target_compile_definitions(Shell PRIVATE HAVE_READLINE)
        target_link_libraries(Shell PRIVATE PkgConfig::Readline)
    endif()

    if(TARGET ZLIB::ZLIB)
        target_compile_definitions(Shell PRIVATE SQLITE_HAVE_ZLIB)
        target_link_libraries(Shell PRIVATE ZLIB::ZLIB)
    endif()
endif()

# -- Testing -------------------------------------------------------------------

if(BUILD_TESTING AND SQLITE_BUILD_SHELL)
    add_test(NAME SQLiteShellCheckVersion
        COMMAND SQLite::Shell -version)
    set_tests_properties(SQLiteShellCheckVersion PROPERTIES
        PASS_REGULAR_EXPRESSION "${SQLite3_VERSION}")
endif()

# -- Install -------------------------------------------------------------------

include(GNUInstallDirs)

install(TARGETS SQLite3
    EXPORT SQLite3Targets
    FILE_SET HEADERS)

if(SQLITE_BUILD_SHELL)
    install(TARGETS Shell
        EXPORT SQLite3Targets)
endif()

install(EXPORT SQLite3Targets
    NAMESPACE SQLite::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SQLite3)

export(EXPORT SQLite3Targets
    NAMESPACE SQLite::
    FILE ${CMAKE_CURRENT_BINARY_DIR}/SQLite3Targets.cmake)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    SQLite3Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/SQLite3Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SQLite3
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/SQLite3ConfigVersion.cmake
    COMPATIBILITY SameMajorVersion)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SQLite3Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/SQLite3ConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SQLite3)

# -- Packaging -----------------------------------------------------------------

if(PROJECT_IS_TOP_LEVEL)
    include(InstallRequiredSystemLibraries)
    set(CPACK_PACKAGE_VERSION_MAJOR ${SQLite3_VERSION_MAJOR})
    set(CPACK_PACKAGE_VERSION_MINOR ${SQLite3_VERSION_MINOR})
    set(CPACK_PACKAGE_VERSION_PATCH ${SQLite3_VERSION_PATCH})
    include(CPack)
endif()
