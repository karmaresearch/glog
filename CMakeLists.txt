cmake_minimum_required (VERSION 2.8.12)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
if (NOT "$ENV{CC}" STREQUAL "")
    set(CMAKE_C_COMPILER "$ENV{CC}")
endif()
if (NOT "$ENV{CXX}" STREQUAL "")
    set(CMAKE_CXX_COMPILER "$ENV{CXX}")
endif()

project (vlog)
include(ExternalProject)
set(CMAKE_MACOSX_RPATH 1)

#Set to release mode if nothing is specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

#Add macPorts
# Detect if the "port" command is valid on this system; if so, return full path
EXECUTE_PROCESS(COMMAND which port RESULT_VARIABLE DETECT_MACPORTS OUTPUT_VARIABLE MACPORTS_PREFIX ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

IF (${DETECT_MACPORTS} EQUAL 0)
    # "/opt/local/bin/port" doesn't have libs, so we get the parent directory
    GET_FILENAME_COMPONENT(MACPORTS_PREFIX ${MACPORTS_PREFIX} DIRECTORY)

    # "/opt/local/bin" doesn't have libs, so we get the parent directory
    GET_FILENAME_COMPONENT(MACPORTS_PREFIX ${MACPORTS_PREFIX} DIRECTORY)

    # "/opt/local" is where MacPorts lives, add `/lib` suffix and link
    LINK_DIRECTORIES(${MACPORTS_PREFIX}/lib)
    INCLUDE_DIRECTORIES(${MACPORTS_PREFIX}/include)
ENDIF()

#Add homebrew
EXECUTE_PROCESS(COMMAND brew --prefix RESULT_VARIABLE exitstatus OUTPUT_VARIABLE BREW_PREFIX ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

IF (${exitstatus} EQUAL 0)
    LINK_DIRECTORIES(${BREW_PREFIX}/lib)
    INCLUDE_DIRECTORIES(${BREW_PREFIX}/include)
ENDIF()

#add files to compile
file(GLOB vlog_SRC
    "src/vlog/*.cpp"
    "src/vlog/common/*.cpp"
    "src/vlog/backward/*.cpp"
    "src/vlog/forward/*.cpp"
    "src/vlog/magic/*.cpp"
    "src/vlog/incremental/*.cpp"
    "src/vlog/embeddings/*.cpp"
    "src/vlog/inmemory/*.cpp"
    "src/vlog/text/*.cpp"
    "src/vlog/trident/*.cpp"
    "src/vlog/cycles/*.cpp"
    "src/vlog/ml/*.cpp"
    "src/vlog/trigger/*.cpp"
    "src/launcher/vloglayer.cpp"
    "src/launcher/vlogscan.cpp"
    "src/vlog/prediction/*.cpp"
    "src/vlog/clique/*.cpp"
    "src/vlog/util/*.cpp"
    "src/glog/*.cpp"
    )

#ZLIB
find_package(ZLIB REQUIRED)

#Trident
IF(TRIDENT_LIB)
    message("Using KOGNAC provided installation: ${KOGNAC_LIB} ${KOGNAC_INC}")
    include_directories(${KOGNAC_INC})
    LINK_DIRECTORIES("${KOGNAC_LIB}")
    add_library(kognac-core SHARED IMPORTED)
    set_property(TARGET kognac-core PROPERTY IMPORTED_LOCATION ${KOGNAC_LIB}/${CMAKE_SHARED_LIBRARY_PREFIX}kognac-core${CMAKE_SHARED_LIBRARY_SUFFIX})

    message("Using Trident provided installation: ${TRIDENT_LIB} ${TRIDENT_INC}")
    include_directories(${TRIDENT_INC})
    include_directories(${TRIDENT_INC}/layers)
    include_directories(${TRIDENT_INC}/../rdf3x/include)
    add_library(trident SHARED IMPORTED)
    set_property(TARGET trident PROPERTY IMPORTED_LOCATION ${TRIDENT_LIB}/libtrident-core${CMAKE_SHARED_LIBRARY_SUFFIX})
    add_library(trident-sparql SHARED IMPORTED)
    set_property(TARGET trident-sparql PROPERTY IMPORTED_LOCATION ${TRIDENT_LIB}/libtrident-sparql${CMAKE_SHARED_LIBRARY_SUFFIX})
    LINK_DIRECTORIES("${TRIDENT_LIB}/")
ELSE()
    ExternalProject_Add(git-trident
        DOWNLOAD_COMMAND git clone https://github.com/karmaresearch/trident.git
        DOWNLOAD_DIR .
        SOURCE_DIR trident
        BUILD_IN_SOURCE 1
        CMAKE_ARGS "-DSPARQL=1" "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
        INSTALL_COMMAND ""
        )
    ExternalProject_Get_Property(git-trident SOURCE_DIR)
    include_directories(${SOURCE_DIR}/include)
    include_directories(${SOURCE_DIR}/include/layers)
    include_directories(${SOURCE_DIR}/rdf3x/include)
    # TRIDENT
    add_library(trident SHARED IMPORTED)
    add_dependencies(trident git-trident)
    set_property(TARGET trident PROPERTY IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/trident/libtrident-core${CMAKE_SHARED_LIBRARY_SUFFIX}")
    LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/trident/")
    # TRIDENT SPARQL
    add_library(trident-sparql SHARED IMPORTED)
    add_dependencies(trident-sparql git-trident)
    set_property(TARGET trident-sparql PROPERTY IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/trident/libtrident-sparql${CMAKE_SHARED_LIBRARY_SUFFIX}")
    # KOGNAC
    set(KOGNAC_LIB ${CMAKE_BINARY_DIR}/trident/kognac)
    add_library(kognac-core SHARED IMPORTED)
    add_dependencies(kognac-core git-trident)
    include_directories(${SOURCE_DIR}/kognac/include)
    LINK_DIRECTORIES("${KOGNAC_LIB}")
    set_property(TARGET kognac-core PROPERTY IMPORTED_LOCATION ${KOGNAC_LIB}/${CMAKE_SHARED_LIBRARY_PREFIX}kognac-core${CMAKE_SHARED_LIBRARY_SUFFIX})
ENDIF()

#SPARSEHASH
find_path(sparsehash sparsehash/sparse_hash_map)
IF (${sparsehash} STREQUAL "sparsehash-NOTFOUND")
    message("Could not find sparsehash. I'm going to use the version from kognac ...")
    include_directories(${KOGNAC_LIB}/external/sparsehash/src/)
ENDIF()

message("I'm going to use LZ4 the version from kognac")
include_directories(${KOGNAC_LIB}/external/lz4/lib/)
LINK_DIRECTORIES(${KOGNAC_LIB}/external/lz4/lib)
add_library(lz4 STATIC IMPORTED)
set_property(TARGET lz4 PROPERTY IMPORTED_LOCATION ${KOGNAC_LIB}/external/lz4/lib/${CMAKE_STATIC_LIBRARY_PREFIX}lz4${CMAKE_STATIC_LIBRARY_SUFFIX})

#Webserver support
IF(WEBINTERFACE)
    file(GLOB vlog_webSRC "src/vlog/web/*.cpp")
    set(vlog_SRC ${vlog_webSRC} ${vlog_SRC})
    set(COMPILE_FLAGS "${COMPILE_FLAGS} -DWEBINTERFACE=1")
ENDIF()

#Create both a library and the executable program
add_library(vlog-core SHARED ${vlog_SRC})
add_executable(vlog src/launcher/main.cpp)

#PTHREADS
find_package(Threads REQUIRED)
if (THREADS_HAVE_PTHREAD_ARG)
    target_compile_options(PUBLIC vlog-core "-pthread")
endif()
if (CMAKE_THREAD_LIBS_INIT)
    target_link_libraries(vlog-core "${CMAKE_THREAD_LIBS_INIT}")
endif()

#Set compiler options
if((${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU"))
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        set(COMPILE_FLAGS "${COMPILE_FLAGS} -c -fPIC -MD -std=c++11 -DCOMPRPROOFS=1 -DPRUNING_QSQR=1 ${CLANG_ONLY_COMPILE_FLAGS}")
        if(LINUX)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lstdc++")
        endif()
    else()
        set(COMPILE_FLAGS "${COMPILE_FLAGS} -c -fPIC -MD -std=c++11 -DCOMPRPROOFS=1 -DPRUNING_QSQR=1 ${GNUCXX_ONLY_COMPILE_FLAGS}")
    endif()
endif()

set_target_properties(vlog-core PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
set_target_properties(vlog PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}" OUTPUT_NAME "glog")

#standard include
include_directories(include/)
TARGET_LINK_LIBRARIES(vlog-core trident trident-sparql ${ZLIB_LIBRARIES} kognac-core lz4)
TARGET_LINK_LIBRARIES(vlog vlog-core)
