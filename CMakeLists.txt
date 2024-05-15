cmake_minimum_required(VERSION 3.22.1)

# define project name, version
project(PSEMolDyn_GroupD VERSION 0.0.1 LANGUAGES CXX)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")


cmake_policy(SET CMP0135 NEW)

# set Release as the default build type if it is not yet set.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif (NOT CMAKE_BUILD_TYPE)


# collect all cpp files
file(GLOB_RECURSE MY_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/*.cpp"
        # we use cxx extension to exclude generated files from clang-tidy
        "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/*.cxx"
        # header don't need to be included but this might be necessary for some IDEs
        "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/*.h"
        # we use hxx extension to exclude generated files from clang-tidy
        "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/*.hxx"
)


# create make target
add_executable(MolSim src/bin/MolSim.cpp ${MY_SRC})


# set cxx standard. You may raise this if you want.
target_compile_features(MolSim
        PRIVATE
        cxx_std_23
)

target_include_directories(MolSim
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/spdlog
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)


target_link_libraries(MolSim
        # stuff that is used in headers and source files
        PUBLIC
        xerces-c
        spdlog::spdlog
)

include(CheckIPOSupported)
check_ipo_supported(RESULT result)
if (result)
    set_property(TARGET MolSim PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
        $<$<CXX_COMPILER_ID:Clang>:-Wall>
)

# Spdlog configuration
add_subdirectory(libs/spdlog)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs/spdlog/include)

# Configure Debug and Release specific options
if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    target_compile_options(MolSim
            PRIVATE
            $<$<CXX_COMPILER_ID:Clang>: -g -O0 -fno-omit-frame-pointer>
    )
elseif (CMAKE_BUILD_TYPE STREQUAL "RELEASE")
    target_compile_options(MolSim
            PRIVATE
            $<$<CXX_COMPILER_ID:Clang>: -O3 -Rpass=loop-vectorize -Rpass-missed=loop-vectorize >
    )
endif ()


# Setup Google Test
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
FetchContent_MakeAvailable(googletest)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# Discover tests
include(GoogleTest)


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

include(doxygen)

# for Googletest Suite
enable_testing()

add_subdirectory(test)

