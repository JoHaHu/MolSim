cmake_minimum_required(VERSION 3.22.1)

# define project name, version
project(PSEMolDyn_GroupD VERSION 0.0.1 LANGUAGES CXX)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif (NOT CMAKE_BUILD_TYPE)


# collect all cpp files
file(GLOB_RECURSE MY_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
        # we use cxx extension to exclude generated files from clang-tidy
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx"
        # header don't need to be included but this might be necessary for some IDEs
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
        # we use hxx extension to exclude generated files from clang-tidy
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hxx"
)


# create make target
add_executable(MolSim ${MY_SRC})


# set cxx standard. You may raise this if you want.
target_compile_features(MolSim
        PRIVATE
        cxx_std_23
)

target_include_directories(MolSim
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/boost
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/spdlog
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)


target_link_libraries(MolSim
        # stuff that is used in headers and source files
        PUBLIC
        xerces-c
        ${Boost_LIBRARIES}
        spdlog::spdlog
)

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
        $<$<CXX_COMPILER_ID:Clang>:-Wall>
)

# Boost configuration
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
if (Boost_FOUND)
    target_include_directories(MolSim PRIVATE ${Boost_INCLUDE_DIRS})
    target_link_libraries(MolSim PRIVATE ${Boost_LIBRARIES})
endif ()

# Spdlog configuration
find_package(spdlog REQUIRED)
if (spdlog_FOUND)
    target_link_libraries(MolSim PRIVATE spdlog::spdlog)
endif ()

# Configure Debug and Release specific options
if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    target_compile_options(MolSim
            PRIVATE
            $<$<CXX_COMPILER_ID:Clang>: -g -O1 -fno-omit-frame-pointer>
    )
elseif (CMAKE_BUILD_TYPE STREQUAL "RELEASE")
    target_compile_options(MolSim
            PRIVATE
            $<$<CXX_COMPILER_ID:Clang>: -O3 -Rpass=loop-vectorize -Rpass-missed=loop-vectorize >
    )
endif (CMAKE_BUILD_TYPE EQUAL "DEBUG")


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
include(doxygen)

# for Googletest Suite
enable_testing()
add_subdirectory(test)

