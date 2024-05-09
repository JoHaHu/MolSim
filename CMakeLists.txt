cmake_minimum_required(VERSION 3.22.1)

# Define project name, version
project(PSEMolDyn_GroupD VERSION 0.0.1 LANGUAGES CXX)

# Let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# Set Release as the default build type if it is not yet set.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif ()

# Collect all cpp files
file(GLOB_RECURSE MY_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

#TODO Do not rely on system wide installed libraries and do not create extra installation steps for the user to perform manually (SPDLOG AND BOOST) submodule?

# Create make target
add_executable(MolSim ${MY_SRC})

# Set CXX standard
target_compile_features(MolSim PRIVATE cxx_std_23)

# Include directories
target_include_directories(MolSim
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/boost
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/spdlog
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Link libraries
target_link_libraries(MolSim
        PUBLIC
        xerces-c
        ${Boost_LIBRARIES}
        spdlog::spdlog
)

# Compiler options
target_compile_options(MolSim
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall>
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
endif ()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

include(doxygen)

# for Googletest Suite
enable_testing()
add_subdirectory(test)

