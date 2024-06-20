cmake_minimum_required(VERSION 3.24.1)

# define project name, version
project(PSEMolDyn_GroupD VERSION 0.0.1 LANGUAGES CXX)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel;Profile")

cmake_policy(SET CMP0135 NEW)

# set Release as the default build type if it is not yet set.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif (NOT CMAKE_BUILD_TYPE)


##### LIBRARY #####

# collect all cpp files
file(GLOB_RECURSE MY_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/*.cpp"
        # we use cxx extension to exclude generated files from clang-tidy
        "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/*.cxx"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/*.h"
        # we use hxx extension to exclude generated files from clang-tidy
        "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/*.hxx"
)

add_library(simulator
        STATIC
        ${MY_SRC}
)

target_include_directories(simulator
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/lib"
)

target_link_libraries(simulator
        PUBLIC
        #        asan
        spdlog::spdlog
        XercesC::XercesC
        range-v3::range-v3

)

target_compile_features(simulator
        PRIVATE
        cxx_std_20
)


target_compile_options(simulator
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall -march=native>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -march=native -Xhost -wd383,981,1418,1572,2259>
        $<$<CXX_COMPILER_ID:Clang>:-Wall -g -march=x86-64-v3  #[[-march=native]] >
)
get_target_property(SIMULATOR_INCLUDES simulator INCLUDE_DIRECTORIES)


####  EXECUTABLE #####


# create make target
add_executable(MolSim src/bin/MolSim.cpp)


# set cxx standard. You may raise this if you want.
target_compile_features(MolSim
        PRIVATE
        cxx_std_20
)

target_include_directories(MolSim
        PUBLIC
        ${SIMULATOR_INCLUDES}
)
target_link_libraries(MolSim
        PRIVATE
        simulator
)

include(CheckIPOSupported)
check_ipo_supported(RESULT result)
if (result)
    set_property(TARGET MolSim PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET simulator PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()


target_compile_options(MolSim
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall -march=native>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -march=native -Xhost -xCORE-AVX512 -wd383,981,1418,1572,2259>
        $<$<CXX_COMPILER_ID:Clang>:-Wall -march=x86-64-v3 #[[-march=native]]   #[[ -fsanitize=address,leak]] >
)


# Configure Debug and Release specific options
if (CMAKE_BUILD_TYPE STREQUAL "Debug")

    target_compile_options(simulator
            PRIVATE
            $<$<CXX_COMPILER_ID:Clang>: -O0 -g -fno-omit-frame-pointer -fdebug-info-for-profiling -fsanitize=address >
    )
    target_compile_options(MolSim
            PRIVATE
            $<$<CXX_COMPILER_ID:Clang>: -O0 -g -fno-omit-frame-pointer -fdebug-info-for-profiling -fsanitize=address >
    )
    target_compile_definitions(MolSim PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)
    target_compile_definitions(simulator PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)
    target_compile_definitions(MolSim PUBLIC DEBUG=true)
    target_compile_definitions(simulator PUBLIC DEBUG=true)
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(simulator
            PRIVATE
            $<$<CXX_COMPILER_ID:Intel>: -O3 -ipo -fast -mprefer-vector-width=512 -fvectorize -fprof-use=../profiles/default.profdata>
            $<$<CXX_COMPILER_ID:Clang>: -O3 -mprefer-vector-width=512 #[[-fsanitize=address]] -fprofile-instr-use=../profiles/default.profdata -fvectorize -Rpass=vectorize >
            $<$<CXX_COMPILER_ID:GNU>: -O3 -ftree-vectorize -mprefer-vector-width=512>
    )
    target_compile_options(MolSim
            PRIVATE
            $<$<CXX_COMPILER_ID:Intel>: -O3 -ipo -fast -mprefer-vector-width=512 -fvectorize -fprof-use=../profiles/default.profdata>
            $<$<CXX_COMPILER_ID:Clang>: -O3 -mprefer-vector-width=512  #[[-fsanitize=address]]  -fprofile-instr-use=../profiles/default.profdata -fvectorize -Rpass=vectorize >
            $<$<CXX_COMPILER_ID:GNU>: -O3 -ftree-vectorize -mprefer-vector-width=512>
    )
    target_compile_definitions(MolSim PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_INFO)
    target_compile_definitions(simulator PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_INFO)

    target_link_options(simulator
            PRIVATE
            $<$<CXX_COMPILER_ID:Clang>: -Wl -fprofile-use=../profiles/default.profdata>

    )
    target_link_options(MolSim
            PRIVATE
            $<$<CXX_COMPILER_ID:Clang>: -Wl -fprofile-use=../profiles/default.profdata >

    )


elseif (CMAKE_BUILD_TYPE STREQUAL "Profile")
    target_compile_options(simulator
            PRIVATE
            $<$<CXX_COMPILER_ID:Clang>: -O3 -Rpass=loop-vectorize -mprefer-vector-width=512 -fcs-profile-generate=../profiles -fvectorize -Rpass-missed=loop-vectorize -fno-omit-frame-pointer -fdebug-info-for-profiling>
            $<$<CXX_COMPILER_ID:Intel>: -O3 -ipo -fast Rpass=loop-vectorize -mprefer-vector-width=512 -fprof-gen=../profiles -fvectorize -Rpass-missed=loop-vectorize -fno-omit-frame-pointer -fdebug-info-for-profiling>
            $<$<CXX_COMPILER_ID:GNU>: -O3 -ftree-vectorize -mprefer-vector-width=512>
    )
    target_compile_options(MolSim
            PRIVATE
            $<$<CXX_COMPILER_ID:Clang>: -O3 -Rpass=loop-vectorize -mprefer-vector-width=512 -fcs-profile-generate=../profiles -fvectorize -Rpass-missed=loop-vectorize -fno-omit-frame-pointer -fdebug-info-for-profiling>
            $<$<CXX_COMPILER_ID:Intel>: -O3 -ipo -fast -Rpass=loop-vectorize -mprefer-vector-width=512 -fprof-gene=../profiles -fvectorize -Rpass-missed=loop-vectorize -fno-omit-frame-pointer -fdebug-info-for-profiling>
            $<$<CXX_COMPILER_ID:GNU>: -O3 -ftree-vectorize -mprefer-vector-width=512>
    )
    target_compile_definitions(MolSim PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_INFO)
    target_compile_definitions(simulator PUBLIC SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_INFO)

    target_link_options(simulator
            PRIVATE
            $<$<CXX_COMPILER_ID:Clang>: -Wl -fcs-profile-generate=../profiles >
            $<$<CXX_COMPILER_ID:Intel>: -Wl -fprof-gen=../profiles >

    )
    target_link_options(MolSim
            PRIVATE
            $<$<CXX_COMPILER_ID:Clang>: -Wl -fcs-profile-generate=../profiles >
            $<$<CXX_COMPILER_ID:Intel>: -Wl -fprof-gen=../profiles >

    )
endif ()


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
include(doxygen)
include(spdlog)
include(g_test)
