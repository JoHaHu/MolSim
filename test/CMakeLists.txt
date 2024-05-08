cmake_minimum_required(VERSION 3.22.1)
#project(Tests)

# Setup Google Test
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
FetchContent_MakeAvailable(googletest)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Define the test executable
add_executable(ParticleTest ParticleTest.cpp)
target_link_libraries(ParticleTest GTest::gtest_main)


# Specify the target directory (should be automatically done by FetchContent but just for good measure)
target_include_directories(ParticleTest
        PRIVATE
        ${gtest_SOURCE_DIR}/include ${gtest_BINARY_DIR})

# Discover tests
include(GoogleTest)
gtest_discover_tests(ParticleTest)