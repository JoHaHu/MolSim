file(GLOB_RECURSE MY_TESTS
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

add_executable(Test Test.cpp ${MY_TESTS})

include(CheckIPOSupported)
check_ipo_supported(RESULT result)
if (result)
    set_property(TARGET Test PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

target_link_libraries(Test
        PUBLIC
        GTest::gtest_main
        simulator
)

target_compile_features(Test
        PRIVATE
        cxx_std_23
)

target_compile_options(Test
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall -march=native>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
        $<$<CXX_COMPILER_ID:Clang>:-Wall -march=native>
)


target_link_options(
        Test
        PRIVATE
        $<$<CXX_COMPILER_ID:Clang>: -fuse-ld=lld>
)


target_include_directories(Test
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/test/
        ${gtest_SOURCE_DIR}/include
        ${gtest_BINARY_DIR}
        ${SIMULATOR_INCLUDES}
)
gtest_discover_tests(Test)
