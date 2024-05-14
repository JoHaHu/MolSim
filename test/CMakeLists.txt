


add_executable(ParticleTest ParticleTest.cpp ${MY_SRC})
target_link_libraries(ParticleTest
        PUBLIC
        GTest::gtest_main
        xerces-c
)

target_compile_features(ParticleTest
        PRIVATE
        cxx_std_23
)

target_include_directories(ParticleTest
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/../libs/libxsd
        PRIVATE
        ${gtest_SOURCE_DIR}/include ${gtest_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../src
)


gtest_discover_tests(ParticleTest)