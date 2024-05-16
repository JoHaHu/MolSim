add_executable(Test Test.cpp ParticleTest.cpp PlanetaryCalculationTest.cpp LennardJonesPotentialTest.cpp)

target_link_libraries(Test
        PUBLIC
        GTest::gtest_main
        simulator
)

target_compile_features(Test
        PRIVATE
        cxx_std_23
)

target_include_directories(Test
        PRIVATE
        ${gtest_SOURCE_DIR}/include ${gtest_BINARY_DIR}
        ${SIMULATOR_INCLUDES}
)

gtest_discover_tests(Test)
