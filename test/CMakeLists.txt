add_executable(ParticleTest ParticleTest.cpp )
add_executable(PlanetaryCalculationTest PlanetaryCalculationTest.cpp )
add_executable(LennardJonesPotentialTest LennardJonesPotentialTest.cpp )

target_link_libraries(ParticleTest
        PUBLIC
        GTest::gtest_main
        simulator
)
target_link_libraries(PlanetaryCalculationTest
        PUBLIC
        GTest::gtest_main
        simulator
)
target_link_libraries(LennardJonesPotentialTest
        PUBLIC
        GTest::gtest_main
        simulator
)

target_compile_features(ParticleTest
        PRIVATE
        cxx_std_23
)
target_compile_features(PlanetaryCalculationTest
        PRIVATE
        cxx_std_23
)
target_compile_features(LennardJonesPotentialTest
        PRIVATE
        cxx_std_23
)

target_include_directories(ParticleTest
        PRIVATE
        ${gtest_SOURCE_DIR}/include ${gtest_BINARY_DIR}
        ${SIMULATOR_INCLUDES}
)
target_include_directories(PlanetaryCalculationTest
        PRIVATE
        ${gtest_SOURCE_DIR}/include ${gtest_BINARY_DIR}
        ${SIMULATOR_INCLUDES}
)
target_include_directories(LennardJonesPotentialTest
        PRIVATE
        ${gtest_SOURCE_DIR}/include ${gtest_BINARY_DIR}
        ${SIMULATOR_INCLUDES}
)

gtest_discover_tests(ParticleTest)
gtest_discover_tests(PlanetaryCalculationTest)
gtest_discover_tests(LennardJonesPotentialTest)
