add_executable(ParticleTest ParticleTest.cpp ${MY_SRC})
add_executable(PlanetaryCalculationTest PlanetaryCalculationTest.cpp ${MY_SRC})
add_executable(LennardJonesPotentialTest LennardJonesPotentialTest.cpp ${MY_SRC})

target_link_libraries(ParticleTest
        PUBLIC
        GTest::gtest_main
        xerces-c
)
target_link_libraries(PlanetaryCalculationTest
        PUBLIC
        GTest::gtest_main
        xerces-c
)
target_link_libraries(LennardJonesPotentialTest
        PUBLIC
        GTest::gtest_main
        xerces-c
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
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/../libs/libxsd
        PRIVATE
        ${gtest_SOURCE_DIR}/include ${gtest_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../src
)
target_include_directories(PlanetaryCalculationTest
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/../libs/libxsd
        PRIVATE
        ${gtest_SOURCE_DIR}/include ${gtest_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../src
)
target_include_directories(LennardJonesPotentialTest
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/../libs/libxsd
        PRIVATE
        ${gtest_SOURCE_DIR}/include ${gtest_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

gtest_discover_tests(ParticleTest)
gtest_discover_tests(PlanetaryCalculationTest)
gtest_discover_tests(LennardJonesPotentialTest)
