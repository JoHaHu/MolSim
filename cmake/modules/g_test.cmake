include(FetchContent)
FetchContent_Declare(
        g_test
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
FetchContent_MakeAvailable(g_test)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# Discover tests
include(GoogleTest)

# for Googletest Suite
enable_testing()
add_subdirectory(test)