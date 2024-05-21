#### Ranges library ###
include(FetchContent)
fetchcontent_declare(ranges
        URL https://github.com/ericniebler/range-v3/archive/refs/tags/0.12.0.tar.gz
        OVERRIDE_FIND_PACKAGE)
fetchcontent_makeavailable(ranges)
find_package(ranges REQUIRED)