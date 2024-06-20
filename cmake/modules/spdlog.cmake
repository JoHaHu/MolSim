#### SPDLOG ###
include(FetchContent)
fetchcontent_declare(spdlog
        URL https://github.com/gabime/spdlog/archive/refs/tags/v1.14.1.tar.gz
        OVERRIDE_FIND_PACKAGE)
fetchcontent_makeavailable(spdlog)
find_package(spdlog REQUIRED)
### Xerces
find_package(XercesC 3.2.1 REQUIRED)


### ranges-v3
include(FetchContent)
fetchcontent_declare(ranges-v3
        URL https://github.com/ericniebler/range-v3/archive/refs/tags/0.12.0.tar.gz
        OVERRIDE_FIND_PACKAGE)
fetchcontent_makeavailable(ranges-v3)
find_package(ranges-v3 REQUIRED)
