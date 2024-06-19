#### SPDLOG ###
include(FetchContent)
fetchcontent_declare(spdlog
        URL https://github.com/gabime/spdlog/archive/refs/tags/v1.14.1.tar.gz
        OVERRIDE_FIND_PACKAGE)
fetchcontent_makeavailable(spdlog)
find_package(spdlog REQUIRED)

find_package(XercesC 3.2.1 REQUIRED)
