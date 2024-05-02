find_package(Doxygen
        REQUIRED dot)

if (DOXYGEN_FOUND)
    doxygen_add_docs(doc_doxygen
            ${PROJECT_SOURCE_DIR}
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            COMMENT "Building Documentation"
            CONFIG_FILE ${PROJECT_SOURCE_DIR}/Doxyfile
    )
endif (DOXYGEN_FOUND)
