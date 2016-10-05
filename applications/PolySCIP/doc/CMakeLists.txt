# add a target to generate API documentation with Doxygen
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(BUILD_DOCUMENTATION)
    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()

    set(doxyfile polyscip.dxy)
    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_SOURCE_DIR}/${DOC}/${doxyfile}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/${DOC}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()
