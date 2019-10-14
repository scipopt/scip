#
# @author Gregor Hendel
#
# this directory adds a doc target
#
#
#

#
# a couple of programs are required on the path
#
find_program(PYTHON python)
find_program(PHP php)
find_program(BASH bash)
find_program(DOXYGEN doxygen)


if(NOT (PYTHON AND PHP AND BASH AND DOXYGEN))
    message(WARNING "Could not find all of doxygen, python, php, and bash in system path, which are required for generation of documentation")
else()

    # change this variable to 'php' the web page

    set(DOC_HTML_FILE_EXTENSION html)

    # Shell Tutorial Is Built In 2 Steps
    #
    # 1a. create the raw shell tutorial using commands for SCIP
    #
    add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/inc/shelltutorial/shelltutorialraw.tmp
        COMMAND $<TARGET_FILE:scip> < doc/inc/shelltutorial/commands | tee doc/inc/shelltutorial/shelltutorialraw.tmp
        COMMAND ${BASH} -c "rm stein27.lp stein27.sol settingsfile.set"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../
    )

    #
    # 1b. create the shelltutorial that contains annotations for doxygen
    #
    add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/inc/shelltutorial/shelltutorialannotated.tmp
        COMMAND ${PYTHON} inc/shelltutorial/insertsnippetstutorial.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/inc/shelltutorial/shelltutorialraw.tmp
    )

    #
    # 2. create the FAQ
    #
    add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/inc/faq.inc
        COMMAND ${PYTHON} parser.py --linkext ${DOC_HTML_FILE_EXTENSION}  && ${PHP} localfaq.php > faq.inc
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/inc/faq
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR/inc/faq/faqtext.txt}
        )

    #
    # 3. parameter file generation
    #
    add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/inc/parameters.set
        COMMAND $<TARGET_FILE:scip> -c "set default set save inc/parameters.set quit"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    #
    # 4. front page tutorial tutorial
    #
    # bin/scip -c "read doc/inc/simpleinstance/simple.lp optimize quit" > doc/inc/simpleinstance/output.log
    add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/inc/simpleinstance/output.log
        COMMAND $<TARGET_FILE:scip> -c "read inc/simpleinstance/simple.lp optimize quit" -l inc/simpleinstance/output.log
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )


    add_custom_target(doc
        ${DOXYGEN} scip.dxy
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/inc/shelltutorial/shelltutorialannotated.tmp
                ${CMAKE_CURRENT_SOURCE_DIR}/inc/faq.inc
                ${CMAKE_CURRENT_SOURCE_DIR}/inc/parameters.set
                ${CMAKE_CURRENT_SOURCE_DIR}/inc/simpleinstance/output.log
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generating doxygen documentation in ${CMAKE_CURRENT_SOURCE_DIR}/html/ ..."
                )
endif()
