set(PROJECT_NAME PolySCIP)
project(${PROJECT_NAME})

find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})

# set variables
set(EXECUTABLE_NAME "polyscip")
set(POLYSCIP_VERSION_MAJOR 2)
set(POLYSCIP_VERSION_MINOR 0)
set(POLYSCIP_MAX_NUMBER_OBJS "30" CACHE STRING "Maximum number of objectives for PolySCIP instances")
set(SRC src)
set(DOC doc)
set(CMAKE cmake)
set(RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/${CMAKE})
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic -Wno-unused-but-set-variable -Wno-unused-variable")


# configure a header file to pass CMake settings to source code
configure_file("${PROJECT_SOURCE_DIR}/${CMAKE}/PolySCIPConfig.h.in" "${PROJECT_SOURCE_DIR}/${SRC}/PolySCIPConfig.h" ESCAPE_QUOTES)

set(HEADERS_FILES
        ${SRC}/cmd_line_args.h
        ${SRC}/global_functions.h
        ${SRC}/polyscip.h
        ${SRC}/polyscip_types.h
        ${SRC}/double_description_method.h
        ${SRC}/prob_data_objectives.h
        ${SRC}/ReaderMOP.h
        ${SRC}/weight_space_facet.h
        ${SRC}/weight_space_polyhedron.h
        ${SRC}/weight_space_vertex.h
        )

set(SOURCE_FILES
        ${SRC}/cmd_line_args.cpp
        ${SRC}/main.cpp
        ${SRC}/polyscip.cpp
        ${SRC}/double_description_method.cpp
        ${SRC}/prob_data_objectives.cpp
        ${SRC}/ReaderMOP.cpp
        ${SRC}/weight_space_facet.cpp
        ${SRC}/weight_space_polyhedron.cpp
        ${SRC}/weight_space_vertex.cpp
        )

include_directories(${PROJECT_SOURCE_DIR}/${SRC})

add_executable(${EXECUTABLE_NAME} ${SOURCE_FILES})
target_link_libraries(polyscip ${SCIP_LIBRARIES})      
add_dependencies(applications ${EXECUTABLE_NAME})



