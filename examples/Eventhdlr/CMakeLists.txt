project(Eventhdlr)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

add_executable(eventhdlr
   src/cmain.c
   src/event_bestsol.c
   src/event_boundwriting.c)

target_link_libraries(eventhdlr libscip)
add_dependencies(examples eventhdlr)
