cmake_minimum_required(VERSION 3.9)

# Set the project name
project(nauty)

# Add the source files
set(SOURCE_FILES
    gtools.c
    naugroup.c
    naugstrings.c
    naurng.c
    nausparse.c
    nautaux.c
    nautinv.c
    nauty.c
    nautycliquer.c
    planarity.c
    schreier.c
    traces.c
)

# Add the library target
add_library(nauty STATIC ${SOURCE_FILES})
