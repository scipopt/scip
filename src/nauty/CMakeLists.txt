cmake_minimum_required(VERSION 3.9)

# Set the project name
project(nauty)

# Add the source files
set(SOURCE_FILES
    gtools.c
    naugroup.c
    naurng.c
    nausparse.c
    nauty.c
    nautycliquer.c
    schreier.c
    nautil.c
)

# Add the library target
add_library(nauty STATIC ${SOURCE_FILES})
