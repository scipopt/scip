cmake_minimum_required(VERSION 3.9)

# Set the project name
project(nauty)

# Add the source files
set(SOURCE_FILES
    nauty.c
    nautil.c
    nausparse.c
    schreier.c
    naurng.c
)

# fPIC
if(SHARED)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Add the library target
add_library(nauty STATIC ${SOURCE_FILES})
