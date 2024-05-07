cmake_minimum_required(VERSION 3.0)
project(velocityTest)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

# Aggiungi gli header alla lista degli include
include_directories(src/GNSSlibrary)

# Aggiungi i file sorgente della libreria GNSSlibrary
set(LIB_SOURCES
    src/GNSSlibrary/utility.cpp
    src/GNSSlibrary/atmo.cpp
    src/GNSSlibrary/broadcastEph.cpp
    src/GNSSlibrary/coresky.cpp
    src/GNSSlibrary/rinex3Nav.cpp
    
)

# Crea una libreria con i file sorgente della libreria GNSSlibrary
add_library(GNSSlibrary ${LIB_SOURCES})


add_executable(velocityTest
    scripts/velocityTest.cpp
    src/GNSSlibrary/utility.cpp
    src/GNSSlibrary/atmo.cpp
    src/GNSSlibrary/coresky.cpp
)

add_executable(rinex
    src/GNSSlibrary/utility.cpp
    src/GNSSlibrary/rinex3Nav.cpp
    src/GNSSlibrary/broadcastEph.cpp
)

add_executable(preciseEph
    scripts/core.cpp
    src/GNSSlibrary/utility.cpp
    src/GNSSlibrary/coresky.cpp
    src/GNSSlibrary/utility.cpp
)

add_executable(atmoMain
    scripts/main.cpp
    src/GNSSlibrary/utility.cpp
    src/GNSSlibrary/atmo.cpp
)

# Link della libreria e delle librerie necessarie
target_link_libraries(velocityTest GNSSlibrary armadillo)
target_link_libraries(rinex GNSSlibrary armadillo)
target_link_libraries(preciseEph GNSSlibrary armadillo)
target_link_libraries(atmoMain GNSSlibrary armadillo)
