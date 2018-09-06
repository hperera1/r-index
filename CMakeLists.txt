cmake_minimum_required(VERSION 2.6)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/.. )

project (s-rlbwt)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/internal)
include_directories(${PROJECT_SOURCE_DIR}/sdsl-lite/include/)
include_directories(${PROJECT_SOURCE_DIR}/sdsl-lite/external/libdivsufsort/include)
# include_directories(~/include) #SDSL headers are here

LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/sdsl-lite/lib)
LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/sdsl-lite/external/libdivsufsort/lib)
# LINK_DIRECTORIES(~/lib) #SDSL lib are here

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS "--std=c++17")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -march=ivybridge")
set(CMAKE_CXX_FLAGS_RELEASE "-g -Ofast -fstrict-aliasing -march=ivybridge -DNDEBUG")
# set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fstrict-aliasing -march=ivyvridge -DNDEBUG -DLF_STATS")  # Uncomment this line in order to enable LF_STATS.
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=ivybridge")


#add_custom_command(TARGET sdsl-lite/lib/libsdsl.a
#        PRE_BUILD
#        COMMAND "cd sdsl-lite && cmake . && make && cd .."
#        )
# I'm trying to force it to build sdsl-lite, but having troubles.
# Will fix later.

add_executable(ri-build ri-build.cpp)
TARGET_LINK_LIBRARIES(ri-build sdsl)
TARGET_LINK_LIBRARIES(ri-build divsufsort)
TARGET_LINK_LIBRARIES(ri-build divsufsort64)

add_executable(ri-buildfasta ri-buildfasta.cpp)
TARGET_LINK_LIBRARIES(ri-buildfasta sdsl)
TARGET_LINK_LIBRARIES(ri-buildfasta divsufsort)
TARGET_LINK_LIBRARIES(ri-buildfasta divsufsort64)

add_executable(ri-locate ri-locate.cpp)
TARGET_LINK_LIBRARIES(ri-locate sdsl)
TARGET_LINK_LIBRARIES(ri-locate divsufsort)
TARGET_LINK_LIBRARIES(ri-locate divsufsort64)

add_executable(ri-count ri-count.cpp)
TARGET_LINK_LIBRARIES(ri-count sdsl)
TARGET_LINK_LIBRARIES(ri-count divsufsort)
TARGET_LINK_LIBRARIES(ri-count divsufsort64)

add_executable(ri-space ri-space.cpp)
TARGET_LINK_LIBRARIES(ri-space sdsl)
TARGET_LINK_LIBRARIES(ri-space divsufsort)
TARGET_LINK_LIBRARIES(ri-space divsufsort64)

add_executable(ri-align ri-align.cpp)
TARGET_LINK_LIBRARIES(ri-align sdsl)
TARGET_LINK_LIBRARIES(ri-align divsufsort)
TARGET_LINK_LIBRARIES(ri-align divsufsort64)
TARGET_LINK_LIBRARIES(ri-align z)

add_executable(check_bwtsa bwtsa.cpp)
TARGET_LINK_LIBRARIES(check_bwtsa sdsl)
