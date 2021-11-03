cmake_minimum_required(VERSION 3.0)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/.. )

project (s-rlbwt)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/internal)
# add_subdirectory(${PROJECT_SOURCE_DIR}/Big-BWT)

set(CMAKE_INSTALL_INCLUDEDIR "include") # This is an hack because include(GUIInstallDirs) doesn't work
include(FetchContent)
## Add Big-BWT
FetchContent_Declare(
  bigbwt
  GIT_REPOSITORY https://github.com/alshai/Big-BWT.git
  )

FetchContent_GetProperties(bigbwt)
if(NOT bigbwt_POPULATED)
  FetchContent_Populate(bigbwt)
  add_subdirectory(${bigbwt_SOURCE_DIR} ${bigbwt_BINARY_DIR})

  endif()

# check for SDSL
find_library(SDSL_LIB sdsl)
if (NOT SDSL_LIB)
    message("sdsl-lite not found, building as external content")

    ## Add sdsl
    FetchContent_Declare(
      sdsl
      GIT_REPOSITORY https://github.com/simongog/sdsl-lite
    )

    FetchContent_GetProperties(sdsl)
    if(NOT sdsl_POPULATED)
      FetchContent_Populate(sdsl)

      set(GENERATE_DOC OFF CACHE BOOL "Do not generate doxygen for sdsl-lite")

      add_subdirectory(${sdsl_SOURCE_DIR} ${sdsl_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()

    ## Add divsuffsort
    FetchContent_Declare(
      divsufsort
      GIT_REPOSITORY https://github.com/simongog/libdivsufsort.git
      GIT_TAG        2.0.1
    )

    FetchContent_GetProperties(divsufsort)
    if(NOT divsufsort_POPULATED)
      FetchContent_Populate(divsufsort)

      set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not build a shared library for libdivsufsort")
      set(BUILD_EXAMPLES OFF CACHE BOOL "Do not build libdivsufsort example")
      set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build libdivsufsort in 64-bits mode")

      add_subdirectory(${divsufsort_SOURCE_DIR} ${divsufsort_BINARY_DIR} EXCLUDE_FROM_ALL)

      target_include_directories(divsufsort PUBLIC "${divsufsort_BINARY_DIR}/include")
      target_include_directories(divsufsort64 PUBLIC "${divsufsort_BINARY_DIR}/include")
    endif()
else()
    find_library(DSS divsufsort)
    find_library(DSS64 divsufsort64)
    message("sdsl: " ${SDSL_LIB})
    message("divsufsort: " ${DSS})
    message("divsufsort64: " ${DSS64})
endif()

## Add klib
FetchContent_Declare(
  klib
  GIT_REPOSITORY https://github.com/attractivechaos/klib
)

FetchContent_GetProperties(klib)
if(NOT klib_POPULATED)
  FetchContent_Populate(klib)

  # add_subdirectory(${klib_SOURCE_DIR} ${klib_BINARY_DIR} EXCLUDE_FROM_ALL)
  add_library(klib INTERFACE)

  target_include_directories(klib INTERFACE ${klib_SOURCE_DIR})
endif()




message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS "--std=c++17")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-g -Ofast -fstrict-aliasing -march=native -DNDEBUG")
# set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fstrict-aliasing -march=native -DNDEBUG -DLF_STATS")  # Uncomment this line in order to enable LF_STATS.
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")


add_executable(ri-build EXCLUDE_FROM_ALL ri-build.cpp)
TARGET_LINK_LIBRARIES(ri-build sdsl)
TARGET_LINK_LIBRARIES(ri-build divsufsort)
TARGET_LINK_LIBRARIES(ri-build divsufsort64)

add_executable(ri-locate EXCLUDE_FROM_ALL ri-locate.cpp)
TARGET_LINK_LIBRARIES(ri-locate sdsl)
TARGET_LINK_LIBRARIES(ri-locate divsufsort)
TARGET_LINK_LIBRARIES(ri-locate divsufsort64)

add_executable(ri-count EXCLUDE_FROM_ALL ri-count.cpp)
TARGET_LINK_LIBRARIES(ri-count sdsl)
TARGET_LINK_LIBRARIES(ri-count divsufsort)
TARGET_LINK_LIBRARIES(ri-count divsufsort64)

add_executable(ri-space EXCLUDE_FROM_ALL ri-space.cpp)
TARGET_LINK_LIBRARIES(ri-space sdsl)
TARGET_LINK_LIBRARIES(ri-space divsufsort)
TARGET_LINK_LIBRARIES(ri-space divsufsort64)
TARGET_LINK_LIBRARIES(ri-space z)
TARGET_LINK_LIBRARIES(ri-space klib)

add_executable(ri-tree-test ri-rasa-tree-test.cpp)
TARGET_LINK_LIBRARIES(ri-tree-test sdsl)
TARGET_LINK_LIBRARIES(ri-tree-test divsufsort)
TARGET_LINK_LIBRARIES(ri-tree-test divsufsort64)
TARGET_LINK_LIBRARIES(ri-tree-test z)
TARGET_LINK_LIBRARIES(ri-tree-test klib)

add_executable(ri-query-rasa ri-query-rasa.cpp)
TARGET_LINK_LIBRARIES(ri-query-rasa sdsl)
TARGET_LINK_LIBRARIES(ri-query-rasa divsufsort)
TARGET_LINK_LIBRARIES(ri-query-rasa divsufsort64)
TARGET_LINK_LIBRARIES(ri-query-rasa z)
TARGET_LINK_LIBRARIES(ri-query-rasa klib)

add_executable(ri-query ri-query.cpp)
TARGET_LINK_LIBRARIES(ri-query sdsl)
TARGET_LINK_LIBRARIES(ri-query divsufsort)
TARGET_LINK_LIBRARIES(ri-query divsufsort64)
TARGET_LINK_LIBRARIES(ri-query z)
TARGET_LINK_LIBRARIES(ri-query klib)

add_executable(ri-buildfasta ri-buildfasta.cpp)
TARGET_LINK_LIBRARIES(ri-buildfasta sdsl)
TARGET_LINK_LIBRARIES(ri-buildfasta divsufsort)
TARGET_LINK_LIBRARIES(ri-buildfasta divsufsort64)
TARGET_LINK_LIBRARIES(ri-buildfasta z)
TARGET_LINK_LIBRARIES(ri-buildfasta klib)

add_executable(ri-align ri-align.cpp)
TARGET_LINK_LIBRARIES(ri-align sdsl)
TARGET_LINK_LIBRARIES(ri-align divsufsort)
TARGET_LINK_LIBRARIES(ri-align divsufsort64)
TARGET_LINK_LIBRARIES(ri-align z)
TARGET_LINK_LIBRARIES(ri-align klib)


install(TARGETS ri-align DESTINATION bin)
install(TARGETS ri-buildfasta DESTINATION bin)
