cmake_minimum_required(VERSION 3.2)

project(htd)

set(CPACK_GENERATOR "STGZ;TGZ;TZ;ZIP")
set(CPACK_SOURCE_GENERATOR "STGZ;TGZ;TZ;ZIP")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A small but efficient C++ library for computing (customized) tree and hypertree decompositions.")
set(CPACK_PACKAGE_VENDOR "Michael Abseher (abseher@dbai.tuwien.ac.at)")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "2")
set(CPACK_PACKAGE_VERSION_PATCH "0")

include(CPack)
include(CTest)
include(CheckIncludeFileCXX)
include(WriteCompilerDetectionHeader)
include(ExternalProject)

find_package(Git)
if(GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} update-index --refresh WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --always WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" RESULT_VARIABLE HTD_GIT_RETURN_VALUE OUTPUT_VARIABLE HTD_GIT_CURRENT_COMMIT_ID)
    if(NOT ${HTD_GIT_RETURN_VALUE} EQUAL 0)
        set(HTD_GIT_COMMIT_ID "")
        message(WARNING "Could not derive Git commit ID. Build will not contain Git revision info.")
    endif()

    string(REPLACE "\n" "" HTD_GIT_COMMIT_ID "${HTD_GIT_CURRENT_COMMIT_ID}")
else()
    set(HTD_GIT_COMMIT_ID "")
    message(WARNING "Could not find Git installation. Build will not contain Git revision info.")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Please choose the type of build. Options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE )
endif()

if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON)
endif()

if(DEFINED HTD_DEBUG_OUTPUT)
    if(HTD_DEBUG_OUTPUT)
        message("Debugging output is enabled!")

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHTD_DEBUG_OUTPUT" )
    endif()
endif()

if(DEFINED HTD_USE_EXTENDED_IDENTIFIERS)
    if(HTD_USE_EXTENDED_IDENTIFIERS)
        message("Extended identifiers will be used!")

        set(HTD_ID_TYPE "std::size_t" )
    else()
        set(HTD_ID_TYPE "std::uint_least32_t" )
    endif()
else()
    set(HTD_ID_TYPE "std::uint_least32_t" )
endif()

if(NOT DEFINED BUILD_TESTING)
    message("Tests are disabled!")

    set(BUILD_TESTING OFF)
endif()

if (NOT DEFINED BUILD_EXTERNAL_SOLVERS)
    ExternalProject_Add(mrprajesh
            GIT_REPOSITORY https://github.com/mrprajesh/pacechallenge
            GIT_TAG 6c29c143d72856f649de99846e91de185f78c15f
            CONFIGURE_COMMAND ""
            BUILD_COMMAND make -s
            INSTALL_COMMAND ""
            TEST_COMMAND ""
            BUILD_IN_SOURCE true
            SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/solvers/mrprajesh
            )
    get_filename_component(mrprajesh_path "${CMAKE_CURRENT_BINARY_DIR}/solvers/mrprajesh" ABSOLUTE)
    add_definitions(-DMRPRAJESH_PATH="${mrprajesh_path}")

    ExternalProject_Add(BZTreewidth
            GIT_REPOSITORY https://github.com/TomvdZanden/BZTreewidth
            GIT_TAG bfbf736a98eea412c7ba2de4a21280816d9b8bbb
            CONFIGURE_COMMAND ""
            BUILD_COMMAND make -s
            INSTALL_COMMAND ""
            TEST_COMMAND ""
            BUILD_IN_SOURCE true
            SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/solvers/BZTreewidth
            )
    get_filename_component(BZTreewidth_path "${CMAKE_CURRENT_BINARY_DIR}/solvers/BZTreewidth" ABSOLUTE)
    add_definitions(-DBZTREEWIDTH_PATH="${BZTreewidth_path}")

    ExternalProject_Add(Jdrasil
            GIT_REPOSITORY https://github.com/maxbannach/Jdrasil
            GIT_TAG fa7855e4c9f33163606a0677485a9e51d26d7b0a
            CONFIGURE_COMMAND ""
            BUILD_COMMAND make -s
            INSTALL_COMMAND ""
            TEST_COMMAND ""
            BUILD_IN_SOURCE true
            SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/solvers/Jdrasil
            )
    get_filename_component(Jdrasil_path "${CMAKE_CURRENT_BINARY_DIR}/solvers/Jdrasil" ABSOLUTE)
    add_definitions(-DJDRASIL_PATH="${Jdrasil_path}")

    ExternalProject_Add(mfjones
            GIT_REPOSITORY https://github.com/mfjones/pace2016
            GIT_TAG 2b7f289e4d182799803a014d0ee1d76a4de70c1f
            CONFIGURE_COMMAND ""
            BUILD_COMMAND make -s
            INSTALL_COMMAND ""
            TEST_COMMAND ""
            BUILD_IN_SOURCE true
            SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/solvers/mfjones
            )
    get_filename_component(mfjones_path "${CMAKE_CURRENT_BINARY_DIR}/solvers/mfjones" ABSOLUTE)
    add_definitions(-DMFJONES_PATH="${mfjones_path}")

    ExternalProject_Add(minfill_mrs
            GIT_REPOSITORY https://github.com/td-mrs/minfill_mrs
            GIT_TAG c1aa4f32ab3350d027ed0e932183f5f649b2ea88
            CONFIGURE_COMMAND ""
            BUILD_COMMAND make -s
            INSTALL_COMMAND ""
            TEST_COMMAND ""
            BUILD_IN_SOURCE true
            SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/solvers/minfill_mrs
            )
    get_filename_component(minfill_mrs_path "${CMAKE_CURRENT_BINARY_DIR}/solvers/minfill_mrs" ABSOLUTE)
    add_definitions(-DMINFILL_MRS_PATH="${minfill_mrs_path}")
else()

    add_definitions(-DBUILD_EXTERNAL_SOLVERS)
endif ()

subdirs(src/htd)
subdirs(src/htd_io)
subdirs(src/htd_cli)
subdirs(src/htd_main)

subdirs(test)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/templates/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/templates/Doxyfile.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile" @ONLY)

    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generate API documentation" VERBATIM
)
endif(DOXYGEN_FOUND)

if(MSVC OR CYGWIN)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/build/debug")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/build/release")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/build/debug")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/build/release")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/build/debug")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/build/release")
endif(MSVC OR CYGWIN)

add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
