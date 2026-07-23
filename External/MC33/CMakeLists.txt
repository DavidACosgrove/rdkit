message("MC33 source dir ${CMAKE_CURRENT_SOURCE_DIR}")
if(NOT DEFINED MC33_DIR)
    set(MC33_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    set(fileToCheck "${MC33_DIR}/mc33/mc33cpp.h")
    set(needDownload "TRUE")
    if(EXISTS "${fileToCheck}")
        set(needDownload "FALSE")
    endif()
else()
    string(REGEX REPLACE "\\\\" "/" MC33_DIR ${MC33_DIR})
    set(needDownload "FALSE")
endif()

message("MC33 need download : ${needDownload}")

if(needDownload)
    if(NOT DEFINED MC33_URL)
        set(MC33_URL "https://github.com/dvega68/MC33_cpp_header-only/archive/main.tar.gz")
    endif()
    if(NOT DEFINED MC33_MD5SUM)
        set(MC33_MD5SUM "cfa234fc12fdb95cf911b2c274b0f793")
    endif()
    if(NOT DEFINED MC33_BASE)
        string(REGEX REPLACE "^.*/" "" MC33_BASE "${MC33_URL}")
    endif()
    downloadAndCheckMD5(${MC33_URL} "${MC33_DIR}/${MC33_BASE}" ${MC33_MD5SUM})
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf
            ${CMAKE_CURRENT_SOURCE_DIR}/main.tar.gz
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    execute_process(COMMAND mv ${CMAKE_CURRENT_SOURCE_DIR}/MC33_cpp_header-only-main mc33
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    execute_process(COMMAND rm ${CMAKE_CURRENT_SOURCE_DIR}/main.tar.gz
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
