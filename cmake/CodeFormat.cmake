## C++ format #################################################################
find_program(CLANG_FORMAT_EXECUTABLE clang-format)
if(CLANG_FORMAT_EXECUTABLE)
    message("-- Found clang-format: ${CLANG_FORMAT_EXECUTABLE}")
else()
    message(FATAL_ERROR "-- clang-format not found")
endif()
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/diff-clang-format.py
    DEPENDS ${CMAKE_SOURCE_DIR}/diff-clang-format.py
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy ${CMAKE_SOURCE_DIR}/diff-clang-format.py
        ${CMAKE_BINARY_DIR}/diff-clang-format.py
    COMMENT "Copying diff-clang-format.py"
)
add_custom_target(
    diff-clang-format.py
    DEPENDS ${CMAKE_BINARY_DIR}/diff-clang-format.py
)
add_custom_target(format-cpp
    ${CMAKE_BINARY_DIR}/diff-clang-format.py
        --file-extension='.h'
        --file-extension='.cc'
        --style=file
        --config=${CMAKE_SOURCE_DIR}/.clang-format
        --apply-patch
        ${CMAKE_SOURCE_DIR}/cpp
    DEPENDS
        ${CMAKE_BINARY_DIR}/diff-clang-format.py
)
file(WRITE
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/check_format_cpp.sh
    "#!/usr/bin/env bash\n"
    "\n"
    "${PYTHON_EXECUTABLE} "
    "${CMAKE_BINARY_DIR}/diff-clang-format.py "
    "--file-extension='.h' --file-extension='.cc' "
    "--style=file "
    "--config=${CMAKE_SOURCE_DIR}/.clang-format "
    "${CMAKE_SOURCE_DIR}/cpp"
)
file(COPY
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/check_format_cpp.sh
    DESTINATION
        ${CMAKE_BINARY_DIR}
    FILE_PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
)
add_test(
    NAME check_format
    COMMAND ${CMAKE_BINARY_DIR}/check_format_cpp.sh
)

## Python format ##############################################################
find_program(AUTOPEP8_EXECUTABLE autopep8)
if(AUTOPEP8_EXECUTABLE)
    message("-- Found autopep8: ${AUTOPEP8_EXECUTABLE}")
else()
    message(FATAL_ERROR "-- autopep8 not found")
endif()
add_custom_target(format-python
    ${AUTOPEP8_EXECUTABLE} 
    --diff ${CMAKE_SOURCE_DIR}/python/source/*.py
    COMMAND
    ${AUTOPEP8_EXECUTABLE} 
    --diff ${CMAKE_SOURCE_DIR}/python/test/*.py
)
