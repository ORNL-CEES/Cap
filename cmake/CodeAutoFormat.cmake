## C++ format #################################################################
find_program(CLANG_FORMAT_EXECUTABLE clang-format)
if(CLANG_FORMAT_EXECUTABLE)
    message("-- Found clang-format: ${CLANG_FORMAT_EXECUTABLE}")
else()
    message(FATAL_ERROR "-- clang-format not found")
endif()
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/.clang-format
    DEPENDS ${CMAKE_SOURCE_DIR}/.clang-format
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy ${CMAKE_SOURCE_DIR}/.clang-format
        ${CMAKE_BINARY_DIR}/.clang-format
    COMMENT "Copying .clang-format"
)
add_custom_target(
    .clang-format
    DEPENDS ${CMAKE_BINARY_DIR}/.clang-format
)
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
    ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/diff-clang-format.py
        --file-extension='.h'
        --file-extension='.cc'
        ${CMAKE_SOURCE_DIR}/cpp
    DEPENDS
        ${CMAKE_BINARY_DIR}/.clang-format
        ${CMAKE_BINARY_DIR}/diff-clang-format.py
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
