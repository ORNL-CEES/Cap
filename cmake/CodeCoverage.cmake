## C++ coverage ###############################################################
find_program(LCOV_EXECUTABLE lcov)
if(LCOV_EXECUTABLE)
    message("-- Found lcov: ${LCOV_EXECUTABLE}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
else()
    message(FATAL_ERROR "-- lcov not found")
endif()
find_program(GENHTML_EXECUTABLE genhtml)
if(GENHTML_EXECUTABLE)
    message("-- Found genhtml: ${GENHTML_EXECUTABLE}")
else()
    message(FATAL_ERROR "-- genhtml not found")
endif()
set(CPP_COVERAGE_FILE ${CMAKE_BINARY_DIR}/lcov.info)
set(CPP_COVERAGE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/htmlcov-cpp)
add_custom_target(coverage-cpp
    COMMAND ${LCOV_EXECUTABLE} 
    --capture
    --directory ${CMAKE_BINARY_DIR}
    --output-file=${CPP_COVERAGE_FILE}
    COMMAND ${GENHTML_EXECUTABLE} ${CPP_COVERAGE_FILE}
    --output-directory ${CPP_COVERAGE_OUTPUT_DIRECTORY}
)

## Python coverage ############################################################
if(ENABLE_PYTHON)
    find_program(COVERAGE_EXECUTABLE coverage)
    if(COVERAGE_EXECUTABLE)
        message("-- Found coverage: ${COVERAGE_EXECUTABLE}")
    else()
        message(FATAL_ERROR "-- coverage not found")
    endif()
    set(PYTHON_COVERAGE_FILE ${CMAKE_BINARY_DIR}/coverage.xml)
    set(PYTHON_COVERAGE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/htmlcov-python)
    add_custom_target(coverage-python
        COMMAND ${COVERAGE_EXECUTABLE} xml
        --ignore-errors -o ${PYTHON_COVERAGE_FILE}
        COMMAND ${COVERAGE_EXECUTABLE} html
        --ignore-errors -d ${PYTHON_COVERAGE_OUTPUT_DIRECTORY}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/python/test
    )
endif()
