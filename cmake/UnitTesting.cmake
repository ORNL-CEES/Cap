#### C++ #####################################################################
function(Cap_ADD_BOOST_TEST TEST_NAME)
    add_executable(${TEST_NAME}.exe ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.cc)
    target_link_libraries(${TEST_NAME}.exe Cap)
    target_link_libraries(${TEST_NAME}.exe ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    target_link_libraries(${TEST_NAME}.exe ${Boost_SYSTEM_LIBRARY})
    target_link_libraries(${TEST_NAME}.exe ${Boost_TIMER_LIBRARY})
    target_link_libraries(${TEST_NAME}.exe ${Boost_CHRONO_LIBRARY})
    set_target_properties(${TEST_NAME}.exe PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED ON
    )
    if(ARGN)
        set(NUMBER_OF_PROCESSES_TO_EXECUTE ${ARGN})
    else()
        set(NUMBER_OF_PROCESSES_TO_EXECUTE 1)
    endif()
    foreach(NPROC ${NUMBER_OF_PROCESSES_TO_EXECUTE})
        add_test(
            NAME ${TEST_NAME}_cpp_${NPROC}
            COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} ./${TEST_NAME}.exe
        )
        set_tests_properties(${TEST_NAME}_cpp_${NPROC} PROPERTIES
            PROCESSORS ${NPROC}
        )
    endforeach()
endfunction()

function(Cap_ADD_CPP_EXAMPLE EXAMPLE_NAME)
    add_executable(${EXAMPLE_NAME}.exe ${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE_NAME}.cc)
    target_link_libraries(${EXAMPLE_NAME}.exe Cap)
    set_target_properties(${EXAMPLE_NAME}.exe PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED ON
    )
endfunction()

#### Python ##################################################################
function(Cap_ADD_PYTHON_TEST TEST_NAME)
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.py
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.py ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py
        COMMENT "Copying ${TEST_NAME}.py"
    )
    add_custom_target(
        ${TEST_NAME}.py ALL
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py
    )
    if(ARGN)
        set(NUMBER_OF_PROCESSES_TO_EXECUTE ${ARGN})
    else()
        set(NUMBER_OF_PROCESSES_TO_EXECUTE 1)
    endif()
    foreach(NPROC ${NUMBER_OF_PROCESSES_TO_EXECUTE})
        if(ENABLE_COVERAGE)
            add_test(
                NAME ${TEST_NAME}_py_${NPROC}
                COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} ${COVERAGE_EXECUTABLE} run --append ${TEST_NAME}.py
            )
        else()
            add_test(
                NAME ${TEST_NAME}_py_${NPROC}
                COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} ${PYTHON_EXECUTABLE} ${TEST_NAME}.py
            )
        endif()
        set_tests_properties(${TEST_NAME}_py_${NPROC} PROPERTIES
            ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}"
            PROCESSORS ${NPROC}
        )
    endforeach()
endfunction()

#### Copy input files ########################################################
function(Cap_COPY_INPUT_FILE INPUT_FILE PATH_TO_FILE)
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${INPUT_FILE}
        DEPENDS ${CMAKE_SOURCE_DIR}/${PATH_TO_FILE}/${INPUT_FILE}
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy ${CMAKE_SOURCE_DIR}/${PATH_TO_FILE}/${INPUT_FILE} ${CMAKE_CURRENT_BINARY_DIR}/${INPUT_FILE}
        COMMENT "Copying ${INPUT_FILE}"
    )
    string(REGEX REPLACE "/" "_" DUMMY ${CMAKE_CURRENT_BINARY_DIR}/${INPUT_FILE})
    add_custom_target(
       ${DUMMY} ALL
       DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${INPUT_FILE}
    )
endfunction()
