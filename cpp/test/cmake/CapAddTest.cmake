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
    endiF()
    foreach(PROCS ${NUMBER_OF_PROCESSES_TO_EXECUTE})
        add_test(
            NAME ${TEST_NAME}_cpp_${PROCS}_procs
            COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ./${TEST_NAME}.exe
        )
    endforeach()
endfunction()

FUNCTION(COPY_CAP_INPUT_FILE INPUT_FILE)
    ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${INPUT_FILE}
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/data/${INPUT_FILE}
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/data/${INPUT_FILE} ${CMAKE_CURRENT_BINARY_DIR}/${INPUT_FILE}
        COMMENT "Copying ${INPUT_FILE}"
    )
    ADD_CUSTOM_TARGET(
        ${INPUT_FILE} ALL
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${INPUT_FILE}
    )
ENDFUNCTION()
