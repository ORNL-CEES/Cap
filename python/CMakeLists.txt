SET(PYTHON_INCLUDE_DIRECTORIES /opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7)
SET(PYTHON_LIBRARIES           /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/libpython2.7.dylib)
SET(PYTHON_EXECUTABLE          /opt/local/bin/python2.7)

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRECTORIES})
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "")
ADD_LIBRARY(pycap MODULE
    python_wrappers.cc
    )

SET_TARGET_PROPERTIES(pycap PROPERTIES 
    PREFIX ""
    SUFFIX ".so"
    )

TARGET_LINK_LIBRARIES(pycap ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} cap deal_II)
INSTALL(TARGETS pycap DESTINATION ${CMAKE_INSTALL_PREFIX}/python)

ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/test_python_wrappers.py
    DEPENDS ${PROJECT_SOURCE_DIR}/python/test_python_wrappers.py
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy ${PROJECT_SOURCE_DIR}/python/test_python_wrappers.py ${CMAKE_CURRENT_BINARY_DIR}/test_python_wrappers.py
    COMMENT "Copying test_python_wrappers.py"
)
ADD_CUSTOM_TARGET(
    test_python_wrappers.py ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/test_python_wrappers.py
    )
ADD_TEST(test_python_wrappers ${PYTHON_EXECUTABLE} test_python_wrappers.py)
