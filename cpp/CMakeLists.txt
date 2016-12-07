add_subdirectory(source)
add_subdirectory(test)
add_subdirectory(example)

include_directories(${Cap_INCLUDE_DIRS})
add_library(Cap ${Cap_SOURCES})
target_include_directories(Cap SYSTEM PUBLIC ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(Cap PUBLIC ${MPI_CXX_LIBRARIES})
target_include_directories(Cap SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(Cap PUBLIC ${Boost_MPI_LIBRARY})
target_link_libraries(Cap PUBLIC ${Boost_SERIALIZATION_LIBRARY})
target_link_libraries(Cap PUBLIC ${Boost_CHRONO_LIBRARY})
target_link_libraries(Cap PUBLIC ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries(Cap PUBLIC ${Boost_REGEX_LIBRARY})
if(ENABLE_DEAL_II)
    target_include_directories(Cap SYSTEM PUBLIC ${DEAL_II_INCLUDE_DIRS})
    target_link_libraries(Cap PUBLIC ${DEAL_II_LIBRARIES})
endif()
set_target_properties(Cap PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
    VERSION ${Cap_VERSION}
)

install(FILES ${Cap_HEADERS}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/cap
)
install(TARGETS Cap EXPORT CapTargets
    LIBRARY DESTINATION lib
    INCLUDES DESTINATION include
)
install(EXPORT CapTargets
    FILE CapTargets.cmake
    DESTINATION lib/cmake/Cap
)
