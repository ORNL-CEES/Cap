include(${CMAKE_SOURCE_DIR}/cmake/UnitTesting.cmake)

include_directories(${CMAKE_SOURCE_DIR}/cpp/source/dummy)
include_directories(${CMAKE_SOURCE_DIR}/cpp/source/deal.II/dummy)

Cap_ADD_CPP_EXAMPLE(scaling)

Cap_COPY_INPUT_FILE(super_capacitor.info cpp/example)
