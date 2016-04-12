#### cap-data ################################################################
if(NOT DEFINED CAP_DATA_DIR)
    message(FATAL_ERROR "\n"
        "You must pass a flag -DCAP_DATA_DIR=/path/to/dir/with/meshes\n"
        "to cmake.\n"
        )
endif()

#### Message Passing Interface (MPI) #########################################
find_package(MPI REQUIRED)

#### Boost ###################################################################
if(DEFINED BOOST_DIR)
    set(BOOST_ROOT ${BOOST_DIR})
endif()
set(Boost_COMPONENTS
    mpi
    serialization
    unit_test_framework
    chrono
    timer
    system
)
if(ENABLE_PYTHON)
    set(Boost_COMPONENTS python ${Boost_COMPONENTS})
endif()
find_package(Boost 1.59.0 REQUIRED COMPONENTS ${Boost_COMPONENTS})

#### Python ##################################################################
if(ENABLE_PYTHON)
    find_package(PythonInterp REQUIRED)
    find_package(PythonLibs REQUIRED)
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -c "from mpi4py import get_include; print(get_include())"
        OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR
    )
    find_path(MPI4PY_INCLUDE_DIR mpi4py/mpi4py.h 
        PATH ${MPI4PY_INCLUDE_DIR}
        NO_DEFAULT
    )
    if(MPI4PY_INCLUDE_DIR)
        message("MPI4PY_INCLUDE_DIR=${MPI4PY_INCLUDE_DIR}")
    else()
        message(FATAL_ERROR "mpi4py not found.")
    endif()
endif()

#### deal.II #################################################################
if(ENABLE_DEAL_II)
    find_package(deal.II 8.4 REQUIRED PATHS ${DEAL_II_DIR})
    add_definitions(-DWITH_DEAL_II)
endif()

#### GNU Scientific Library ##################################################
if(ENABLE_GSL)
    if(DEFINED GSL_DIR)
        set(GSL_ROOT_DIR ${GSL_DIR})
    endif()
    find_package(GSL REQUIRED)
    add_definitions(-DWITH_GSL)
ENDIF()
