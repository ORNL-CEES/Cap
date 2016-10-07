set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/.ycm_extra_conf.py.in
    ${CMAKE_SOURCE_DIR}/.ycm_extra_conf.py
    @ONLY
)
