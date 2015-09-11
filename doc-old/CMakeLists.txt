FIND_PACKAGE(LATEX)
IF(PDFLATEX_COMPILER) 
    INCLUDE(${PROJECT_SOURCE_DIR}/cmake/UseLATEX.cmake)

    ADD_LATEX_DOCUMENT(doc.tex BIBFILES mybib.bib DEFAULT_PDF)
    FILE(COPY
        ${CAP_DATA_DIR}/figures/series_of_constant_current_charge_and_constant_voltage_discharge.png
        ${CAP_DATA_DIR}/figures/constant_current_cycling_current_and_voltage.png
        ${CAP_DATA_DIR}/figures/constant_current_cycling_power_density.png
        ${CAP_DATA_DIR}/figures/constant_current_cycling_temperature.png
        ${CAP_DATA_DIR}/figures/constant_current_charge_constant_voltage_discharge_current_and_voltage.png
        ${CAP_DATA_DIR}/figures/constant_current_charge_constant_voltage_discharge_power_density.png
        ${CAP_DATA_DIR}/figures/constant_current_charge_constant_voltage_discharge_temperature.png
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/figures
    )
    MESSAGE("${CAP_DATA_DIR}/figures/")
ENDIF()
