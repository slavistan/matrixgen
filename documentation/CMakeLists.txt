# Project: Caesr
# Target: documentation
# Author: Stanislaw Hüll

# Generate doxygen config file
set(DOCCFG Doxyfile)
add_custom_target(doc-generate-doxyfile
  DEPENDS matrixgen
  COMMAND echo PROJECT_NAME=\"${CMAKE_PROJECT_NAME}\" > ${DOCCFG}
  COMMAND echo PROJECT_NUMBER=\"${CMAKE_PROJECT_VERSION}\" >> ${DOCCFG}
  COMMAND echo INPUT=\"${PROJECT_SOURCE_DIR}/include\" >> ${DOCCFG}
  COMMAND echo RECURSIVE=yes >> ${DOCCFG}
  COMMAND echo OUTPUT_DIRECTORY=\"${CMAKE_CURRENT_BINARY_DIR}\" >> ${DOCCFG}
  COMMAND echo TAB_SIZE=2 >> ${DOCCFG}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Run doxygen on the config file
add_custom_target(doc
  DEPENDS doc-generate-doxyfile
  COMMAND doxygen ${DOCCFG}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
