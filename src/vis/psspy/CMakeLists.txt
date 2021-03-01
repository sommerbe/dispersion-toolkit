

set(NAME "vis${PROGRAM_DELIMITER}psspy")

message(STATUS "Building: ${NAME}")

configure_file(pss.py ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME}.py COPYONLY)

##############################################################
# Deployment
##############################################################
install(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME}.py
        DESTINATION bin
)