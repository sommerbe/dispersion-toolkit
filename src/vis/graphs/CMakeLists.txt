

set(NAME "vis${PROGRAM_DELIMITER}graphs")

message(STATUS "Building: ${NAME}")

configure_file(graphs.py ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME}.py COPYONLY)

##############################################################
# Deployment
##############################################################
install(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME}.py
        DESTINATION bin
)