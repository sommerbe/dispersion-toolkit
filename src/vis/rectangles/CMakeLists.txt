

set(NAME "vis${PROGRAM_DELIMITER}rectangles")

message(STATUS "Building: ${NAME}")

configure_file(rectangles.py ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME}.py COPYONLY)

##############################################################
# Deployment
##############################################################
install(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME}.py
        DESTINATION bin
)