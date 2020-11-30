

set(NAME "psspy")

message(STATUS "Building: ${NAME}")

configure_file(manpage.1 ${CMAKE_MAN_OUTPUT_DIRECTORY}/${NAME}.1 COPYONLY)
configure_file(pss.py ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pss.py COPYONLY)

##############################################################
# Deployment
##############################################################
install(PROGRAMS pss.py
        DESTINATION bin
)
install(FILES manpage.1 TYPE MAN)