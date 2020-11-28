

set(NAME "psspy")

message(STATUS "Building: ${NAME}")

# add_custom_target(${NAME} SOURCES pss.py)

configure_file(pss.py ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pss.py COPYONLY)

##############################################################
# Deployment
##############################################################
install(PROGRAMS pss.py
        DESTINATION bin
)