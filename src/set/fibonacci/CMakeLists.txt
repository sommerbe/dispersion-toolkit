
set(NAME "set${PROGRAM_DELIMITER}fibonacci")

set(SOURCES
  ${src}/io/argparse.cpp
  ${src}/io/opointset.cpp
  ${src}/io/ostream.cpp
)

set(HEADERS
  ${src}/io/argparse.hpp
  ${src}/io/opointset.hpp
  ${src}/io/ostream.hpp
  ${src}/math/pointset.hpp
  ${src}/math/types.hpp
  manpage.hpp
)

message(STATUS "Building: ${NAME}")

configure_file(manpage.1 ${CMAKE_MAN_OUTPUT_DIRECTORY}/${NAME}.1 COPYONLY)

add_executable(${NAME} main.cpp ${HEADERS} ${SOURCES})

# target_link_libraries(${NAME} -lquadmath)

##############################################################
# Deployment
##############################################################
install(TARGETS ${NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
)
install(FILES manpage.1 TYPE MAN)