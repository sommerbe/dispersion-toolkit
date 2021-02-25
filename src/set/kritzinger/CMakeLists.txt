
set(NAME "set${PROGRAM_DELIMITER}kritzinger")

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
)

message(STATUS "Building: ${NAME}")


add_executable(${NAME} main.cpp ${HEADERS} ${SOURCES})
target_link_libraries_common(${NAME})

##############################################################
# Deployment
##############################################################
install(TARGETS ${NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
)