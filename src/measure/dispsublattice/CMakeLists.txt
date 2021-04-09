
set(NAME "disp${PROGRAM_DELIMITER}sublattice")

set(SOURCES
  ${src}/io/argparse.cpp
  ${src}/io/ipointset.cpp
  ${src}/io/opointset.cpp
  ${src}/io/ostream.cpp
  ${src}/math/arithmetic.cpp
  ${src}/index/pointset_dsorted.cpp
)

set(HEADERS
  ${src}/io/argparse.hpp
  ${src}/io/ipointset.hpp
  ${src}/io/opointset.hpp
  ${src}/io/ostream.hpp
  ${src}/math/pointset.hpp
  ${src}/math/arithmetic.hpp
  ${src}/math/types.hpp
  ${src}/index/pointset_dsorted.hpp
)

message(STATUS "Building: ${NAME}")


add_executable(${NAME} main.cpp ${HEADERS} ${SOURCES})
target_link_libraries_common(${NAME})

if (OpenMP_FOUND)
  target_link_libraries(${NAME} OpenMP::OpenMP_CXX)
endif()

##############################################################
# Deployment
##############################################################
install(TARGETS ${NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
)