
set(NAME "dispoptgs")

set(SOURCES
  ${src}/index/pointset_dsorted.cpp
  ${src}/io/argparse.cpp
  ${src}/io/ipointset.cpp
  ${src}/io/opointset.cpp
  ${src}/io/ostream.cpp
)

set(HEADERS
  ${src}/index/pointset_dsorted.hpp
  ${src}/io/argparse.hpp
  ${src}/io/ipointset.hpp
  ${src}/io/opointset.hpp
  ${src}/io/ostream.hpp
  ${src}/math/pointset.hpp
  ${src}/math/types.hpp
)

message(STATUS "Building: ${NAME}")

add_executable(${NAME} main.cpp ${HEADERS} ${SOURCES})
# target_link_libraries(${NAME} -lquadmath)

##############################################################
# Deployment
##############################################################
install(TARGETS ${NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
)