set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${PROJECT_SOURCE_DIR}/temp_libs)
 set(LIBCELL_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
 FIND_LIBRARY(LIBCELL_LIBRARY libcell)
 IF(LIBCELL_INCLUDE_DIR AND LIBCELL_LIBRARY)
 SET(LIBCELL_FOUND TRUE)
 ENDIF(LIBCELL_INCLUDE_DIR AND LIBCELL_LIBRARY)
