 IF (WIN32)
 	INCLUDE_DIRECTORIES(C:/Program Files/PCL/include/pcl-1.8)
 ELSEIF (UNIX)
 	INCLUDE_DIRECTORIES(/usr/local/include/pcl-1.8)
 ENDIF ()
 
 ADD_SUBDIRECTORY(thin_plate_splines)
 ADD_SUBDIRECTORY(optimal_nonrigid_icp)

 SET(SRC_LIST main.cpp)

 ADD_EXECUTABLE(Non-Rigid-Registar ${SRC_LIST})
 SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

 #TARGET_LINK_LIBRARIES(Non-Rigid-Registar)

 INSTALL(TARGETS Non-Rigid-Registar RUNTIME DESTINATION bin)