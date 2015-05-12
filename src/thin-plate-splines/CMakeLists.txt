ADD_EXECUTABLE(tpsdemo main.cpp linalg3d.h ludecomposition.h)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

IF (WIN32)
INCLUDE_DIRECTORIES(D:/boost_1_55_0 "C:/Program Files (x86)/Microsoft SDKs/Windows/v7.0A/Include" D:/glut-3.7/include)
LINK_DIRECTORIES("C:/Program Files (x86)/Microsoft SDKs/Windows/v7.0A/Lib/x64" "D:/glut-3.7/lib/x64")
TARGET_LINK_LIBRARIES(tpsdemo GLU32 OPENGL32 glut32)
ELSEIF (UNIX)
TARGET_LINK_LIBRARIES(tpsdemo GLU /usr/lib/nvidia-331/libGL.so glut)
ENDIF ()

ADD_EXECUTABLE(2d-morph 2d-morph.cpp 2d-morph.h linalg3d.h ludecomposition.h)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

ADD_EXECUTABLE(3d-morph 3d-morph.cpp 3d-morph.h linalg3d.h ludecomposition.h)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

INSTALL(TARGETS tpsdemo 2d-morph 3d-morph RUNTIME DESTINATION bin)