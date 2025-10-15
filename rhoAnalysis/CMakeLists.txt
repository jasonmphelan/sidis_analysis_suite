#include_directories(${INCLUDES} ${LIBRARIES})
include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})

add_executable(mergeChunks mergeChunks.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/CLAS.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp)
target_link_libraries(mergeChunks ${ROOT_LIBRARIES} clashit_tree) 

add_executable(rotateRho rotateRho_old.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/CLAS.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp)
target_link_libraries(rotateRho ${ROOT_LIBRARIES} clashit_tree) 

add_executable(rotateRhoNew rotateRho.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/CLAS.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp)
target_link_libraries(rotateRhoNew ${ROOT_LIBRARIES} clashit_tree) 

add_executable(analyzeSidebands analyzeSidebands.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/CLAS.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp)
target_link_libraries(analyzeSidebands ${ROOT_LIBRARIES} clashit_tree) 

add_executable(plotXS plotDependence.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp)
target_link_libraries(plotXS ${ROOT_LIBRARIES} clashit_tree) 

add_executable(norms computeNorms.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp)
target_link_libraries(norms ${ROOT_LIBRARIES} clashit_tree) 

add_executable(studyBack studyBack.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/CLAS.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp)
target_link_libraries(studyBack ${ROOT_LIBRARIES} clashit_tree) 

add_executable(studyWeights studyWeights.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/CLAS.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp)
target_link_libraries(studyWeights ${ROOT_LIBRARIES} clashit_tree) 

