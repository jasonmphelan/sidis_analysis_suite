#include_directories(${INCLUDES} ${LIBRARIES})
include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})

add_executable(mergeChunks mergeChunks.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/CLAS.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp)
target_link_libraries(mergeChunks ${ROOT_LIBRARIES} clashit_tree) 

add_executable(rotateRho rotateRho_old.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/CLAS.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp)
target_link_libraries(rotateRho ${ROOT_LIBRARIES} clashit_tree) 

add_executable(rotateRhoSym rotateRho_sym.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/CLAS.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp)
target_link_libraries(rotateRhoSym ${ROOT_LIBRARIES} clashit_tree) 


add_executable(plotXS plotDependence.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp)
target_link_libraries(plotXS ${ROOT_LIBRARIES} clashit_tree) 

add_executable(norms computeNorms.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp)
target_link_libraries(norms ${ROOT_LIBRARIES} clashit_tree) 

add_executable(norms_sym computeNormsSym.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp)
target_link_libraries(norms_sym ${ROOT_LIBRARIES} clashit_tree) 

add_executable(studyBack studyBack.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/CLAS.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp  ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp)
target_link_libraries(studyBack ${ROOT_LIBRARIES} clashit_tree) 

add_executable(studyWeights studyWeights.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/CLAS.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp)
target_link_libraries(studyWeights ${ROOT_LIBRARIES} clashit_tree) 

add_executable(mcRhoCorrections mcRhoCorrections.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp 
							${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp 
							${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp
							${CMAKE_SOURCE_DIR}/classes/genPion.cpp ${CMAKE_SOURCE_DIR}/classes/genElectron.cpp)

target_link_libraries(mcRhoCorrections ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree genhit_tree)

add_executable(norms_3D computeNorms3D.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp)
target_link_libraries(norms_3D ${ROOT_LIBRARIES} clashit_tree)

add_executable(makeRhoKinematics makeRhoKinematicPlots.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp 
							${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp 
							${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp )
target_link_libraries(makeRhoKinematics ${ROOT_LIBRARIES} clashit_tree) 
