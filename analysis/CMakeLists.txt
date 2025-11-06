include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})
#include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})

#add_executable(kaons makeKaonPlots.cpp ${CMAKE_SOURCE_DIR}/classes/clashit.cpp)
#target_link_libraries(kaons ${ROOT_LIBRARIES}) 

#add_executable(makeRatio makeRatioBinned.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp )
#target_link_libraries(makeRatio ${ROOT_LIBRARIES} clashit_tree) 

add_executable(makeRatio makeRatioBinned_new.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp)
target_link_libraries(makeRatio ${ROOT_LIBRARIES} clashit_tree) 

add_executable(studyZ studyZ.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp)
target_link_libraries(studyZ ${ROOT_LIBRARIES} clashit_tree) 


add_executable(makeRatioBinned makeRatioBinned3D.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp)
target_link_libraries(makeRatioBinned ${ROOT_LIBRARIES} clashit_tree) 

add_executable(makeXsRatio makeXsRatioBinned.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp)
target_link_libraries(makeXsRatio ${ROOT_LIBRARIES} clashit_tree) 


add_executable(makeKinematics makeKinematicPlots.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp 
							${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp 
							${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp )
target_link_libraries(makeKinematics ${ROOT_LIBRARIES} clashit_tree) 

add_executable(calcCharge calcCharge.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp 
							${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp 
							${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp
							${CMAKE_SOURCE_DIR}/classes/genPion.cpp ${CMAKE_SOURCE_DIR}/classes/genElectron.cpp)

target_link_libraries(calcCharge ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree genhit_tree)

add_executable(studyRuns studyRuns.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp 
							${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp 
							${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp
							${CMAKE_SOURCE_DIR}/classes/genPion.cpp ${CMAKE_SOURCE_DIR}/classes/genElectron.cpp)

target_link_libraries(studyRuns ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree genhit_tree)

add_executable(mergeEnergy mergeEnergy.cpp)

target_link_libraries(mergeEnergy ${ROOT_LIBRARIES})
