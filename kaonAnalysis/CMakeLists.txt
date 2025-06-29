#include_directories(${INCLUDES} ${LIBRARIES})
include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})

add_executable(makeKaonPlots makeKaonPlots.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp )
target_link_libraries(makeKaonPlots ${ROOT_LIBRARIES} clashit_tree) 

add_executable(calcKaonCorr calcKaonCuts.cpp  )
target_link_libraries(calcKaonCorr ${ROOT_LIBRARIES} ) 

add_executable(fitKaonCorr fitKaonCorrection.cpp  )
target_link_libraries(fitKaonCorr ${ROOT_LIBRARIES} ) 

add_executable(makeDataRatios calcDataRatio.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp )
target_link_libraries(makeDataRatios ${ROOT_LIBRARIES} clashit_tree) 

add_executable(evalKaonCorrection evalCorrectionVals.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp )
target_link_libraries(evalKaonCorrection ${ROOT_LIBRARIES} clashit_tree ) 

add_executable(studyCorrections studyCorrections.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp)
target_link_libraries(studyCorrections ${ROOT_LIBRARIES} clashit_tree) 

add_executable(acceptanceMatchingKaons doAccMatchingKaon.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp)

target_link_libraries(acceptanceMatchingKaons ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree)

add_executable(mcKaonCorrections mcKaonCorrections.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp )
target_link_libraries(mcKaonCorrections ${ROOT_LIBRARIES} clashit_tree) 

add_executable(kaonSys kaonSys.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp )
target_link_libraries(kaonSys ${ROOT_LIBRARIES} clashit_tree) 
