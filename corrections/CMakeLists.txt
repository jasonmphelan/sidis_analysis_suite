#include_directories(${INCLUDES} ${LIBRARIES})
include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})

add_executable(fitCorrection fitCorrection.cpp  )
target_link_libraries(fitCorrection ${ROOT_LIBRARIES} ) 

#add_executable(fitWeight fitWeight.cpp  )
#target_link_libraries(fitWeight ${ROOT_LIBRARIES} ) 

#add_executable(evalCorrection evalCorrectionVals.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp )
#target_link_libraries(evalCorrection ${ROOT_LIBRARIES} clashit_tree ) 

#add_executable(calcCorrections calcCorrections.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/genElectron.cpp ${CMAKE_SOURCE_DIR}/classes/genPion.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp) 
#target_link_libraries(calcCorrections ${ROOT_LIBRARIES} clashit_tree genhit_tree) 

#add_executable(testCorrections testCorrections.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp )
#target_link_libraries(testCorrections ${ROOT_LIBRARIES} clashit_tree ) 

add_executable(calcCorrections calcCorrections.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp 
							${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp 
							${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp
							${CMAKE_SOURCE_DIR}/classes/genPion.cpp ${CMAKE_SOURCE_DIR}/classes/genElectron.cpp)

target_link_libraries(calcCorrections ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree genhit_tree)
