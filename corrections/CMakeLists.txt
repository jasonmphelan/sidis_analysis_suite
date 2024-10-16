#include_directories(${INCLUDES} ${LIBRARIES})
include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})

add_executable(fitCorrection fitCorrection.cpp  )
target_link_libraries(fitCorrection ${ROOT_LIBRARIES} ) 

add_executable(fitWeight fitWeight.cpp  )
target_link_libraries(fitWeight ${ROOT_LIBRARIES} ) 

add_executable(evalCorrection evalCorrectionVals.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp )
target_link_libraries(evalCorrection ${ROOT_LIBRARIES} clashit_tree ) 

#add_executable(calcCorrections calcCorrections.cpp)
#target_link_libraries(calcCorrections ${ROOT_LIBRARIES} clashit_tree) 

add_executable(testCorrections testCorrections.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp )
target_link_libraries(testCorrections ${ROOT_LIBRARIES} clashit_tree ) 

