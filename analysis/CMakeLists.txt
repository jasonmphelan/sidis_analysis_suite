include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})
#include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})

#add_executable(kaons makeKaonPlots.cpp ${CMAKE_SOURCE_DIR}/classes/clashit.cpp)
#target_link_libraries(kaons ${ROOT_LIBRARIES}) 

#add_executable(makeRatio makeRatioBinned.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp )
#target_link_libraries(makeRatio ${ROOT_LIBRARIES} clashit_tree) 

add_executable(makeRatio makeRatioBinned_new.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp)
target_link_libraries(makeRatio ${ROOT_LIBRARIES} clashit_tree) 

add_executable(makeRatioBinned makeRatioBinned3D.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp)
target_link_libraries(makeRatioBinned ${ROOT_LIBRARIES} clashit_tree) 

add_executable(makeXsRatio makeXsRatioBinned.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/correctionTools.cpp)
target_link_libraries(makeXsRatio ${ROOT_LIBRARIES} clashit_tree) 


add_executable(makeKinematics makeKinematicPlots.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp )
target_link_libraries(makeKinematics ${ROOT_LIBRARIES} clashit_tree) 
