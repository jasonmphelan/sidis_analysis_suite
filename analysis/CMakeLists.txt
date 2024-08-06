include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})
#include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})

#add_executable(kaons makeKaonPlots.cpp ${CMAKE_SOURCE_DIR}/classes/clashit.cpp)
#target_link_libraries(kaons ${ROOT_LIBRARIES}) 

add_executable(makeRatio makeRatioBinned.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp )
target_link_libraries(makeRatio ${ROOT_LIBRARIES} clashit_tree) 

