#include_directories(${INCLUDES} )
include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})


add_executable(makeThetaPhi makeThetaPhiPlots.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp)

target_link_libraries(makeThetaPhi ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree)

add_executable(makeAcceptanceMap makeAcceptanceMap.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp)

target_link_libraries(makeAcceptanceMap ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree)
