#include_directories(${INCLUDES} )
include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})

add_executable(detector_skimmer test_skimmer.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp 
							${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp 
							${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp
							${CMAKE_SOURCE_DIR}/classes/genPion.cpp ${CMAKE_SOURCE_DIR}/classes/genElectron.cpp
							)

target_link_libraries(detector_skimmer ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree genhit_tree)

add_executable(final_skimmer final_skimmer.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp)

target_link_libraries(final_skimmer ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree)

add_executable(rho_skimmer rho_skimmer.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp ${CMAKE_SOURCE_DIR}/classes/genElectron.cpp ${CMAKE_SOURCE_DIR}/classes/genPion.cpp )

target_link_libraries(rho_skimmer ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree genhit_tree)

add_executable(skimmer_plotter skimmer_plots.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp 
							${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp 
							${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp)
target_link_libraries(skimmer_plotter ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree )


add_executable(kaon_skimmer kaon_skimmer.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp 
							${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp 
							${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp
							${CMAKE_SOURCE_DIR}/classes/genPion.cpp ${CMAKE_SOURCE_DIR}/classes/genElectron.cpp)

target_link_libraries(kaon_skimmer ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree genhit_tree)

add_executable(gen_skimmer gen_skimmer.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp ${CMAKE_SOURCE_DIR}/classes/genElectron.cpp ${CMAKE_SOURCE_DIR}/classes/genPion.cpp )

target_link_libraries(gen_skimmer ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree genhit_tree)

add_executable(inclusiveStats inclusive_stats.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp 
							${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp 
							${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp
							${CMAKE_SOURCE_DIR}/classes/genPion.cpp ${CMAKE_SOURCE_DIR}/classes/genElectron.cpp)

target_link_libraries(inclusiveStats ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree genhit_tree)

add_executable(reskim reskim.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp)

target_link_libraries(reskim ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree)

add_executable(getCharge getCharge.cpp ${CMAKE_SOURCE_DIR}/classes/electron.cpp ${CMAKE_SOURCE_DIR}/classes/pion.cpp 
							${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp ${CMAKE_SOURCE_DIR}/classes/e_pid.cpp 
							${CMAKE_SOURCE_DIR}/classes/analyzer.cpp ${CMAKE_SOURCE_DIR}/classes/reader.cpp
							${CMAKE_SOURCE_DIR}/classes/genPion.cpp ${CMAKE_SOURCE_DIR}/classes/genElectron.cpp)
						
target_link_libraries(getCharge ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root clashit_tree genhit_tree)
