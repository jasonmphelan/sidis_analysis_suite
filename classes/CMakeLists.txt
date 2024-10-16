# find source and header files, set paths, includes...
SET(LIBNAME "clashit_tree")
FILE(GLOB SRCS electron.cpp pion.cpp)

FILE(GLOB HEADERS ${CMAKE_SOURCE_DIR}/include/electron.h ${CMAKE_SOURCE_DIR}/include/pion.h)
set(INCLUDES ${ROOT_INCLUDE_DIR} "${CMAKE_SOURCE_DIR}/include/"  )
include_directories(${CLAS12ROOT}/Clas12Banks ${CLAS12ROOT}/Clas12Root ${INCLUDES} ${LIBRARIES})

# create root dictionary
SET(LINKDEF_e ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef_e.h)
SET(LINKDEF_gen_e ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef_gen_e.h)
SET(LINKDEF_gen_pi ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef_gen_pi.h)
SET(LINKDEF_pi ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef_pi.h)
SET(LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h)

SET(DICTIONARY ${LIBNAME}Dict.cpp)
ROOT_GENERATE_DICTIONARY(G__electron ${INCLUDES}/electron.h MODULE electrons LINKDEF ${LINKDEF_e})
ROOT_GENERATE_DICTIONARY(G__genElectron ${INCLUDES}/genElectron.h MODULE genElectrons LINKDEF ${LINKDEF_gen_e})
ROOT_GENERATE_DICTIONARY(G__genPion ${INCLUDES}/genPion.h MODULE genPions LINKDEF ${LINKDEF_gen_pi})
ROOT_GENERATE_DICTIONARY(G__pion ${INCLUDES}/pion.h MODULE pions LINKDEF ${LINKDEF_pi})
ROOT_GENERATE_DICTIONARY(G__vectors TLorentzVector.h MODULE vectors_tree LINKDEF ${LINKDEF})
#LIST(APPEND SRCS ${DICTIONARY})

# build the library
ADD_LIBRARY(${LIBNAME} SHARED ${SRCS} G__electron G__pion G__vectors ${ROOT_LIBRARIES}  )
ADD_LIBRARY(genhit_tree SHARED ${SRCS} G__genElectron G__genPion ${ROOT_LIBRARIES}  )
TARGET_LINK_LIBRARIES(${LIBNAME} ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root   )
TARGET_LINK_LIBRARIES(genhit_tree ${ROOT_LIBRARIES} -lEG PkgConfig::hipo4 -L${CLAS12ROOT}/lib -lClas12Banks -lClas12Root   )



#INSTALL ( TARGETS ${LIBNAME} DESTINATION "${CMAKE_SOURCE_DIR}/lib")

#INSTALL ( FILES ${CMAKE_CURRENT_BINARY_DIR}/${LIBNAME}Dict_rdict.pcm DESTINATION ${CMAKE_SOURCE_DIR}/lib)
