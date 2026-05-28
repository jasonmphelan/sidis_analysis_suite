# richMisid/CMakeLists.txt
#
# Builds two executables:
#
#   richKaonMisidHists  – Stage 1: fill m²_RICH histograms from ePi tree
#   richKaonMisidCorr   – Stage 2: fit pion peak, output TH2F + diagnostic PDFs
#
# Add this subdirectory from the top-level CMakeLists.txt with:
#   add_subdirectory(richMisid)

include_directories(
    ${ROOT_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/include
    ${CLAS12ROOT}/Clas12Banks
    ${CLAS12ROOT}/Clas12Root
    ${INCLUDES}
)

# ── Stage 1: histogram filling ─────────────────────────────────────────────────
# Needs pion/electron classes (read from ePi TTree) → links clashit_tree
add_executable(richKaonMisidHists
    richKaonMisidHists.cpp
    ${CMAKE_SOURCE_DIR}/classes/pion.cpp
    ${CMAKE_SOURCE_DIR}/classes/electron.cpp
)
target_link_libraries(richKaonMisidHists
    ${ROOT_LIBRARIES}
    clashit_tree
)

# ── Stage 2: fitting and correction output ─────────────────────────────────────
# Pure ROOT analysis; no custom classes needed at runtime.
add_executable(richKaonMisidCorr
    richKaonMisidCorr.cpp
)
target_link_libraries(richKaonMisidCorr
    ${ROOT_LIBRARIES}
)

# ── Stage 3: pidRICH pion-fraction histograms vs probICH cut ──────────────────
add_executable(pidRICHHists
    pidRICHHists.cpp
    ${CMAKE_SOURCE_DIR}/classes/pion.cpp
    ${CMAKE_SOURCE_DIR}/classes/electron.cpp
    ${CMAKE_SOURCE_DIR}/classes/analyzer.cpp
    ${CMAKE_SOURCE_DIR}/classes/DCfid_SIDIS.cpp
)
target_link_libraries(pidRICHHists
    ${ROOT_LIBRARIES}
    clashit_tree
)

# ── Stage 4: plot corrections vs probICH and momentum ──────────────────────────
add_executable(plotCorrVsProb
    plotCorrVsProb.cpp
    ${CMAKE_SOURCE_DIR}/classes/pion.cpp
    ${CMAKE_SOURCE_DIR}/classes/electron.cpp
)
target_link_libraries(plotCorrVsProb
    ${ROOT_LIBRARIES}
    clashit_tree
)
