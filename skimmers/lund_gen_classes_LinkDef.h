//==============================================================================
//  lund_gen_classes_LinkDef.h
//
//  Dictionary directives for the standalone genElectron / genPion classes used
//  by gen_skimmer_lund.cpp.  The "+" requests the new-style (TStreamerInfo)
//  I/O, exactly as the original LinkDef_gen_e.h / LinkDef_gen_pi.h do, so the
//  generated StreamerInfo (class version 1) matches the original dictionary.
//==============================================================================
#if defined(__CINT__) || defined(__CLING__) || defined(__ROOTCLING__)

#pragma link off all classes;
#pragma link off all globals;
#pragma link off all functions;

#pragma link C++ class genElectron+;
#pragma link C++ class genPion+;
#pragma link C++ class vector<genElectron>+;
#pragma link C++ class vector<genPion>+;

#endif
