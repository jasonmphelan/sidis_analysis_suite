#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TRatioPlot.h"

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

void makeCanvas(TString inFileName_1,TString outFileName, TString var, TString varTit);

void makeComparisonCanvas(TString inFileName_1, TString inFileName_2,TString outFileName, TString var, TString varTit,  TString pionType_1, TString pionType_2, TString tit_1, TString tit_2);

void makeComparisonCanvas(TH1F* h1, TH1F* h2,TString outFileName, TString var, TString varTit,  TString pionType_1, TString pionType_2, TString tit_1, TString tit_2){

void makeComparisonCanvas(TString inFileName_1, TString inFileName_2,TString outFileName, TString var, TString varTit,  TString pionType, TString tit_1, TString tit_2);

void makeAllComparisonChargeCanvas(TString inFileName_1, TString inFileName_2, int nBinsQ2, int nBinsXb, TString outFileBase);

void makeAllComparisonCanvas(TString inFileName_1, TString inFileName_2, int nBinsQ2, int nBinsXb, TString outFileBase, TString tit1, TString tit2);

void make2DCanvas(TString inFileName, TString outFileName, TString histName, TString xAxis, TString yAxis);

void make2DCanvasLog(TString inFileName, TString outFileName, TString histName, TString xAxis, TString yAxis);

void loopThrough2DCanvas(TString inFileName, TString outFileName);

void makeAllDetectorPlots();
