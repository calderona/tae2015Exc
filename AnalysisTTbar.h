#ifndef AnalysisTTbar_h
#define AnalysisTTbar_h

#include "AnalysisBase.h"

#include <TH1.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream>
#include <iostream>


class AnalysisTTbar : public AnalysisBase
{
 public :

 AnalysisTTbar(TTree *tree=0) : AnalysisBase(tree) {}

  
  void  FillHistograms   (int     ichannel,
			  int     icut);
  void  Loop             (TString sample);
  void  Summary          (TString precision,
			  TString title);
};

#endif
