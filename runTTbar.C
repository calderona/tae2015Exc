#include "AnalysisTTbar.C"


void runTTbar(TString sample, TString filename)
{
  TFile* file = TFile::Open(filename);

  TTree* inputTree = (TTree*)file->Get("Tree");

  AnalysisTTbar att(inputTree);

  att.Loop(sample);
}


# ifndef __CINT__
int main()
{
  /*
  runTTbar("Data",           "../nTuples/pseudoData.root");

  runTTbar("TTbar_Pow",           "../nTuples/TTbar_Pow.root");
  runTTbar("tW",                  "../nTuples/tW.root");
  runTTbar("DY",                  "../nTuples/DY.root");
  runTTbar("WW",                  "../nTuples/WW.root");
  runTTbar("WJets",               "../nTuples/WJets.root");
  runTTbar("WZ",                  "../nTuples/WZ.root");
  runTTbar("ZZ",                  "../nTuples/ZZ.root");
  */

  runTTbar("TTbar_Mad",           "../nTuples/TTbar_Mad.root");
  runTTbar("TTbar_MCNLO",           "../nTuples/TTbar_MCNLO.root");
 

  return 0;
}
# endif
