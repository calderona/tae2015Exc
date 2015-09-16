////////////////////////////////////////////////////
// Universidad de Cantabria                       // 
// TAE 2015, Benasque, Spain                      //
// Authors: Jonatan Piedra, Alicia Calderon       //
// Septiembre 2015 - Version v1.0                 // 
////////////////////////////////////////////////////


#define AnalysisTTbar_cxx
#include "AnalysisTTbar.h"


//==============================================================================
//
// Constants, enums and structs
//
//==============================================================================
const int  verbosity    = 1;


const float ELECTRON_MASS =  0.000511;  // [GeV]
const float MUON_MASS     =  0.106;     // [GeV]
const float TAU_MASS      =  1.777;     // [GeV]
const float Z_MASS        = 91.188;     // [GeV]

using namespace std; 

const int nchannel = 4;

enum {
  ee,
  mm,
  em,
  all
};

const TString schannel[nchannel] = {
  "mm",
  "ee",
  "em",
  "all"
};


const int ncut = 3;

enum {
  Exactly2Leptons,
  Has2Jets,
  OneBJet,
};

const TString scut[ncut] = {
  "Exactly2Leptons",
  "Has2Jets",
  "OneBJet",
};


//==============================================================================
//
// Data members
//
//==============================================================================

float               luminosity;
Float_t             event_weight;
int                 selected_channel;
ofstream            txt_summary;
ofstream            txt_events_ee;
ofstream            txt_events_mm;
ofstream            txt_events_em;
TFile*              root_output;
TString             filename;

TH1F*               h_counter_raw[nchannel][ncut];
TH1F*               h_counter_lum[nchannel][ncut];

TH1F*               h_met[nchannel][ncut];

//------------------------------------------------------------------------------
// Loop
//------------------------------------------------------------------------------
void AnalysisTTbar::Loop(TString sample)
{
  TH1::SetDefaultSumw2();

  luminosity = 1.0;  // fb-1 luminosity to be used 
  filename   = sample;

  gSystem->mkdir("rootfiles", kTRUE);
  gSystem->mkdir("txt",       kTRUE);

  root_output = new TFile("rootfiles/" + filename + ".root", "recreate");

 


  // Initialize histograms
  //----------------------------------------------------------------------------

  for (int i=0; i<nchannel; i++) {
    for (int j=0; j<ncut; j++) {

      h_counter_raw[i][j] = new TH1F("h_counter_raw_" + schannel[i] + "_" + scut[j], "",   3, 0,   3);
      h_counter_lum[i][j] = new TH1F("h_counter_lum_" + schannel[i] + "_" + scut[j], "",   3, 0,   3);
      h_met[i][j] = new TH1F("h_met_" + schannel[i] + "_" + scut[j], "",  150, 0, 150 );

    }
  }


  // Loop over events
  //----------------------------------------------------------------------------
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();


  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);

    if (ientry < 0) break;

    fChain->GetEntry(jentry);

    if (verbosity == 1 && jentry%10000 == 0) std::cout << "." << std::flush;

    // to weight each MC sample to the expected number of events for a luminosity value. 
    event_weight = evtweight;
    event_weight *= luminosity;

    if (verbosity > 0 && jentry==0)
      {
      printf("\n");
      printf(" [%s]\n", sample.Data());
      printf(" Will run on %lld events\n", nentries);
      printf(" event_weight: %f\n", event_weight);
      printf("\n");
      }


    // Classify the channels
    //--------------------------------------------------------------------------
   
    selected_channel = -1;
 
    if (channel == 0)  selected_channel = mm; 
    if (channel == 1)  selected_channel = ee; 
    if (channel == 2 || channel == 3)  selected_channel = em; 
   
 
    // 1. Loop over leptons
    //--------------------------------------------------------------------------
    // The leptons should have sufficient transverse momentum and should not be too forward
    
    if (pt2 <= 20) continue; 
    if (pt1 <= 20) continue;
    if (ch1*ch2 > 0) continue;

 
    // Fill histograms
    //--------------------------------------------------------------------------
    FillHistograms(selected_channel, Exactly2Leptons);


    // 2. Loop over jets
    //--------------------------------------------------------------------------
    // Require that there are at least four central jets with pt > 30 GeV

    unsigned nGoodJets = 0;

    Int_t nJet1 = 0, nJet2 = 0, nJet3 = 0, nJet4 = 0, nJet5 = 0; 
  
    if ( jetpt1 >= 30  && fabs(jeteta1 < 2.4) ) nJet1++;
    if ( jetpt2 >= 30  && fabs(jeteta2 < 2.4) ) nJet2++;
    if ( jetpt3 >= 30  && fabs(jeteta3 < 2.4) ) nJet3++;
    if ( jetpt4 >= 30  && fabs(jeteta4 < 2.4) ) nJet4++;

    nGoodJets = nJet1 + nJet2 + nJet3 + nJet4; 
    
    // how is the njet spectrum ? --> try to plot the njet distribution. 

    if ( nGoodJets < 2 )  continue;


    // Fill histograms
    //--------------------------------------------------------------------------
    FillHistograms(selected_channel, Has2Jets);


    // 3. Require b-tagging: at least one Good jet passing the b-jet tagger
    //--------------------------------------------------------------------------

    unsigned nBJets = 0;

    Int_t nBJet1 = 0, nBJet2 = 0, nBJet3 = 0, nBJet4 = 0, nBJet5 = 0;

    if ( nJet1 > 0 && jetbtag1 > 2.1 ) nBJet1++;
    if ( nJet2 > 0 && jetbtag2 > 2.1 ) nBJet2++;
    if ( nJet3 > 0 && jetbtag3 > 2.1 ) nBJet3++;
    if ( nJet4 > 0 && jetbtag4 > 2.1 ) nBJet4++;

    nBJets = nBJet1 + nBJet2 + nBJet3 + nBJet4; 

    if (nBJets < 1) continue;

  
    FillHistograms(selected_channel, OneBJet);

  
  }


  if (verbosity >  0) printf("\n");
  if (verbosity == 1) printf("\n");



  //----------------------------------------------------------------------------
  // Summary
  //----------------------------------------------------------------------------
  txt_summary.open("txt/" + filename + ".txt");

  txt_summary << Form("\n%20s results with %.0f fb\n", filename.Data(), luminosity);

  Summary("11.0", "raw yields");
  Summary("11.2", "predicted yields");

  txt_summary.close();

  root_output->cd();

  root_output->Write("", TObject::kOverwrite);

  root_output->Close();
}



//------------------------------------------------------------------------------
// FillHistograms
//------------------------------------------------------------------------------
void AnalysisTTbar::FillHistograms(int ichannel, int icut)
{
  h_counter_raw[ichannel][icut]->Fill(1);
  h_counter_lum[ichannel][icut]->Fill(1, event_weight);
  h_met[ichannel][icut]->Fill(metpt, event_weight);
			      //h_invMass2Lep[ichannel][icut]->Fill(mll, event_weight);

  if (ichannel != all) FillHistograms(all, icut);
}


//------------------------------------------------------------------------------
// Summary
//------------------------------------------------------------------------------
void AnalysisTTbar::Summary(TString precision, TString title)
{
  txt_summary << Form("\n%20s", title.Data());

  for (int i=0; i<nchannel; i++) txt_summary << Form("%11s    %11s", schannel[i].Data(), " ");

  txt_summary << Form("\n---------------------\n");

  for (int i=0; i<ncut; i++) {
      
    txt_summary << Form("%20s", scut[i].Data());

    for (int j=0; j<nchannel; j++) {

      TH1F* h_counter = h_counter_raw[j][i];

      if (title.Contains("predicted")) h_counter = h_counter_lum[j][i];

      float yield = h_counter->Integral();
      float error = sqrt(h_counter->GetSumw2()->GetSum());

      TString yield_error = Form("%s%sf +- %s-%sf",
				 "%", precision.Data(),
				 "%", precision.Data());

      txt_summary << Form(yield_error.Data(), yield, error);
    }
      
    txt_summary << "\n";
  }
  
  txt_summary << "\n";
}
