#include "Plotter.h"

# ifndef __CINT__
int main()
{
    // Create a plotter
    Plotter plotter("rootfiles/");
    
    // Specify the data histogram, along with a label to be used in the legend
    plotter.AddDataHist("Data", " Data ");
    
    // Specify MC histograms, along with their colours and legend labels
    plotter.AddMCHist("TTbar_Pow", kOrange + 1, " t#bar{t} ");
    plotter.AddMCHist("tW", kRed + 1, " tW ");
    plotter.AddMCHist("WJets", kGray + 1, " W+Jets");
    plotter.AddMCHist("ZZ", kBlue, " ZZ ");  
    plotter.AddMCHist("WZ", kBlue, " WZ ");
    plotter.AddMCHist("DY", kGreen+2, " Z/#gamma* ");
    plotter.AddMCHist("WW", kCyan, " WW ");
    
    // Produce files with pictures
    plotter.Plot("h_met", "em", "Met (GeV/c)",1);
    
    
    return 0;
}
# endif
