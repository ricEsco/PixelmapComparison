#include <string>
#include <fstream>
#include <algorithm> // for std::fill_n and std::copy_n
#include <iterator>  // for std::__to_address

using namespace std;

void xtalk(){

    // ---- ---- ---- ---- Declaration of variables ---- ---- ---- ----
    string module = "RH00026";
    string injtype1_run = "Run000059"; // Run # of Pixelalive with injection type 1
    string injtype5_run = "Run000060"; // Run # of Pixelalive with injection type 5
    string injtype6_run = "Run000061"; // Run # of Pixelalive with injection type 6
    float alive_eff = 0.9;             // Threshold (as a fraction) for a pixel to be considered alive
    float coupled_eff = 0.5;           // Threshold (as a fraction) for a pixel to be considered coupled
    float uncoupled_eff = 0.3;         // Threshold (as a fraction) for a pixel to be considered uncoupled
   
   // ---- ---- ---- ---- Retrieving files and histograms ---- ---- ---- ----
    std::string filename1 = "input/" + injtype1_run + "_PixelAlive.root";
    std::string filename5 = "input/" + injtype5_run + "_PixelAlive.root";
    std::string filename6 = "input/" + injtype6_run + "_PixelAlive.root";

    TFile *f_injtype1 = new TFile(filename1.c_str());
    TFile *f_injtype5 = new TFile(filename5.c_str());
    TFile *f_injtype6 = new TFile(filename6.c_str());

    string baseDir = "Detector/Board_0/OpticalGroup_0/Hybrid_0/Chip_";
    string shortBaseDir = "D_B(0)_O(0)_H(0)_";

    for (int ch = 12; ch >= 12; ch--){

    string chip = to_string(ch);
    string plotDir = "plots/";
    
    TCanvas *c0_pixelalive1 = (TCanvas*) f_injtype1->Get((baseDir+chip+"/"+shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
    TH2F *h_pixelalive1 = (TH2F*)c0_pixelalive1->GetPrimitive((shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
    TCanvas *c0_pixelalive5 = (TCanvas*) f_injtype5->Get((baseDir+chip+"/"+shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
    TH2F *h_pixelalive5 = (TH2F*)c0_pixelalive5->GetPrimitive((shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
    TCanvas *c0_pixelalive6 = (TCanvas*) f_injtype6->Get((baseDir+chip+"/"+shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
    TH2F *h_pixelalive6 = (TH2F*)c0_pixelalive6->GetPrimitive((shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
    
    int nColumns = h_pixelalive5->GetXaxis()->GetNbins();
    int nRows = h_pixelalive5->GetYaxis()->GetNbins();
    
    cout << "nRows: " << nRows << "  nColumns: " << nColumns << endl;
    
    

    // ---- ---- ---- ---- analysis ---- ---- ---- ----
    
    vector<int> dead_row,dead_col,suspicious_row,suspicious_col,confirmed_row,confirmed_col;
    TH2F *h_suspicious2D = (TH2F*)h_pixelalive5->Clone("h_suspicious2D");
    h_suspicious2D->SetTitle(("suspicious channels of chip "+chip).c_str());
    TH2F *h_confirmed2D = (TH2F*)h_pixelalive5->Clone("h_confirmed2D");
    h_confirmed2D->SetTitle(("confirmed disconnected channels of chip "+chip).c_str());
    
    
    for (int i=0; i<nRows; i++){
        for (int j=0; j<nColumns; j++){
            bool detectable = true;
            if (h_pixelalive1->GetBinContent(j+1,i+1)<alive_eff){
                dead_row.push_back(i);
                dead_col.push_back(j);
                detectable = false;
            }
            if (detectable){
                if ((i==0 && j%2 == 0) || (i==nRows-1 && j%2 == 1)){
                    detectable = false;
                }
            }
            if (detectable && h_pixelalive1->GetBinContent(j+1,i+1)>alive_eff && h_pixelalive5->GetBinContent(j+1,i+1)<coupled_eff && h_pixelalive6->GetBinContent(j+1,i+1)<uncoupled_eff){
                confirmed_row.push_back(i);
                confirmed_col.push_back(j);
                h_confirmed2D->SetBinContent(j+1,i+1,1);
            }
            else {
                h_confirmed2D->SetBinContent(j+1,i+1,0);
            }
        }
    }
    
    cout << "chip " << ch << ":" << endl;
    cout << "    alive_eff, coupled_eff, uncoupled_eff: " << alive_eff << ", " << coupled_eff << ", " << uncoupled_eff << endl;
    cout << "    dead:        " << dead_row.size() << endl;
    cout << "    suspicious:  " << suspicious_row.size() << endl;
    cout << "    confirmed:   " << confirmed_row.size() << endl;
    
    
    // ---- ---- ---- ---- plotting ---- ---- ---- ----
    
    gStyle->SetOptStat(0);
    // TCanvas arguments are: name, title, x, y, width, height    
    TCanvas *c_confirmed2D = new TCanvas(("c_confirmed2D_"+chip).c_str(),("Confirmed disconnected channels of chip "+chip).c_str(),800,768); // Necessary plot
    h_confirmed2D->Draw("colz");
    //c_confirmed2D->SaveAs((plotDir+module+"_confirmed2D_chip"+chip+".png").c_str());

    //Create root file out of missing bumps histogram
    TFile out_file(("output/xtalk_m-"+module+"_c-"+chip+".root").c_str(),"RECREATE");
    h_confirmed2D->Write();
    out_file.Close();
     
     
     
     
   } //loop on chips
  
}
