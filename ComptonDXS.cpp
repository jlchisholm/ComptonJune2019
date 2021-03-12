///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  ComptonDXS.cpp                                                           //
//                                                                           //
//  Jenna Chisholm                                                           //
//  March 11, 2021                                                           //
//                                                                           //
//  Last Updated: March 11, 2021                                             //
//                                                                           //
//  Calculates differential cross sections at 150, 175, 200, 225, 250, 275,  //
//  300, 325, and 350MeV for Compton scattering.                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "physics.h"
#include "math.h"
using namespace std;



// ----------------------- TAGGER HISTOGRAM CONVERSION ----------------------------- //


// A method written by Phil which takes in a 1D histogram that is a function of Tagger
// channel and converts it to a function of photon energy based on a given
// Tagger conversion file for the beamtime
// Input: conversion file name, histogram to be converted, Output: converted histogram
TH1D* TaggChanToEnergy(TString sFile, TH1D *hCrossSec)
{
    TTree *tTagg = new TTree("tTagg","tTagg");            // TTree for conversion?
    tTagg->ReadFile(sFile, "Channel:Energy");             // read conversion file into tree
    int iTagg = tTagg->Draw("Channel:Energy","","goff");  // get the number of channels/energies?
    double *dTaggCh = tTagg->GetV1();                     // first column in file/tree is channel
    double *dTaggEn = tTagg->GetV2();                     // second column in file/tree is energy

    // I think we're checking the order of things
    bool bRevrCh = (dTaggCh[1] < dTaggCh [0]);   // check if channels are in descending order
    bool bRevrEn = (dTaggEn[1] < dTaggEn[0]);    // check if energies are in descending order
    bool bSumw2 = (hCrossSec->GetSumw2N() > 0);  // sum of weights(???) positive?

    double *dTaggBn;                  // for lowest value of tagger bins
    dTaggBn = new double[iTagg+1];    // set it to be an array of length # of bins

    // So we're going through and getting the lowest value/edge in each bin?
    for (int i=0; i<=iTagg; i++)
    {
        if (bRevrEn) // if energy is in descending order, we'll need to go backwards through indices (i.e. iTagg-i)
        {
            if (i==0) dTaggBn[0] = (dTaggEn[iTagg-1] + 0.5*(dTaggEn[iTagg-1] - dTaggEn[iTagg-2]));
            else if (i==iTagg) dTaggBn[iTagg] = (dTaggEn[0] + 0.5*(dTaggEn[0] - dTaggEn[1]));
            else dTaggBn[i] = 0.5*(dTaggEn[iTagg-1-i] + dTaggEn[iTagg-i]);
        }
        else
        {
          if (i==0) dTaggBn[0] = (dTaggEn[0] + 0.5*(dTaggEn[0] - dTaggEn[1]));
          else if (i==iTagg) dTaggBn[iTagg] = (dTaggEn[iTagg-1] + 0.5*(dTaggEn[iTagg-1] - dTaggEn[iTagg-2]));
          else dTaggBn[i] = 0.5*(dTaggEn[i-1] + dTaggEn[i]);  // low edge bin energy is the average of this energy and the previous energy
        }
    }

    // Plot
    TH1D *hTaggCS = new TH1D("hTaggCS", "Pi0 Production Total Cross Section from 150MeV to 350MeV", iTagg, dTaggBn);

    for (int i=0; i<iTagg; i++)
    {
        if (bRevrCh == bRevrEn) // if both channels and energies are both descending order or both are in ascending order
        {
            hTaggCS->SetBinContent(i+1, hCrossSec->GetBinContent(i+1));
            if (bSumw2) hTaggCS->SetBinError(i+1, hCrossSec->GetBinError(i+1));
        }
        else // basically if channels are ascending but energies are descending or vice versa
        {
            hTaggCS->SetBinContent(i+1, hCrossSec->GetBinContent(iTagg-i));
            if (bSumw2) hTaggCS->SetBinError(i+1, hCrossSec->GetBinError(iTagg-i));
        }
    }

    return hTaggCS;
}

// Another method written by Phil which takes in a 3D histogram that is a function of Tagger
// channel and converts it to a function of photon energy based on a given
// Tagger conversion file for the beamtime
// Input: conversion file name, histogram to be converted, axis that Tagger channel is on, Output: converted histogram
TH3* Tagg_Chan_To_Ener(TString sFile, TH3 *hChan, TString sAxis="X")
{
  TTree *tTagg = new TTree("tTagg", "tTagg");
  tTagg->ReadFile(sFile, "Channel:Energy");
  Int_t iTagg = tTagg->Draw("Channel:Energy", "", "goff");
  Double_t *dTaggCh = tTagg->GetV1();
  Double_t *dTaggEn = tTagg->GetV2();

  Bool_t bRevrCh = (dTaggCh[1] < dTaggCh[0]);
  Bool_t bRevrEn = (dTaggEn[1] < dTaggEn[0]);
  Bool_t bSumw2 = (hChan->GetSumw2N() > 0);

  Double_t *dTaggBn;
  dTaggBn = new Double_t[iTagg+1];

  for (Int_t i=0; i<=iTagg; i++)
  {
      if (bRevrEn)
      {
          if (i==0) dTaggBn[0] = (dTaggEn[iTagg-1] + 0.5*(dTaggEn[iTagg-1] - dTaggEn[iTagg-2]));
          else if (i==iTagg) dTaggBn[iTagg] = (dTaggEn[0] + 0.5*(dTaggEn[0] - dTaggEn[1]));
          else dTaggBn[i] = 0.5*(dTaggEn[iTagg-1-i] + dTaggEn[iTagg-i]);
      }
      else
      {
          if (i==0) dTaggBn[0] = (dTaggEn[0] + 0.5*(dTaggEn[0] - dTaggEn[1]));
          else if (i==iTagg) dTaggBn[iTagg] = (dTaggEn[iTagg-1] + 0.5*(dTaggEn[iTagg-1] - dTaggEn[iTagg-2]));
          else dTaggBn[i] = 0.5*(dTaggEn[i-1] + dTaggEn[i]);
      }
  }

  TH3 *hEner = (TH3*)hChan->Clone("hEner");
  hEner->Reset();

  TAxis *axis;
  if (sAxis == "x" || sAxis == "X") axis = hEner->GetXaxis();
  else if (sAxis == "y" || sAxis == "Y") axis = hEner->GetYaxis();
  else if (sAxis == "z" || sAxis == "Z") axis = hEner->GetZaxis();

  if (axis->GetNbins() != iTagg) cout << "WARNING - Different number of bins between config file and histogram!" << endl << "This could lead to unexpected results!" << endl;
  axis->Set(iTagg, dTaggBn);

  Int_t iOld, iNew;
  for (Int_t iX=1; iX<=hChan->GetNbinsX(); iX++)
  {
      for (Int_t iY=1; iY<=hChan->GetNbinsY(); iY++)
      {
          for (Int_t iZ=1; iZ<=hChan->GetNbinsZ(); iZ++)
          {
              iOld = hChan->GetBin(iX, iY, iZ);
              if (bRevrCh == bRevrEn) iNew = hEner->GetBin(iX, iY, iZ);
              else
              {
                  if (sAxis == "x" || sAxis == "X") iNew = hEner->GetBin(iTagg-iX+1, iY, iZ);
                  else if (sAxis == "y" || sAxis == "Y") iNew = hEner->GetBin(iX, iTagg-iY+1, iZ);
                  else if (sAxis == "z" || sAxis == "Z") iNew = hEner->GetBin(iX, iY, iTagg-iZ+1);
              }
              hEner->SetBinContent(iNew, hChan->GetBinContent(iOld));
              if (bSumw2) hEner->SetBinError(iNew, hChan->GetBinError(iOld));
          }
      }
  }

  return hEner;
}



// ------------------ TOTAL DETECTION EFFICIENCY VALUE AT CHANNEL/ENERGY ------------------------- //


// A function that gets the detection efficiency from the detection efficiency files
// Input: incident photon energy, output: total detection efficiency for that photon energy
TH1D* GetTotDetEff(int ke)
{
    // Open file
    TFile *edetFile = TFile::Open(Form("DetEff/DetEff_%d.root",ke));
    if(!edetFile->IsOpen())
    {
            cout << "Error: detection efficiency file could not be opened.";
            exit(-1);
    }

    // Get detection efficiency histogram
    TH1D* hEdet = (TH1D*)edetFile->Get("hDetEff");

    return hEdet;
}



// ---------------------- SCALING FACTOR FOR EMPTY TARGET SUBTRACTION --------------------------- //


// A function which draws the tagger counts histograms for full and empty target
// data and makes a histogram for the scaling ratio
// Input: photon energy, output: scaling factor for that channel
double GetScalingFactor(int ke, TFile *ftFile, TFile *etFile)
{
    // Get tagger count histograms
    TH1D *hFull = (TH1D*)ftFile->Get("He4Compton/h_ScalarCounts");
    TH1D *hEmpty = (TH1D*)etFile->Get("He4Compton/h_ScalarCounts");

    // Divide the full target data by empty target data to get a ratio
    TH1D *hRatio = (TH1D*) hFull->Clone("hRatio");
    hRatio->Divide(hEmpty);

    // Convert to a function of photon energy
    TH1D *hRatio_conv = TaggChanToEnergy("Tagger_Conversions_by_k_2019_06.txt",hRatio);

    // Get factor for this energy
    double factor = hRatio_conv->GetBinContent(ke);

    return factor;
}




//int GetCounts(int ke, int th, TH3* hFull_ME_3D_conv, TH3* hEmpty_ME_3D_conv, double scale)
int GetCounts(int ke)
{

    // --------------------------- Open Files --------------------------- //

    // Full target file
    TFile *ftFile = TFile::Open("Full_Target_Compton_mikes.root");
    if(!ftFile->IsOpen()){cout << "Error: full target file could not be opened."; exit(-1);}

    // Empty target file
    TFile *etFile = TFile::Open("Empty_Target_Compton_mikes.root");
    if(!etFile->IsOpen()){cout << "Error: empty target file could not be opened.";exit(-1);}

    // Tagging Efficiency File
    TFile *tagFile = TFile::Open("~/TaggingEfficienciesJune2019/TaggEff-MOELLER-20.root");
    if(!tagFile->IsOpen()){cout << "Error: tagging efficiency file could not be opened.";exit(-1);}


    // Get 3D, 1 particle event ME Histograms
    TH3D *hFull_ME_3D = (TH3D*)ftFile->Get("He4Compton/h3D_ME1");
    TH3D *hEmpty_ME_3D = (TH3D*)etFile->Get("He4Compton/h3D_ME1");

    // Convert the 3D histograms to be a function of photon energy instead of Tagger channel
    TH3 *hFull_ME_3D_conv = Tagg_Chan_To_Ener("Tagger_Conversions_by_k_2019_06.txt",hFull_ME_3D,"Z");
    TH3 *hEmpty_ME_3D_conv = Tagg_Chan_To_Ener("Tagger_Conversions_by_k_2019_06.txt",hEmpty_ME_3D,"Z");


    // Get scaling factor
    double scale = GetScalingFactor(ke, ftFile, etFile);



    // Project 3D histograms at input Tagger channel
    TH1D *hFull_ME_projx = hFull_ME_3D_conv->ProjectionX(Form("hME_projx_%dMeV",ke));  // just all of theta and ke for a moment now
    TH1D *hEmpty_ME_projx = hEmpty_ME_3D_conv->ProjectionX(Form("hEmpty_ME_projx_%dMeV",ke));

    // Subtract empty from full with scaling
    TH1D *hSubtracted = (TH1D*) hFull_ME_projx->Clone(Form("hSubtracted_%dMeV",ke));
    hSubtracted->Add(hEmpty_ME_projx, -scale);

    //Create a scaled version of the empty target histogram -- for plotting only
    TH1D *hEmpty_ME_scaled = (TH1D*) hEmpty_ME_projx->Clone(Form("hEmpty_ME_scaled_%dMeV",ke));
    hEmpty_ME_scaled->Scale(scale);


    // Create a canvas and draw the differential cross sections
    TCanvas *ctemp = new TCanvas("ctemp","",200,10,750,750);
    ctemp->cd();
    hFull_ME_projx->SetLineColor(1);
    hFull_ME_projx->Rebin(4);
    hFull_ME_projx->Draw("SAME HIST");
    hEmpty_ME_projx->SetLineColor(2);
    hEmpty_ME_projx->Rebin(4);
    hEmpty_ME_projx->Draw("SAME HIST");
    hEmpty_ME_scaled->SetLineColor(3);
    hEmpty_ME_scaled->Rebin(4);
    hEmpty_ME_scaled->Draw("SAME HIST");
    hSubtracted->SetLineColor(4);
    hSubtracted->Rebin(4);
    hSubtracted->Draw("SAME HIST");

    // Create a legend for the histogram
    TLegend *legend = new TLegend(0.1,0.75,0.3,0.9);
    legend->AddEntry(hFull_ME_projx, "Full Target Data","l");
    legend->AddEntry(hEmpty_ME_projx, "Empty Target Data","l");
    legend->AddEntry(hEmpty_ME_scaled, "Scaled Empty Target Data","l");
    legend->AddEntry(hSubtracted, "Subtracted Data","l");
    legend->Draw();


    int counts = 0;

    return counts;
}




//TH1D* GetDXS(int ke=200)
//{
//    return DXS;
//}




//void PlotDXS()
//{
//    // Calculate the nine differential cross sections
//    TH1D* hDiffCrossSec[9];
//    int ke=150;
//    for(int i=0; i<=8; i++)
//    {
//        hDiffCrossSec[i] = GetDXS(ke);
//        hDiffCrossSec[i]->SetLineColor(1);
//        ke=ke+25;
//    }

//    // Create a canvas and draw the differential cross sections
//    TCanvas *c1 = new TCanvas("c1","",200,10,750,750);
//    c1->Divide(3,3);
//    c1->cd(1);
//    hDiffCrossSec[0]->Draw("pfc");
//    c1->cd(2);
//    hDiffCrossSec[1]->Draw("pfc");
//    c1->cd(3);
//    hDiffCrossSec[2]->Draw("pfc");
//    c1->cd(4);
//    hDiffCrossSec[3]->Draw("pfc");
//    c1->cd(5);
//    hDiffCrossSec[4]->Draw("pfc");
//    c1->cd(6);
//    hDiffCrossSec[5]->Draw("pfc");
//    c1->cd(7);
//    hDiffCrossSec[6]->Draw("pfc");
//    c1->cd(8);
//    hDiffCrossSec[7]->Draw("pfc");
//    c1->cd(9);
//    hDiffCrossSec[8]->Draw("pfc");
//}
