#ifndef _PEAKFINDER_H
#define _PEAKFINDER_H

#include  <TFile.h>
#include  <TGraph.h>
#include  <TGraphErrors.h>
#include  <TClonesArray.h>
#include  <TF1.h>
#include  <TH1.h>
#include  <TH2.h>
#include  <TH1F.h>
#include  <TH2F.h>
#include  <TH1D.h>
#include  <TH2D.h>
#include  <TCanvas.h>
#include  <TLegend.h>
#include  <TPaveText.h>
#include  <TStyle.h>
#include  <TString.h>

TFile* infile;    //data file
TFile* outfile;   //output file
TString fname;    //name of the data file

static const Int_t nBUS = 7;   //last BUS #6
static const Int_t nSDD = 65;  //last SDD #64

static const Double_t PI = TMath::Pi();

Int_t rebinFactor;

//Energies of X-ray emmission lines
static const Double_t eTiKa = 4508.83;
static const Double_t eTiKb = 4931.81;
static const Double_t eCuKa = 8041.05;
static const Double_t eCuKb = 8905.29;
static const Double_t eMnKa = 5895.23;
static const Double_t eFeKa = 6399.47;
static const Double_t eFeKb = 7058.0;
static const Double_t eWKa = 8631.10;
static const Double_t eWKg = 11285.9;
static const Double_t eBrKa = 11908.26;
static const Double_t eBrKb = 13292.0;
static const Double_t eBiKa = 10828.1; //?
static const Double_t eBiKb = 13011.6;
static const Double_t eBiKg = 15247.7;
static const Double_t eSrKa = 14142.04;
static const Double_t eSrKb = 15836.0;
static const Double_t ePbKa = 10541.39; //?
static const Double_t ePbKb = 12616.14;
static const Double_t ePbKg = 14764.4;

TH1F  *hADC[nBUS][nSDD],*hADCfit[nBUS][nSDD];

TF1* fBkgnd[nBUS][nSDD];
TF1 *fPre[nBUS][nSDD];

TCanvas* myCanvas[7][65];

#endif
