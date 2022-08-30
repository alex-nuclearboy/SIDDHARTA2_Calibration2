/***********************************************
* SIDDHARTA-2 Experiment
* Aleksander K.                 2022-08
* Licensed under the Apache License, Version 2.0
***********************************************/

// Peak Finder

#include  <TROOT.h>
#include  <Riostream.h>
#include  <stdlib.h>
#include  <TTree.h>
#include  <TBranch.h>
#include  <TLeaf.h>
#include  <TMinuit.h>
#include  <TFitResult.h>
#include  <TFitResultPtr.h>
#include  "PeakFinder.h"

#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

void PeakFinder() {

  rebinFactor = 4;

  // Read ROOT file
  TString path = "/home/nuclearboy/SIDDHARTA2/rootfiles/SIDDHARTA2_xrays/output/20220602/";
  TSystemDirectory dir(path,path);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file; TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(".root")) {
        infile = new TFile(path+fname.Data(),"read");
        cout << fname.Data() << endl;
      }
    }
  }

  fname.Remove(0,7);
  fname.Remove(fname.Length()-5);

  // Take SDDs histograms
  for(Int_t iBUS=1;iBUS<nBUS;iBUS++) {
    for(Int_t iSDD=1;iSDD<nSDD;iSDD++) {
      hADC[iBUS][iSDD] = (TH1F*)infile->Get(Form("bus%d_sdd%d",iBUS,iSDD));
      if(!hADC[iBUS][iSDD]) continue;
      hADC[iBUS][iSDD]->Rebin(rebinFactor);
      hADCfit[iBUS][iSDD] = (TH1F*)hADC[iBUS][iSDD]->Clone(Form("hADCfit_bus%d_sdd%d",iBUS,iSDD));
      hADC[iBUS][iSDD]->GetXaxis()->SetTitle("ADC [channel]");
      hADC[iBUS][iSDD]->GetYaxis()->SetTitle("counts / 4 channel");
      hADC[iBUS][iSDD]->SetLineColor(1);
      hADC[iBUS][iSDD]->SetLineWidth(1);
    }
  }

  /////
  const Int_t nPFPeaksMAX = 13; // MAX number of peaks for the Peak Finder
  const Int_t nPFPeaks    = 3;  // DESIRED number of peaks

  // RANGE of ADC for the Peak Finder
  Float_t xminPeakFinder  = 1700;
  Float_t xmaxPeakFinder  = 3700;

  Int_t sigmaPeakFinder   = 20; // sigma for the Peak Finder (in units of ADC channels)
  sigmaPeakFinder = sigmaPeakFinder/rebinFactor;

  Float_t initThresholdPeakFinder = 0.01; // initial threshold parameter for the Peak Finder (std in TSpectrum is 0.05)
  Float_t initTolerance = 0.05;           // tolerance to check that the peak assumption is correct (5%)

  Double_t *xpeaks;
  Double_t *ypeaks;

  Float_t xadc[nPFPeaksMAX] = {};
  Float_t yadc[nPFPeaksMAX] = {};
  Float_t PFPeakE[nPFPeaks] = {};
  TString* PFPeakName[nPFPeaks];

  // Assumption for the Peak Finder peaks:  TiKa - CuKa - CuKb
  PFPeakE[0] = eTiKa;
  PFPeakName[0] = new TString("TiKa");
  PFPeakE[1] = eCuKa;
  PFPeakName[1] = new TString("CuKa");
  PFPeakE[2] = eCuKb;
  PFPeakName[2] = new TString("CuKb");

  // PRECALIBRATION configuration:
  // ============================
  Int_t bkgNparam   = 3;  //number of parameters for the background: p0+exp(p1+p2*x)
  Int_t gausNparam  = 3;  //number of parameters for the Gaussian function
  TString fstr ("");
  const Int_t nPeaksMAX = 99;
  Float_t peakE[nPeaksMAX] = {-1.};
  TString* peakName[nPeaksMAX];
  Float_t AddPeakE[nPeaksMAX] = {-1.};
  TString * AddPeakName[nPeaksMAX];
  Int_t nPeaks = 0;
  Int_t iInt[nPeaksMAX] = {-1}; //Gaussian parameters
  Int_t iMean[nPeaksMAX] = {-1};
  Int_t iSigma[nPeaksMAX] = {-1};
  Int_t ip0 = -1;  //background parameters
  Int_t ip1 = -1;
  Int_t ip2 = -1;
  Float_t xminPre[nPeaksMAX] = {-1.};
  Float_t xmaxPre[nPeaksMAX] = {-1.};

  // Firstly add PeakFinder peaks:
  for(Int_t i=0;i<nPFPeaks;i++) {AddPeakE[nPeaks] = PFPeakE[i]; AddPeakName[nPeaks] = PFPeakName[i]; nPeaks++;}
  //Then add the peaks as needed:
  AddPeakE[nPeaks] = eTiKb; AddPeakName[nPeaks] = new TString("TiKb"); nPeaks++;
  //AddPeakE[nPeaks] = eCaKa; AddPeakName[nPeaks] = new TString("CaKa"); nPeaks++;
  //AddPeakE[nPeaks] = eMnKa; AddPeakName[nPeaks] = new TString("MnKa"); nPeaks++;
  //AddPeakE[nPeaks] = eFeKa; AddPeakName[nPeaks] = new TString("FeKa"); nPeaks++;
  //AddPeakE[nPeaks] = eFeKb; AddPeakName[nPeaks] = new TString("FeKb"); nPeaks++;
  //AddPeakE[nPeaks] = eBiLa; AddPeakName[nPeaks] = new TString("BiLa"); nPeaks++;

  // Ordered list of peaks for precalibration:
  Bool_t skip[nPeaksMAX] = {false};
  Int_t PeakOrder[nPeaksMAX] = {-1};
  Float_t themin = 999999.;
  Int_t imin = 0;
  for(Int_t i=0;i<nPeaks;i++) {
    //find the smallest in the AddPeaks list, write it in PeakOrder and remove it:
    themin = 999999.;
    for(Int_t j=0;j<nPeaks;j++) {
      if(skip[j]) continue;
      if(AddPeakE[j]<themin) {
        themin = AddPeakE[j];
        imin = j;
      }
    }
    PeakOrder[i] = imin;
    cout<<"PeakOder["<<i<<"] = "<<imin<<endl;
    skip[imin]=true;
  }
  // Reorder and print them:
  cout<<"Number of peaks "<<nPeaks<<endl;
  for(Int_t i=0;i<nPeaks;i++) {
    peakE[i] = AddPeakE[PeakOrder[i]];
    peakName[i] = AddPeakName[PeakOrder[i]];
    cout<<"Peak#"<<i<<". PeakName: "<<peakName[i]->Data()<<". Peak Energy = "<<peakE[i]<<endl;
  }

  // Set precalibration fit function: Gaussians + background
  fstr = "";
  for(Int_t i=0;i<nPeaks;i++) {
    iInt[i]   = gausNparam*i; //define order of the parameters
    iMean[i]  = gausNparam*i+1;
    iSigma[i] = gausNparam*i+2;
    if(i!=0) {fstr += "+";}
    fstr += Form("[%i]*exp(-0.5*pow(x-[%i],2)/pow([%i],2))/sqrt(2.*%f)/[%i]",iInt[i],iMean[i],iSigma[i],PI,iSigma[i]);
  }
  ip0 = gausNparam*nPeaks;
  ip1 = gausNparam*nPeaks+1;
  ip2 = gausNparam*nPeaks+2;
  fstr += Form("+exp([%i]+[%i]*x)+[%i]",ip0,ip1,ip2);
  //fstr += Form("+[%i]+x*[%i]",ip0,ip1);
  //cout<<" PRECALIBRATION fit function is: "<<fstr<<endl;

  // PEAK FINDER:
  //============
  for(Int_t iBUS=1;iBUS<nBUS;iBUS++) {
    for(Int_t iSDD=1;iSDD<nSDD;iSDD++) {
      if(!hADC[iBUS][iSDD]) continue;
      TH1F* histoPF;
      histoPF = (TH1F*)hADC[iBUS][iSDD]->Clone("histoPF");
      Int_t minStats = 10000;   //min statistics for calibration
      if(histoPF->GetEntries()<minStats) continue;
      cout<<endl<<"--- PEAK FINDER: BUS#"<<iBUS<<" SDD#"<<iSDD<<". Statistics: "<<histoPF->GetEntries()<<" ---"<<endl<<endl;
      histoPF->SetAxisRange(xminPeakFinder,xmaxPeakFinder,"X");

      // Use TSpectrum to find the peak candidates
      Int_t nfound = 0;
      Float_t thresholdPeakFinder = initThresholdPeakFinder;
      Int_t nPFtries = 0;
      TSpectrum *spectrum = new TSpectrum(nPFPeaksMAX);
      while(nfound<nPFPeaks && nPFtries<15) {
        nfound = spectrum->Search(histoPF,sigmaPeakFinder,"",thresholdPeakFinder);
        printf("Found %d candidate peaks to fit:\n",nfound);
        thresholdPeakFinder = thresholdPeakFinder*.1;   //initial = 0.01. It changes until it finds peaks
        nPFtries++;
      }
      if(nPFtries>=15) {
        cout<<"PEAK FINDER DOES NOT WORK, ==> CONTINUE"<<endl;
        continue;
      }

      xpeaks = spectrum->GetPositionX();  // array with X-positions of the centroids found by TSpectrum
      ypeaks = spectrum->GetPositionY();  // array with Y-positions of the centroids found by TSpectrum

      // Reorder in adc counts and check compatibility with the maximum peaks assumption
      Float_t themin = 999999.;
      Int_t imin = 0;
      for(Int_t i=0;i<nfound;i++) {  // find the smallest, write it in xadc and remove it:
        themin = 999999.;
        for(Int_t j=0;j<nfound;j++) {
          if(xpeaks[j]<themin){
            themin = xpeaks[j];
            imin = j;
          }
        }
        xadc[i] = xpeaks[imin];
        yadc[i] = ypeaks[imin];
        xpeaks[imin] = 999999.;
      }
      // print them, now ordered:
      for(Int_t i=0;i<nfound;i++) {cout<<"ADC: "<<xadc[i]<<". Height: "<<yadc[i]<<endl;}

      // Check if the peaks found are compatible with the assumption:
      Float_t eDist10 = PFPeakE[1]-PFPeakE[0];
      Float_t eDist21 = PFPeakE[2]-PFPeakE[1];
      Float_t Erelation = eDist21/eDist10;
      //make a triad from all the peaks found:
      Float_t GPF  = 0.; //gain
      Float_t G0PF = 0.; //offset
      Bool_t TestPassed = false;
      Int_t ipeak0 = -1;
      Int_t ipeak1 = -1;
      Int_t ipeak2 = -1;
      for(Int_t i0=0;i0<nfound;i0++) {
        for(Int_t i1=i0+1;i1<nfound;i1++) {
          for(Int_t i2=i1+1;i2<nfound;i2++) {
            cout<<endl<<"-> Trying the triad: "<<i0<<" "<<i1<<" "<<i2<<endl;
            Float_t xDist10 = xadc[i1]-xadc[i0];
            Float_t xDist21 = xadc[i2]-xadc[i1];
            Float_t ADCrelation = xDist21/xDist10;
            cout<<"Checking assumption: Energy relation "<<Erelation<<" vs ADC relation "<<ADCrelation<<endl;
            //cout<<endl;

            //Define tolerance parameter
            Float_t tol = initTolerance; // 5%
            Bool_t TolerancePass = true;
            if(fabs(1.-(Erelation/ADCrelation))>tol) TolerancePass = false;
            if(!TolerancePass) cout<<"The test was not passed: ERROR!!!"<<endl;

            // Get the Peak Finder calibration offset GPFo and slope GPFs
            Float_t Dadc = xadc[i0]-xadc[i1];
            Float_t De   = PFPeakE[0]-PFPeakE[1];
            GPF  = De/Dadc;
            G0PF = -1.*xadc[i0]*GPF + PFPeakE[0];
            cout<<"Offset G0PF = "<<G0PF<<" and gain GPF = "<<GPF<<endl;

            //define an acceptable gain and offset; check if conditions (tolerance, G and G0) are met
            Float_t mingoodG   = 2.9;    Float_t maxgoodG  = 3.9;
            Float_t mingoodG0  = -3000;  Float_t maxgoodG0 = -1000;
            if(GPF<maxgoodG && GPF>mingoodG && G0PF<maxgoodG0 && G0PF>mingoodG0 && TolerancePass) {
              cout<<" -- TEST PASSED! --"<<endl;
              cout<<endl;
              TestPassed = true;

              ipeak0 = i0;
              ipeak1 = i1;
              ipeak2 = i2;
            }
            if(TestPassed) break;
          }
       if(TestPassed) break;
       }
       if(TestPassed) break;
     }

     // Find also the highest height among the selected ones
     Float_t highestPeak = 0.;
     if(ipeak0>-1) {
       if(yadc[ipeak0]>highestPeak) highestPeak = yadc[ipeak0];
       if(yadc[ipeak1]>highestPeak) highestPeak = yadc[ipeak1];
       if(yadc[ipeak2]>highestPeak) highestPeak = yadc[ipeak2];
     }

     // FIT spectrum:
     //=============
     Float_t sigmaADCguess = 20.; //start with an initial peak resolution in ADC of 20
     for(Int_t i=0;i<nPeaks;i++) { //fit limits
       xminPre[i] = (peakE[i]-G0PF)/GPF - 2.5*sigmaADCguess;
       xmaxPre[i] = (peakE[i]-G0PF)/GPF + 3.*sigmaADCguess;
     }
     //put the errors to infinite outside the boundaries so this area does not fit:
     for(Int_t ibin=1;ibin<=histoPF->GetNbinsX();ibin++) {
       Bool_t keep = false;
       for(Int_t i=0;i<nPeaks;i++) {
         if(histoPF->GetBinCenter(ibin)>xminPre[i] && histoPF->GetBinCenter(ibin)<xmaxPre[i]) {keep=true;}
       }
       if(!keep) histoPF->SetBinError(ibin,0.);
     }

     fPre[iBUS][iSDD] = new TF1(Form("fPre_bus%d_sdd%d",iBUS,iSDD),fstr,xminPre[0],xmaxPre[nPeaks-1]); //fit function
     for(Int_t i=0;i<nPeaks;i++) {  //set initial parameters for Gaussian function
       fPre[iBUS][iSDD]->SetParameter(iInt[i],highestPeak*sigmaADCguess);
       fPre[iBUS][iSDD]->SetParName(iInt[i],Form("gauss%i_Hei",i));
       //cout<<"fPre->SetParameter("<<iInt[i]<<","<<yadc[0]*sigmaADCguess<<");"<<endl;
       fPre[iBUS][iSDD]->SetParameter(iMean[i],(peakE[i]-G0PF)/GPF);
       fPre[iBUS][iSDD]->SetParName(iMean[i],Form("gauss%i_Mea",i));
       //cout<<"fPre->SetParameter("<<iMean[i]<<","<<(PeakE[i]-G0PF)/GPF<<");"<<endl;
       fPre[iBUS][iSDD]->SetParameter(iSigma[i],sigmaADCguess);
       fPre[iBUS][iSDD]->SetParName(iSigma[i],Form("gauss%i_Sig",i));
       //cout<<"fPre->SetParameter("<<iSigma[i]<<","<<sigmaADCguess<<");"<<endl;
     }
     //set parameters for background function
     fBkgnd[iBUS][iSDD] = new TF1(Form("fBkgnd_bus%d_sdd%d",iBUS,iSDD),"expo(0)+pol0(2)",1500,4500);
     hADCfit[iBUS][iSDD]->Fit(Form("fBkgnd_bus%d_sdd%d",iBUS,iSDD),"","",1500,4500);
     fPre[iBUS][iSDD]->SetParameter(ip0,fBkgnd[iBUS][iSDD]->GetParameter(0));
     fPre[iBUS][iSDD]->SetParameter(ip1,fBkgnd[iBUS][iSDD]->GetParameter(1));
     fPre[iBUS][iSDD]->SetParameter(ip2,fBkgnd[iBUS][iSDD]->GetParameter(2));
/*
     fPre[iBUS][iSDD]->SetParameter(ip0,histoPF->GetEntries()/histoPF->GetNbinsX()/20);
     fPre[iBUS][iSDD]->SetParLimits(ip0,0.,histoPF->GetEntries()/histoPF->GetNbinsX());
     //cout<<"fPre->SetParameter("<<ip0<<","<<thehisto->GetEntries()/thehisto->GetNbinsX()/10<<");"<<endl;
     fPre[iBUS][iSDD]->SetParameter(ip1,0.);
     fPre[iBUS][iSDD]->SetParLimits(ip1,-0.01,0.01);
     //fPre[iBUS][iSDD]->SetParameter(ip2,1.);
     //cout<<"fPre->SetParameter("<<ip1<<",0.);"<<endl;
*/
     //histoPF->Fit(fPre[iBUS][iSDD],"","0R",xminPre[0],xmaxPre[nPeaks-1]);
     histoPF->Fit(fPre[iBUS][iSDD],"","R",xminPre[0],xmaxPre[nPeaks-1]);

     //put back original errors:
     for(Int_t ibin=1;ibin<=histoPF->GetNbinsX();ibin++) {
       histoPF->SetBinError(ibin,hADC[iBUS][iSDD]->GetBinError(ibin));
     }

     ofstream  peaksFile;

     peaksFile.open("output/files/peaks_"+fname+".dat", ios::app);
     peaksFile<<Form("%d  %d  %g  %g  %g  %g",iBUS,iSDD,fPre[iBUS][iSDD]->GetParameter(1),fPre[iBUS][iSDD]->GetParError(1),fPre[iBUS][iSDD]->GetParameter(7),fPre[iBUS][iSDD]->GetParError(7))<<endl;
     peaksFile.close();

     histoPF->Delete();

   }
 }

/////////////////////////////////////////HISTOGRAMS/////////////////////////////////////////

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPalette(1,0);
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetTextFont(42);

  TPad *pad[nBUS][nSDD];
  for(Int_t i=1;i<7;i++) {
    for(Int_t j=0;j<4;j++) {
      myCanvas[i][j] = new TCanvas;
      myCanvas[i][j]->Divide(4,4);  //column, line
      for(Int_t k=1;k<17;k++) {
        pad[i][j] = (TPad*)myCanvas[i][j]->cd(k);
        pad[i][j]->SetLogy();
        pad[i][j]->SetGrid();
        if(!hADC[i][k+16*j]) continue;
        hADC[i][k+16*j]->SetTitle(Form("BUS: %d, SDD: %d",i,k+16*j));
        hADC[i][k+16*j]->SetAxisRange(1500,4500,"X");
        hADC[i][k+16*j]->Draw();
        if(!fPre[i][k+16*j]) continue;
        fPre[i][k+16*j]->Draw("same");
        fBkgnd[i][k+16*j]->SetLineWidth(1);
        fBkgnd[i][k+16*j]->SetLineColor(7);
        fBkgnd[i][k+16*j]->Draw("same");
      }
      myCanvas[i][j]->Print(Form("output/plots/20220602/hADC_Ylog_bus%d_sdd%d_%d.png",i,j*16,16*(j+1)), "png");
    }
  }

}
