#ifndef GETRESIDUALJETCORR_H
#define GETRESIDUALJETCORR_H

#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"

#include <string>
#include <iostream>

const std::string residualStr = "RESIDUALCORR";

const Int_t nResCentBins = 4;
const Int_t resCentBins[nResCentBins+1] = {200, 100, 60, 20, 0};

TFile* residualJetCorrFile_p;
TH1F* residualJetCorrHists_p[nResCentBins];
TF1* residualJetCorrFuncs_p[nResCentBins];


//run once before any loops
Bool_t initGetResidualJetCorr(const std::string inResCorrFileName)
{
  if(inResCorrFileName.find(residualStr) == std::string::npos){
    std::cout << "Input residual file \'" << inResCorrFileName << "\' not proper format. Give \'*RESIDUALCORR*\'. Return false." << std::endl;
    return false;
  }

  residualJetCorrFile_p = new TFile(inResCorrFileName.c_str(), "READ");

  if(residualJetCorrFile_p->IsZombie()){
    std::cout << "Input residual file \'" << inResCorrFileName << "\' is zombie. Return false" << std::endl;
    return false;
  }

  for(Int_t iter = 0; iter < nResCentBins; iter++){
    std::cout << Form("resCorr_cent%dto%d_h", resCentBins[iter+1]/2, resCentBins[iter]/2) << std::endl;
    residualJetCorrHists_p[iter] = (TH1F*)residualJetCorrFile_p->Get(Form("resCorr_cent%dto%d_h", resCentBins[iter+1]/2, resCentBins[iter]/2));
    
    residualJetCorrFuncs_p[iter] = (TF1*)residualJetCorrHists_p[iter]->GetFunction("f1_p");
  }

  return true;
}


Float_t getResCorrJetPt(Float_t jtPt, Int_t hiBin)
{
  if(jtPt < 0){
    std::cout << "Input jtpt < 0. return -999" << std::endl;
    return -999;
  }
  else if(hiBin < 0 || hiBin >= 200){
    std::cout << "Input hibin < 0 or >= 200. Return -999" << std::endl;
    return -999;
  }

  Int_t centPos = -1;

  for(Int_t iter = 0; iter < nResCentBins; iter++){
    if(hiBin >= resCentBins[iter+1] && hiBin < resCentBins[iter]){
      centPos = iter;
      break;
    }
  }
  
  Float_t corrJtPt = jtPt/residualJetCorrFuncs_p[centPos]->Eval(jtPt);

  return corrJtPt;
}

#endif
