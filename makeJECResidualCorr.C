#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TNamed.h"
#include "TMath.h"

#include <string>

const std::string resString1 = "RESIDUALHIST.root";
const std::string resString2 = "RESIDUAL2HIST.root";

const Int_t nCentBins = 4;
const Int_t centBins[nCentBins+1] = {100, 50, 30, 10, 0};


void FitResidual1(TH1F* hist_p, TH1F* histPerph_p = NULL)
{
  if(hist_p->GetEntries() == 0) return;

  TF1* f1_p = new TF1("f1_p", Form("[0] + [1]/TMath::Sqrt(x) + [2]/x"), 30, hist_p->GetXaxis()->GetXmax());

  Float_t xVal0 = hist_p->GetBinCenter(1);
  Float_t xVal1 = hist_p->GetBinCenter((Int_t)(hist_p->GetNbinsX()/2));
  Float_t xVal2 = hist_p->GetBinCenter(hist_p->GetNbinsX());

  Float_t yVal0 = hist_p->GetBinContent(1);
  Float_t yVal1 = hist_p->GetBinContent((Int_t)(hist_p->GetNbinsX()/2));
  Float_t yVal2 = hist_p->GetBinContent(hist_p->GetNbinsX());

  Float_t par0 = yVal2;
  Float_t par0Err = 0;
  Float_t par1 = (yVal1 - par0)*TMath::Sqrt(xVal1);
  Float_t par2 = (yVal0 - par0 - par1/TMath::Sqrt(xVal0))*xVal0;

  if(histPerph_p != NULL){
    par0 = histPerph_p->GetFunction("f1_p")->GetParameter(0);
    par0Err = histPerph_p->GetFunction("f1_p")->GetParError(0);
  }

  f1_p->SetParameter(0, par0);
  if(histPerph_p != NULL) f1_p->SetParLimits(0, par0-par0Err, par0+par0Err);
  else f1_p->SetParLimits(0, .9, 1.05);
  f1_p->SetParameter(1, par1);
  f1_p->SetParameter(2, par2);

  hist_p->Fit("f1_p", "Q M", "", 30, hist_p->GetXaxis()->GetXmax());

  delete f1_p;

  return;
}

void makeJECResidualCorr(const std::string inFileName)
{
  Bool_t isRes1 = false;
  if(inFileName.find(resString1.c_str()) != std::string::npos) isRes1 = true;

  Bool_t isRes2 = false;
  if(inFileName.find(resString2.c_str()) != std::string::npos) isRes2 = true;

  if(!isRes1 && !isRes2){
    std::cout << "Input file not residual. return" << std::endl; 
    return;
  } 

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TH1F* meanHist_p[nCentBins];
  TH1F* recoToGenHist_p[nCentBins];

  for(Int_t iter = 0; iter < nCentBins; iter++){
    meanHist_p[iter] = (TH1F*)inFile_p->Get(Form("akPu3PF/jtRecoOverGenVPt_Inc_FitMean_akPu3PF_cent%dto%d_h", centBins[iter+1], centBins[iter]));

    recoToGenHist_p[iter] = (TH1F*)inFile_p->Get(Form("akPu3PF/jtRecoOverGenVRecoPt_Inc_FitMean_akPu3PF_cent%dto%d_h", centBins[iter+1], centBins[iter]));

    std::cout << "iter: " << iter << std::endl;

    FitResidual1(meanHist_p[iter]);
  }

  std::string outFileName = inFileName;
  outFileName.replace(outFileName.find(resString1), resString1.size()+1, "RESIDUALCORR.root");
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "UPDATE");
  TNamed pathStr(Form("pathStr"), inFileName.c_str());
  pathStr.Write("", TObject::kOverwrite);
  
  for(Int_t iter = 0; iter < nCentBins; iter++){
    meanHist_p[iter]->Write(Form("resCorr_cent%dto%d_h", centBins[iter+1], centBins[iter]), TObject::kOverwrite);
    recoToGenHist_p[iter]->Write("", TObject::kOverwrite);
  }
  outFile_p->Close();
  delete outFile_p;

  inFile_p->Close();
  delete inFile_p;

  return;  
}
