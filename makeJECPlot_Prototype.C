#include "TFile.h"
#include "TDirectoryFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TKey.h"
#include "TLine.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TStyle.h"
#include "TF1.h"

#include <string>
#include <iostream>
#include <vector>

#include "include/doGlobalDebug.h"
#include "include/checkMakeDir.h"
#include "include/jecConfigParser.h"

#include "SuperCanvas/include/superCanvas.h"

Bool_t debugMode = false;

Bool_t plotTrue = false;

Bool_t addResCorrStr = true;

const Int_t nHistName = 3;
const std::string inHistName[nHistName] = {"RecoOverGen", "Eff", "DPhi"};
const std::string xAxisLabel[nHistName] = {"Reco./Gen.", "Eff.", "<#Delta#phi>"};

const Int_t nMeanRes = 2;
const std::string meanResStr[nMeanRes] = {"Mean", "Res"};
const std::string meanResStr2[nMeanRes] = {"#mu", "#sigma"};

const Int_t nQG = 3;
const std::string qgStr[nQG] = {"Inc", "Q", "G"};
const std::string qgStr2[nQG] = {"Inc.", "Quarks", "Gluons"};
const Int_t qgCol[nQG] = {kGray+1, kBlue, kRed};

const Int_t nPtEta = 5;
const std::string ptEtaStr[nPtEta] = {"VPt_EtaInc_", "VPt_Eta0p0to0p5_", "VPt_Eta0p5to1p0_", "VPt_Eta1p0to1p6_", "VEta_"};
const std::string ptEtaStr2[nPtEta] = {"Pt", "Pt", "Pt", "Pt", "Eta"};
const std::string ptEtaStr3[nPtEta] = {"p_{T}", "p_{T}", "p_{T}", "p_{T}", "#eta"};


void claverCanvasSaving(TCanvas* c, TString s,TString format="gif"){
  TDatime* date = new TDatime();
  c->SaveAs(Form("%s_%d.%s",s.Data(),date->GetDate(), format.Data()));
  return;
}

Float_t setMaxMinNice(Float_t inMaxMin, Bool_t isMax){
  for(Int_t iter = 0; iter < 5; iter++){
    if(inMaxMin*TMath::Power(10, iter+1) > 10){
      inMaxMin = inMaxMin*TMath::Power(10, iter+1);
      if(isMax) inMaxMin = std::ceil(inMaxMin);
      else inMaxMin = std::floor(inMaxMin);
      inMaxMin = inMaxMin/TMath::Power(10, iter+1);
      if(isMax) inMaxMin -= .01/TMath::Power(10, iter+1);
      else inMaxMin += .01/TMath::Power(10, iter+1);

      break;
    }
  }
  
  return inMaxMin;
}


int makeJECPlotMeanRes(const std::string inFileName, jecConfigParser config, const Int_t inHistNum, const Int_t meanResNum, const Int_t ptEtaNum)
{
  const Bool_t isPbPb = config.GetIsPbPb();

  const Int_t nCentBins = config.GetNCentBins();
  std::vector<unsigned int> centBins = config.GetCentBins();

  if(isPbPb){
    if((unsigned int)(nCentBins+1) != centBins.size()){
      std::cout << "Config nCentBins \'" << nCentBins << "\' != centBins.size \'" << centBins.size() << "\'. return 1" << std::endl;
      return 1;
    }
  }

  std::string centStrings[nCentBins];
  std::string centStrings2[nCentBins];

  for(Int_t iter = 0; iter < nCentBins; iter++){
    if(isPbPb){
      centStrings[iter] = "Cent" + std::to_string(centBins.at(iter+1)) + "to" + std::to_string(centBins.at(iter));
      centStrings2[iter] = std::to_string(centBins.at(iter+1)) + "-" + std::to_string(centBins.at(iter)) + "%";
    }
    else{
      centStrings[iter] = "PP";
      centStrings2[iter] = "PP";
    }
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nContents = inFile_p->GetListOfKeys()->GetEntries();

  Int_t nTH1Temp = 0;
  std::vector<std::string>* th1Names_p = new std::vector<std::string>;

  Int_t nDirTemp = 0;
  std::vector<std::string>* dirNames_p = new std::vector<std::string>;

  for(Int_t iter = 0; iter < nContents; iter++){
    TString className = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetClassName();
    TString name = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetName();

    if(className.Index("TDirectoryFile") >= 0 && name.Index("ak") >= 0){
      nDirTemp++;
      dirNames_p->push_back(name.Data());

      TDirectoryFile* tempDir_p = (TDirectoryFile*)inFile_p->Get(name);

      const Int_t nDirContents = tempDir_p->GetListOfKeys()->GetEntries();
      
      for(Int_t dirIter = 0; dirIter < nDirContents; dirIter++){
	TString name2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetName();
	TString className2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetClassName();

	if(className2.Index("TH1") >= 0){
	  if(name2.Index(Form("%s_", meanResStr[meanResNum].c_str())) >= 0){
	    if(name2.Index(ptEtaStr[ptEtaNum].c_str()) < 0) continue;

	    if(name2.Index("Res") >= 0 && name2.Index("_Q_") >= 0) continue;
	    if(name2.Index("Res") >= 0 && name2.Index("_G_") >= 0) continue;
	    if(name2.Index("Mean") >= 0 && name2.Index("_G_") >= 0) continue;

	    if(!strcmp("Fake", inHistName[inHistNum].c_str()) && (name2.Index("_Q_") >= 0 || name2.Index("_G_") >= 0)) continue;
	    if(!strcmp("Eff", inHistName[inHistNum].c_str()) && (name2.Index("_Q_") >= 0 || name2.Index("_G_") >= 0)) continue;

	    if(name2.Index("RecoOverGen") >= 0 && name2.Index("_Mean_") >= 0) continue;
	    if(name2.Index("RecoOverGen") >= 0 && name2.Index("_Res_") >= 0) continue;

	    if(name2.Index("RecoPtCut5") >= 0) continue;
	    if(name2.Index("RecoPtCut30") >= 0) continue;
	    //	    if(name2.Index("RecoPtCut10") >= 0) continue;

	    if(inHistName[inHistNum].find("DPhi") != std::string::npos && name2.Index("_Q_") >= 0) continue;
	    if(inHistName[inHistNum].find("DPhi") != std::string::npos && name2.Index("_G_") >= 0) continue;

	    if(name2.Index(inHistName[inHistNum].c_str()) >= 0){
	      nTH1Temp++;
	      th1Names_p->push_back(Form("%s/%s", name.Data(), name2.Data()));
	    }
	  }
	}
      }
    }
  }


  const Int_t nTH1 = nTH1Temp;
  const Int_t nDir = nDirTemp;
  TH1F* th1_p[nTH1];  
  
  if(debugMode) std::cout << __LINE__ << std::endl;

  std::vector<superCanvas*> th1Canv_p(nDir);
  //  TCanvas* th1Canv_p[nDir];

  Int_t nXPanel = nCentBins;
  Int_t nYPanel = 1;
  if(!isPbPb){
    nXPanel = 1;
    nYPanel = 1;
  }

  Int_t th1CanvTot[nDir][nXPanel][nYPanel];
  Int_t th1CanvCount[nDir][nXPanel][nYPanel];
  
  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t xIter = 0; xIter < nXPanel; xIter++){
      for(Int_t yIter = 0; yIter < nYPanel; yIter++){
        th1CanvTot[iter][xIter][yIter] = 0;
        th1CanvCount[iter][xIter][yIter] = 0;
      }
    }
  }

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
      if(th1Name.substr(0, th1Name.find("/")).size() != dirNames_p->at(dirIter).size()) continue;
      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    Int_t centPos = -1;
    if(isPbPb){
      for(Int_t centIter = 0; centIter < nCentBins; centIter++){
	std::size_t pos = th1Name.find(centStrings[centIter]);
        if(pos != std::string::npos){
          centPos = centIter;
          break;
        }
      }
    }
    else centPos = 0;

    th1CanvTot[dirPos][centPos][0]++;
  }  

  if(debugMode) std::cout << __LINE__ << std::endl;

  for(Int_t iter = 0; iter < nDir; iter++){
    th1Canv_p[iter] = new superCanvas(nXPanel, nYPanel, th1CanvTot[iter][0][0], 1000, dirNames_p->at(iter));

    th1Canv_p[iter]->SetGlobalMaxMin();
    
    if(meanResStr[meanResNum].find("Mean") != std::string::npos){
      if(inHistName[inHistNum].find("RecoOverGen") != std::string::npos){
	th1Canv_p[iter]->CapGlobalMaxMin(1.5, 0.8);
	th1Canv_p[iter]->UnderCapGlobalMaxMin(1.5, 0.8);
      }
      if(inHistName[inHistNum].find("Eff") != std::string::npos){
	th1Canv_p[iter]->CapGlobalMaxMin(1.05, 0.0);
	th1Canv_p[iter]->UnderCapGlobalMaxMin(1.05, 0.0);
      }
    }
    if(meanResStr[meanResNum].find("Res") != std::string::npos){
      if(inHistName[inHistNum].find("DPhi") != std::string::npos){
	th1Canv_p[iter]->CapGlobalMaxMin(0.2, 0.0);
	th1Canv_p[iter]->UnderCapGlobalMaxMin(0.2, 0.0);
      }
      else th1Canv_p[iter]->CapGlobalMaxMin(0.7, 0.0);
    }

    th1Canv_p[iter]->canv_p->cd();
    th1Canv_p[iter]->canv_p->Draw();

    for(Int_t xIter = 0; xIter < nXPanel; xIter++){
      for(Int_t yIter = 0; yIter < nYPanel; yIter++){
	th1Canv_p[iter]->canv_p->cd();
	th1Canv_p[iter]->pads_p[xIter][yIter]->Draw();
	th1Canv_p[iter]->pads_p[xIter][yIter]->cd();

      }
    }

    //    th1Canv_p[iter] = new TCanvas(Form("%s_%s%s%sc", dirNames_p->at(iter).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str()), Form("%s_%s%s%sc", dirNames_p->at(iter).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str()), nXPanel*300, nYPanel*325);
    //  th1Canv_p[iter]->Divide(nXPanel, nYPanel, 0.0, 0.0);
  }
  /*
  Bool_t isDrawn[nDir][nXPanel*nYPanel];
  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t iter2 = 0; iter2 < nXPanel*nYPanel; iter2++){
      isDrawn[iter][iter2] = false;
    }
    }*/

  if(debugMode) std::cout << __LINE__ << std::endl;

  //  Bool_t legAdded[2][nQG] = {{false}, {false}};

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  //  Float_t dirMax[nDir];
  //  Float_t dirMin[nDir];
  /*
  for(Int_t iter = 0; iter < nDir; iter++){
  dirMax[iter] = -1;
    dirMin[iter] = 100;
    }*/

  Float_t th1XMin = 30;

  if(debugMode) std::cout << __LINE__ << std::endl;

  TLegend* meanLeg_p;
  if(isPbPb){
    if(!strcmp("Eff", inHistName[inHistNum].c_str())) meanLeg_p = new TLegend(.70, .08, .90, .55);
    if(!strcmp("Res", meanResStr[meanResNum].c_str())) meanLeg_p = new TLegend(.70, .08, .90, .35);
    else meanLeg_p = new TLegend(.50, .68, .90, .95);
  }
  else if(!strcmp("Res", meanResStr[meanResNum].c_str())) meanLeg_p = new TLegend(.3, .38, .90, .5);
  else if(!strcmp("Eff", inHistName[inHistNum].c_str())) meanLeg_p = new TLegend(.6, .25, .95, .75);
  else meanLeg_p = new TLegend(.50, .58, .90, .85);

  meanLeg_p->SetBorderSize(0);
  meanLeg_p->SetFillColor(0);
  meanLeg_p->SetFillStyle(0);
  meanLeg_p->SetTextFont(43);
  meanLeg_p->SetTextSize(18);


  TLegend* meanLeg2_p;
  if(!strcmp("Eff", inHistName[inHistNum].c_str())){
    if(isPbPb) meanLeg2_p = new TLegend(.60, .25, .95, .75);
    else  meanLeg2_p = new TLegend(.50, .35, .95, .85);
  }
  else if(!strcmp("Res", meanResStr[meanResNum].c_str())) meanLeg2_p = new TLegend(.50, .38, .90, .65);
  else meanLeg2_p = new TLegend(.50, .68, .90, .95);
  meanLeg2_p->SetBorderSize(0);
  meanLeg2_p->SetFillColor(0);
  meanLeg2_p->SetFillStyle(0);
  meanLeg2_p->SetTextFont(43);
  meanLeg2_p->SetTextSize(18);

  if(debugMode) std::cout << __LINE__ << std::endl;

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    th1_p[iter] = (TH1F*)inFile_p->Get(th1Name.c_str());

    th1_p[iter]->SetMarkerSize(0.5);

    if(iter == 0 && strcmp("Pt", ptEtaStr2[ptEtaNum].c_str()) != 0) th1XMin = th1_p[iter]->GetBinCenter(1);

    Int_t xMinBin = th1_p[iter]->FindBin(th1XMin);

    if(iter == 0) th1XMin = th1_p[iter]->GetBinLowEdge(xMinBin);

    //insert new xmin
    th1XMin = th1_p[iter]->GetXaxis()->GetXmin();
    //    if(!strcmp("Eta", ptEtaStr2[ptEtaNum].c_str())) th1XMin += .1;
    th1_p[iter]->SetAxisRange(th1XMin, th1_p[iter]->GetXaxis()->GetXmax(), "X");

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
      
      //std::cout << dirNames_p->at(dirIter) << ", " << th1Name << std::endl;

      if(th1Name.substr(0, th1Name.find("/")).size() != dirNames_p->at(dirIter).size()) continue;

      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    th1Canv_p[dirPos]->canv_p->cd();
  
    Int_t centPos = -1;
    if(isPbPb){
      for(Int_t centIter = 0; centIter < nCentBins; centIter++){
	std::size_t pos = th1Name.find(centStrings[centIter]);
	if(pos != std::string::npos){
	  centPos = centIter;
	  break;
	}
      }
    }
    else centPos = 0;
    

    std::string legStr = "Incl.";
    if(inHistName[inHistNum].find("RecoOverGen") != std::string::npos){
      if(th1Name.find("_Q_") != std::string::npos) legStr = "Quark";
      else if(th1Name.find("_G_") != std::string::npos) legStr = "Gluon";
    }
    else if(inHistName[inHistNum].find("Eff") != std::string::npos){
      if(th1Name.find("RecoPtCut5") != std::string::npos) legStr = "Reco. p_{T} > 5";
      else if(th1Name.find("RecoPtCut10") != std::string::npos) legStr = "Reco. p_{T} > 10";
      else if(th1Name.find("RecoPtCut15") != std::string::npos) legStr = "Reco. p_{T} > 15";
      else if(th1Name.find("RecoPtCut20") != std::string::npos) legStr = "Reco. p_{T} > 20";
      else if(th1Name.find("RecoPtCut25") != std::string::npos) legStr = "Reco. p_{T} > 25";
      else if(th1Name.find("RecoPtCut30") != std::string::npos) legStr = "Reco. p_{T} > 30";
    }
    else if(inHistName[inHistNum].find("DPhi") != std::string::npos){
      if(th1Name.find("_Mean_") != std::string::npos) legStr = "Hist Mean";
      else if(th1Name.find("_FitMean_") != std::string::npos) legStr = "Fit Mean";
      else if(th1Name.find("_Res_") != std::string::npos) legStr = "Hist Res.";
      else if(th1Name.find("_FitRes_") != std::string::npos) legStr = "Fit Res.";
    }

    th1Canv_p[dirPos]->SetHist(th1_p[iter], centPos, 0, th1CanvCount[dirPos][centPos][0], legStr);
    th1CanvCount[dirPos][centPos][0]++;
  }

  if(debugMode) std::cout << __LINE__ << std::endl;

  for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
    th1Canv_p[dirIter]->canv_p->cd();
    th1Canv_p[dirIter]->MakeHistMaxMinNice();
    th1Canv_p[dirIter]->SetHistMaxMin();

    th1Canv_p[dirIter]->SetPanelWhiteSpace();

    if(debugMode) std::cout << __LINE__ << ", " << dirIter << std::endl;
	  
    for(Int_t iter = 0; iter < nXPanel; iter++){
      for(Int_t iter2 = 0; iter2 < nYPanel; iter2++){
	if(debugMode) std::cout << __LINE__ << ", " << dirIter << ", " << iter << ", " << iter2 << std::endl;

	th1Canv_p[dirIter]->canv_p->cd();
	th1Canv_p[dirIter]->pads_p[iter][iter2]->cd();
	
	for(Int_t histIter = 0; histIter < th1CanvCount[dirIter][iter][iter2]; histIter++){
	  if(debugMode) std::cout << __LINE__ << ", " << dirIter << ", " << iter << ", " << iter2 << ", " << histIter << std::endl;

	  th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->SetMarkerSize(2);
	  th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->SetMarkerStyle(20);

	  if(histIter == 0) th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->DrawCopy("E1 P");
	  else th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->DrawCopy("SAME E1 P");
	}

	if(debugMode) std::cout << __LINE__ << std::endl;

	if(iter == 0) th1Canv_p[dirIter]->DrawLabel1(iter, 0, dirNames_p->at(dirIter));
	if(isPbPb) th1Canv_p[dirIter]->DrawLabel2(iter, 0, centStrings2[iter]);
	else th1Canv_p[dirIter]->DrawLabel2(iter, 0, "PP");
      }
    }
    
    if(debugMode) std::cout << __LINE__ << std::endl;

    //    th1Canv_p[dirIter]->DrawLegend();

    if(debugMode) std::cout << __LINE__ << std::endl;

    if(ptEtaStr[ptEtaNum].find("VPt") != std::string::npos){
      for(Int_t iter = 0; iter < nXPanel; iter++){
	th1Canv_p[dirIter]->SetColumnLogX(iter); 
      }
    }
    
  }

  if(debugMode) std::cout << __LINE__ << std::endl;


  std::string outName = inFileName;
  const std::string inString = "_HIST.root";
  const std::string outString = "_PLOT.root";
  std::size_t strIndex = 0;

  strIndex = outName.find(inString);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, inString.length(), outString);
  }
  //EDITING  

  TF1* csn_p = 0;
  TF1* csn2_p = 0;
  TF1* csn3_p = 0;

  TF1* csn_PP_p = 0;
  TF1* csn2_PP_p = 0;
  TF1* csn3_PP_p = 0;

  if(isPbPb){
    csn_p = new TF1("csn_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))", 30, 150);
    csn_p->FixParameter(0, 0.061);
    csn_p->FixParameter(1, 1.24);
    csn_p->FixParameter(2, 8.08);
    
    csn2_p= new TF1("csn2_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))", 30, 150);
    csn2_p->FixParameter(0, 0.061 - .001);
    csn2_p->FixParameter(1, 1.24 - .04);
    csn2_p->FixParameter(2, 8.08 - .15);
    csn2_p->SetLineStyle(2);
    
    csn3_p = new TF1("csn3_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))", 30, 150);
    csn3_p->FixParameter(0, 0.061 + .001);
    csn3_p->FixParameter(1, 1.24 + .04);
    csn3_p->FixParameter(2, 8.08 + .15);
    csn3_p->SetLineStyle(2);
  }
  else{
    csn_p = new TF1("csn_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x)", 30, 150);
    csn_p->FixParameter(0, 0.061);
    csn_p->FixParameter(1, 1.24);
    
    csn2_p= new TF1("csn2_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x)", 30, 150);
    csn2_p->FixParameter(0, 0.061 - .001);
    csn2_p->FixParameter(1, 1.24 - .04);
    csn2_p->SetLineStyle(2);
    
    csn3_p = new TF1("csn3_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x)", 30, 150);
    csn3_p->FixParameter(0, 0.061 + .001);
    csn3_p->FixParameter(1, 1.24 + .04);
    csn3_p->SetLineStyle(2);

    csn_PP_p = new TF1("csn_PP_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x)", 30, 150);
    csn_PP_p->FixParameter(0, 0.061);
    csn_PP_p->FixParameter(1, .95);
    csn_PP_p->SetLineColor(kBlue);
    csn_PP_p->SetMarkerColor(kBlue);

    csn2_PP_p= new TF1("csn2_PP_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x)", 30, 150);
    csn2_PP_p->FixParameter(0, 0.061 - .001);
    csn2_PP_p->FixParameter(1, .95 - .04);
    csn2_PP_p->SetLineStyle(2);
    csn2_PP_p->SetLineColor(kBlue);
    csn2_PP_p->SetMarkerColor(kBlue);
    
    csn3_PP_p = new TF1("csn3_PP_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x)", 30, 150);
    csn3_PP_p->FixParameter(0, 0.061 + .001);
    csn3_PP_p->FixParameter(1, .95 + .04);
    csn3_PP_p->SetLineStyle(2);
    csn3_PP_p->SetLineColor(kBlue);
    csn3_PP_p->SetMarkerColor(kBlue);
  }

  TFile* outFile_p = new TFile(outName.c_str(), "UPDATE");

  for(Int_t iter = 0; iter < nDir; iter++){
    if(inHistName[inHistNum].find("RecoOverGen") != std::string::npos){
      if(meanResStr[meanResNum].find("Mean") != std::string::npos){
	th1Canv_p[iter]->DrawGlobalHorizontalLine(1);
	th1Canv_p[iter]->DrawGlobalHorizontalLine(.95);
	th1Canv_p[iter]->DrawGlobalHorizontalLine(1.05);
      }
      else if(meanResStr[meanResNum].find("Res") != std::string::npos){
	th1Canv_p[iter]->pads_p[0][0]->cd();
	csn_p->Draw("SAME C");
	csn2_p->Draw("SAME C");
	csn3_p->Draw("SAME C");
	
	if(!isPbPb){
	  csn_PP_p->Draw("SAME C");
	  csn2_PP_p->Draw("SAME C");
	  csn3_PP_p->Draw("SAME C");
	}
      }

      th1Canv_p[iter]->DrawGlobalVerticalLine(30);
    }

    if(inHistName[inHistNum].find("Eff") != std::string::npos && meanResStr[meanResNum].find("Mean") != std::string::npos){
      th1Canv_p[iter]->DrawGlobalHorizontalLine(1);
      th1Canv_p[iter]->DrawGlobalVerticalLine(30);
    }

    if(inHistName[inHistNum].find("DPhi") != std::string::npos){
      th1Canv_p[iter]->DrawGlobalVerticalLine(30);
    }

    std::string resCorrStr = "";
    if(addResCorrStr) resCorrStr = "_RESCORR";

    th1Canv_p[iter]->canv_p->Write(Form("%s_%s%s%sc", dirNames_p->at(iter).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str()), TObject::kOverwrite);
    claverCanvasSaving(th1Canv_p[iter]->canv_p, Form("pdfDir/%s_%s%s%sc", dirNames_p->at(iter).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str()), "pdf");
  }

  outFile_p->Close();
  delete outFile_p;

  delete meanLeg_p;

  delete label_p;

  

  for(Int_t iter = 0; iter < nDir; iter++){
    delete th1Canv_p[iter]->canv_p;
    delete th1Canv_p[iter];
  }

  th1Names_p->clear();
  delete th1Names_p;

  dirNames_p->clear();
  delete dirNames_p;

  inFile_p->Close();
  delete inFile_p;

  return 0;
}


int makeJECPlotMeanPts(const std::string inFileName, jecConfigParser config, const Int_t inHistNum, const Int_t ptEtaNum, const Int_t qgNum)
{
  const Bool_t isPbPb = config.GetIsPbPb();

  const Int_t nCentBins = config.GetNCentBins();
  std::vector<unsigned int> centBins = config.GetCentBins();

  if(isPbPb){
    if((unsigned int)(nCentBins+1) != centBins.size()){
      std::cout << "Config nCentBins \'" << nCentBins << "\' != centBins.size \'" << centBins.size() << "\'. return 1" << std::endl;
      return 1;
    }
  }

  std::string centStrings[nCentBins];
  std::string centStrings2[nCentBins];

  for(Int_t iter = 0; iter < nCentBins; iter++){
    if(isPbPb){
      centStrings[iter] = "Cent" + std::to_string(centBins.at(iter+1)) + "to" + std::to_string(centBins.at(iter));
      centStrings2[iter] = std::to_string(centBins.at(iter+1)) + "-" + std::to_string(centBins.at(iter)) + "%";
    }
    else{
      centStrings[iter] = "PP";
      centStrings2[iter] = "PP";
    }
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nContents = inFile_p->GetListOfKeys()->GetEntries();

  Int_t nTH1Temp = 0;
  std::vector<std::string>* th1Names_p = new std::vector<std::string>;

  Int_t nDirTemp = 0;
  std::vector<std::string>* dirNames_p = new std::vector<std::string>;

  Int_t nPtBinsTemp = 0;
  std::vector<std::string>* ptBins_p = new std::vector<std::string>;
  std::vector<Int_t>* ptInts_p = new std::vector<Int_t>;
  std::vector<Int_t>* ptDec_p = new std::vector<Int_t>;

  for(Int_t iter = 0; iter < nContents; iter++){
    TString className = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetClassName();
    TString name = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetName();

    if(className.Index("TDirectoryFile") >= 0 && name.Index("ak") >= 0){
      nDirTemp++;
      dirNames_p->push_back(name.Data());

      TDirectoryFile* tempDir_p = (TDirectoryFile*)inFile_p->Get(name);

      const Int_t nDirContents = tempDir_p->GetListOfKeys()->GetEntries();
      
      for(Int_t dirIter = 0; dirIter < nDirContents; dirIter++){
	TString name2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetName();
	TString className2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetClassName();

	if(className2.Index("TH1") >= 0){
	  if(name2.Index(ptEtaStr[ptEtaNum].c_str()) < 0) continue;
	  if(name2.Index(Form("_%s_", qgStr[qgNum].c_str())) < 0) continue;


	  if(name2.Index(Form("_%s", ptEtaStr2[ptEtaNum].c_str())) >= 0){
	    if(name2.Index(inHistName[inHistNum].c_str()) >= 0){
	      
	      TString name3 = name2(name2.Index(Form("_%s", ptEtaStr2[ptEtaNum].c_str()))+1, name2.Length()); 
	      TString name4 = name3(0, name3.Index("_"));
	      std::string ptRange = name4.Data();

	      Bool_t addPtRange = true;

	      for(Int_t ptIter = 0; ptIter < (Int_t)ptBins_p->size(); ptIter++){
		if(!strcmp(ptRange.c_str(), ptBins_p->at(ptIter).c_str())){
		  addPtRange = false;
		  break;
		}
	      }

	      if(addPtRange){
		ptBins_p->push_back(ptRange);
		ptInts_p->push_back(0);
		ptDec_p->push_back(0);
		nPtBinsTemp++;
	      }

	      nTH1Temp++;
	      th1Names_p->push_back(Form("%s/%s", name.Data(), name2.Data()));
	    }
	  }
	}
      }
    }
  }

  ptInts_p->push_back(0);
  ptDec_p->push_back(0);

  Int_t sortIter = 0; 

  while(sortIter < (Int_t)(ptBins_p->size() - 1)){
    TString ptBins1 = ptBins_p->at(sortIter).c_str();
    TString ptBins2 = ptBins_p->at(sortIter+1).c_str();

    TString pt1Part1 = ptBins1(ptBins1.Index(ptEtaStr2[ptEtaNum].c_str()) + ptEtaStr2[ptEtaNum].size(), TString(ptBins1(ptBins1.Index(ptEtaStr2[ptEtaNum].c_str())+ptEtaStr2[ptEtaNum].size(), ptBins1.Length() - ptEtaStr2[ptEtaNum].size())).Index("p"));
    TString pt1Part2 = ptBins1(ptBins1.Index("p")+1, 1);
    if(!strcmp(pt1Part2.Data(), "-")) pt1Part2 = ptBins1(ptBins1.Index("p")+1, 2);

    TString pt2Part1 = ptBins2(ptBins2.Index(ptEtaStr2[ptEtaNum].c_str()) + ptEtaStr2[ptEtaNum].size(), TString(ptBins2(ptBins2.Index(ptEtaStr2[ptEtaNum].c_str()) + ptEtaStr2[ptEtaNum].size(), ptBins2.Length() - ptEtaStr2[ptEtaNum].size())).Index("p"));
    TString pt2Part2 = ptBins2(ptBins2.Index("p")+1, 1);
    if(!strcmp(pt2Part2.Data(), "-")) pt2Part2 = ptBins2(ptBins2.Index("p")+1, 2);

    Float_t pt1 = 0;
    Float_t pt2 = 0;
    Bool_t pt1IsNeg = false;
    Bool_t pt2IsNeg = false;

    if(pt1Part1.Index("-") >= 0){
      pt1Part1 = pt1Part1(1, pt1Part1.Length() - 1);
      pt1IsNeg = true;
    }
    if(pt1Part2.Index("-") >= 0){
      pt1Part2 = pt1Part2(1, pt1Part2.Length() - 1);
      pt1IsNeg = true;
    }
    if(pt2Part1.Index("-") >= 0){
      pt2Part1 = pt2Part1(1, pt2Part1.Length() - 1);
      pt2IsNeg = true;
    }
    if(pt2Part2.Index("-") >= 0){
      pt2Part2 = pt2Part2(1, pt2Part2.Length() - 1);
      pt2IsNeg = true;
    }
  
    for(Int_t iter = 0; iter < 10; iter++){
      for(Int_t iter2 = 0; iter2 < Int_t(pt1Part1.Length()); iter2++){
	if(((TString)(pt1Part1(iter2, 1))).Index(Form("%d", iter)) >= 0) pt1 += TMath::Power(10, pt1Part1.Length() - iter2 - 1)*iter;
      }

      for(Int_t iter2 = 0; iter2 < Int_t(pt2Part1.Length()); iter2++){
	if(((TString)(pt2Part1(iter2, 1))).Index(Form("%d", iter)) >= 0) pt2 += TMath::Power(10, pt2Part1.Length() - iter2 - 1)*iter;
      }


      if(((TString)(pt1Part2(0, 1))).Index(Form("%d", iter)) >= 0) pt1 += iter/10.;
      if(((TString)(pt2Part2(0, 1))).Index(Form("%d", iter)) >= 0) pt2 += iter/10.;

    }



    Float_t tempPt1 = pt1;
    Float_t tempPt2 = pt2;

    if(pt1IsNeg) pt1 = -pt1;
    if(pt2IsNeg) pt2 = -pt2;
   
    if(pt1 > pt2){
      std::string tempStr = ptBins_p->at(sortIter); 
      ptBins_p->at(sortIter) = ptBins_p->at(sortIter+1);
      ptBins_p->at(sortIter+1) = tempStr;
      
      ptInts_p->at(sortIter) = std::floor(tempPt1);
      ptDec_p->at(sortIter) = tempPt1*10 - std::floor(tempPt1)*10;

      if(pt1IsNeg) ptInts_p->at(sortIter) = -ptInts_p->at(sortIter);

      ptInts_p->at(sortIter+1) = std::floor(tempPt2);
      ptDec_p->at(sortIter+1) = tempPt2*10 - std::floor(tempPt2)*10;

      if(pt2IsNeg) ptInts_p->at(sortIter+1) = -ptInts_p->at(sortIter+1);

      sortIter = 0;
    }
    else{
      ptInts_p->at(sortIter) = std::floor(tempPt1);
      ptDec_p->at(sortIter) = tempPt1*10 - std::floor(tempPt1)*10;

      if(pt1IsNeg) ptInts_p->at(sortIter) = -ptInts_p->at(sortIter);

      ptInts_p->at(sortIter+1) = std::floor(tempPt2);
      ptDec_p->at(sortIter+1) = tempPt2*10 - std::floor(tempPt2)*10;

      if(pt2IsNeg) ptInts_p->at(sortIter+1) = -ptInts_p->at(sortIter+1);

      sortIter++;
    }
  }

  Int_t nCentBinsTemp = 1;
  if(isPbPb) nCentBinsTemp = nCentBins;
  const Int_t nCentBins2 = nCentBinsTemp;

  const Int_t nTH1 = nTH1Temp;
  const Int_t nDir = nDirTemp;
  const Int_t nPtBins = nPtBinsTemp;
  TH1F* th1_p[nTH1];
  TCanvas* th1Canv_p[nDir][nCentBins2];

  Int_t xBins = 0;
  Int_t yBins = 0;
  Float_t breakBool = false;

  for(Int_t iter = 0; iter < nPtBins; iter++){
    if(breakBool) break;
    for(Int_t iter2 = 0; iter2 < 4; iter2++){
   
      xBins = iter + 4 - iter2;
      yBins = iter;

      //      if(xBins == nPtBins-1) continue;

      if(xBins >= 6){
	if(xBins*yBins > nPtBins && xBins*yBins - nPtBins < 4){
	  breakBool = true;
	  break;
	}
      }
      else if(xBins*yBins > nPtBins && xBins*yBins - nPtBins < ((Float_t)xBins)/2.){
	breakBool = true;
	break;
      }
    }
  }

  if(xBins == nPtBins && nPtBins >= 8){
    xBins = nPtBins/2;
    yBins = 2;

    if(xBins*yBins < nPtBins) yBins = 3;
  }

  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      std::string centStr = centStrings[centIter];
      if(!isPbPb) centStr = "PP";

      th1Canv_p[iter][centIter] = new TCanvas(Form("%s_%s_MeanPts%s%s%s_c", dirNames_p->at(iter).c_str(), qgStr[qgNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str(), centStr.c_str()), Form("%s_%s_MeanPts%s%s%s_c", dirNames_p->at(iter).c_str(), qgStr[qgNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str(), centStr.c_str()), xBins*300, yBins*325);
      //      th1Canv_p[iter][centIter]->Divide(xBins, yBins, 0.0, 0.0);
      th1Canv_p[iter][centIter]->Divide(xBins, yBins, .000005, .000005);
    }
  }

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  Float_t dirMax[nDir][nCentBins2];
  Float_t dirMin[nDir][nCentBins2];

  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      dirMax[iter][centIter] = -1;
      dirMin[iter][centIter] = 100;
    }
  }

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    th1_p[iter] = (TH1F*)inFile_p->Get(th1Name.c_str());

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));

      if(th1Name.substr(0, th1Name.find("/")).size() != dirNames_p->at(dirIter).size()) continue;

      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    Int_t centPos = -1;
    if(isPbPb){
      for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
	std::size_t pos = th1Name.find(centStrings[centIter]);
	if(pos != std::string::npos){
	  centPos = centIter;
	  break;
	}
      }
    }
    else centPos = 0;



    for(Int_t binIter = 0; binIter < th1_p[iter]->GetNbinsX(); binIter++){
      if(th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1) > dirMax[dirPos][centPos]) dirMax[dirPos][centPos] = th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1);

      if(th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1) < dirMin[dirPos][centPos]) dirMin[dirPos][centPos] = th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1);
    }
  }


  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      dirMax[iter][centIter] = setMaxMinNice(dirMax[iter][centIter], true);
      dirMin[iter][centIter] = setMaxMinNice(dirMin[iter][centIter], false);
    }
  }

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));

      if(th1Name.substr(0, th1Name.find("/")).size() != dirNames_p->at(dirIter).size()) continue;

      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    Int_t centPos = -1;
    if(isPbPb){
      for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
	std::size_t pos = th1Name.find(centStrings[centIter]);
	if(pos != std::string::npos){
	  centPos = centIter;
	  break;
	}
      }
    }
    else centPos = 0;

    Int_t ptPos = -1;
    for(Int_t ptIter = 0; ptIter < (Int_t)ptBins_p->size(); ptIter++){
      std::size_t pos = th1Name.find(ptBins_p->at(ptIter).c_str());
      if(pos != std::string::npos){
        ptPos = ptIter;
        break;
      }
    }

    th1Canv_p[dirPos][centPos]->cd(ptPos+1);

    gPad->SetLeftMargin(gPad->GetLeftMargin()/2);
    gPad->SetBottomMargin(gPad->GetBottomMargin()/2);

    gStyle->SetOptStat(0);

    //    th1_p[iter]->SetMaximum(dirMax[dirPos][centPos]);
    //    th1_p[iter]->SetMinimum(dirMin[dirPos][centPos]);

    th1_p[iter]->GetYaxis()->SetNdivisions(404);

    th1_p[iter]->GetXaxis()->CenterTitle();
    th1_p[iter]->GetXaxis()->SetTitleOffset(th1_p[iter]->GetXaxis()->GetTitleOffset() + 2.0);
    th1_p[iter]->GetXaxis()->SetTitleFont(43);
    th1_p[iter]->GetXaxis()->SetTitleSize(26);
    th1_p[iter]->GetYaxis()->CenterTitle();
    th1_p[iter]->GetYaxis()->SetTitleOffset(th1_p[iter]->GetYaxis()->GetTitleOffset() + 2.5);
    
    th1_p[iter]->GetXaxis()->SetTitle(xAxisLabel[inHistNum].c_str());
    th1_p[iter]->GetYaxis()->SetTitle("Events");

    th1_p[iter]->DrawCopy("E1 P");
    
    TLine* meanLine_p = new TLine(th1_p[iter]->GetMean(), th1_p[iter]->GetMinimum(), th1_p[iter]->GetMean(), th1_p[iter]->GetMaximum()*.8);
    meanLine_p->SetLineStyle(2);
    meanLine_p->DrawClone();
    delete meanLine_p;
    
    if(th1_p[iter]->GetListOfFunctions()->GetSize() != 0){
      TLine* meanLine2_p = new TLine(th1_p[iter]->GetFunction("f1_p")->GetParameter(1), th1_p[iter]->GetMinimum(), th1_p[iter]->GetFunction("f1_p")->GetParameter(1), th1_p[iter]->GetMaximum()*.8);
      meanLine2_p->SetLineStyle(2);
      meanLine2_p->SetLineColor(kRed);
      meanLine2_p->DrawClone();
      delete meanLine2_p;
    }

    //    if(centPos == 0) label_p->DrawLatex(.30, .9, Form("#bf{%s}", dirNames_p->at(dirPos).c_str()));
    
    if(ptPos%xBins == 0) label_p->DrawLatex(.30, .88, Form("%d.%d<%s<%d.%d", ptInts_p->at(ptPos), ptDec_p->at(ptPos), ptEtaStr3[ptEtaNum].c_str(), ptInts_p->at(ptPos+1), ptDec_p->at(ptPos+1)));
    else label_p->DrawLatex(.20, .88, Form("%d.%d<%s<%d.%d", ptInts_p->at(ptPos), ptDec_p->at(ptPos), ptEtaStr3[ptEtaNum].c_str(), ptInts_p->at(ptPos+1), ptDec_p->at(ptPos+1)));


    if(ptPos == 0){
      if(isPbPb) label_p->DrawLatex(.30, .78, Form("#bf{#color[2]{%s (%s)}}", dirNames_p->at(dirPos).c_str(), centStrings2[centPos].c_str()));
      else label_p->DrawLatex(.30, .78, Form("#bf{#color[2]{%s (%s)}}", dirNames_p->at(dirPos).c_str(), "PP"));
    }
    if(ptPos == 1) label_p->DrawLatex(.15, .78, Form("#bf{#color[2]{|#eta_{jet}|<1.6}}"));
    if(ptPos == 2 && !strcmp(ptEtaStr2[ptEtaNum].c_str(), "Eta" )) label_p->DrawLatex(.15, .78, Form("#bf{#color[2]{Gen. p_{T}>30}}"));
    if(ptPos == 3) label_p->DrawLatex(.15, .78, Form("#bf{#color[2]{%s}}", qgStr2[qgNum].c_str()));

    //    if(centPos%xBins == 0) label_p->DrawLatex(.68, .9, centStrings2[centPos].c_str());
    //    else label_p->DrawLatex(.60, .9, centStrings2[centPos].c_str());
  }

  std::string outName = inFileName;
  const std::string inString = "_HIST.root";
  const std::string outString = "_PLOT.root";
  std::size_t strIndex = 0;

  strIndex = outName.find(inString);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, inString.length(), outString);
  }


  TFile* outFile_p = new TFile(outName.c_str(), "UPDATE");

  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      th1Canv_p[iter][centIter]->Write("", TObject::kOverwrite);
      claverCanvasSaving(th1Canv_p[iter][centIter], Form("pdfDir/%s", th1Canv_p[iter][centIter]->GetName()), "pdf");
    }
  }

  outFile_p->Close();
  delete outFile_p;

  delete label_p;

  delete ptBins_p;
  delete ptInts_p;
  delete ptDec_p;


  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      delete th1Canv_p[iter][centIter];
    }
  }

  th1Names_p->clear();
  delete th1Names_p;

  dirNames_p->clear();
  delete dirNames_p;

  inFile_p->Close();
  delete inFile_p;

  return 0;
}


int makeJECPlot(const std::string inFileName)
{
  jecConfigParser config;
  if(!config.SetConfigParser(inFileName)) return 1;

  checkMakeDir("pdfDir");
  
  int retVal = 0;

  for(Int_t iter = 0; iter < nHistName; iter++){
    for(Int_t ptEtaIter = 0; ptEtaIter < nPtEta; ptEtaIter++){
      for(Int_t iter2 = 0; iter2 < nMeanRes; iter2++){
	
	if(iter == 0 || (iter == 2 && ptEtaIter < nPtEta-1) || (ptEtaIter < nPtEta-1 && iter2 == 0)) retVal += makeJECPlotMeanRes(inFileName, config, iter, iter2, ptEtaIter);
	
      }
      
      if(iter == 0){
	for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	  if(qgIter == 0) retVal += makeJECPlotMeanPts(inFileName, config, iter, ptEtaIter, qgIter);
	}
      }
    }
  }
  
  
  return retVal;
}

int main(int argc, char *argv[])
{
  if(argc != 2){
    std::cout << "Usage: makeJECPlot_Prototype.exe <inHistFile>" << std::endl;
    std::cout << "Number of args given: " << argc << std::endl;
    for(int iter = 0; iter < argc; iter++){
      std::cout << "  argv[" << iter << "]: " << argv[iter] << std::endl;
    }
    return -1;
  }

  int retVal = makeJECPlot(argv[1]);
  return retVal;
}
