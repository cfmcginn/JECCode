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

Bool_t plotTrue = false;

const Int_t nJetAlgo = 1;
const std::string jetAlgo[nJetAlgo] = {"akPu3PF"/*, "akCs3PF"*/};
//const std::string jetAlgo[nJetAlgo] = {"akVs4Calo", "akPu4Calo", "akVs4PF", "akPu4PF", "akVs3PF", "akPu3PF"};
//const std::string jetAlgo[nJetAlgo] = {"ak4Calo", "ak3Calo", "ak4PF", "ak3PF"};

//const Int_t nCentBins = 4;
//const std::string centStrings[nCentBins] = {"cent50to100", "cent30to50", "cent10to30", "cent0to10"};
//const std::string centStrings2[nCentBins] = {"50-100%", "30-50%", "10-30%", "0-10%"};

const Int_t nCentBins = 2;
const std::string centStrings[nCentBins] = {"cent30to100", "cent0to30"};
const std::string centStrings2[nCentBins] = {"30-100%", "0-30%"};

//const Int_t nCentBins = 8;
//const std::string centStrings[nCentBins] = {"cent70to100", "cent50to70", "cent40to50", "cent30to40", "cent20to30", "cent10to20", "cent5to10", "cent0to5"};
//const std::string centStrings2[nCentBins] = {"70-100%", "50-70%", "40-50%", "30-40%", "20-30%", "10-20%", "5-10%", "0-5%"};

const Int_t nHistName = 5;
const std::string inHistName[nHistName+4+10+1] = {"RecoOverGen", "RecoOverRaw", "RawOverGen", "Eff", "Fake", "RecoGenDR", "RecoVGen", "RecoGenDEta", "RecoGenDPhi", "RecoRho0OverGen", "RecoRho1OverGen", "RecoRho2OverGen", "RecoRho3OverGen", "RecoRho4OverGen", "RecoRho5OverGen", "RecoRho6OverGen", "RecoRho7OverGen", "RecoRho8OverGen", "RecoRho9OverGen", "RecoRho"};
const std::string xAxisLabel[nHistName] = {"Reco./Gen.", "Reco./Raw", "Raw/Gen.", "Eff.", "Fake"};

const Int_t nMeanRes = 2;
const std::string meanResStr[nMeanRes] = {"Mean", "Res"};
const std::string meanResStr2[nMeanRes] = {"#mu", "#sigma"};


const Int_t nQG = 3;
const std::string qgStr[nQG] = {"Inc", "Q", "G"};
const std::string qgStr2[nQG] = {"Inc.", "Quarks", "Gluons"};
const Int_t qgCol[nQG] = {kGray+1, kBlue, kRed};

const Int_t nPtCut = 5;
const std::string ptCutStr[nPtCut] = {"(p_{T}^{reco} > 5)", "(p_{T}^{reco} > 10)", "(p_{T}^{reco} > 15)", "(p_{T}^{reco} > 20)", "(p_{T}^{reco} > 25)"};
const Int_t ptCutCol[nPtCut] = {kGray+1, kBlue, kRed, kYellow+1, kMagenta};

const Int_t nPtEta = 3;
const std::string ptEtaStr[nPtEta] = {"VPt_", "VEta_", "VPtEta_"};
const std::string ptEtaStr2[nPtEta] = {"Pt", "Eta", "PtEta"};
const std::string ptEtaStr3[nPtEta] = {"p_{T}", "#eta", ""};


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


void makeJECPlotMeanRes(const std::string inFileName, const Int_t inHistNum, const Int_t meanResNum, const Int_t ptEtaNum, const Bool_t isPbPb)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nContents = inFile_p->GetListOfKeys()->GetEntries();

  Int_t nTH1Temp = 0;
  std::vector<std::string>* th1Names_p = new std::vector<std::string>;

  Int_t nDirTemp = 0;
  std::vector<std::string>* dirNames_p = new std::vector<std::string>;

  for(Int_t iter = 0; iter < nContents; iter++){
    TString className = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetClassName();
    TString name = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetName();

    if(className.Index("TDirectoryFile") >= 0){
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
	    if(!plotTrue && name2.Index(Form("Fit")) < 0 && strcmp("Eff", inHistName[inHistNum].c_str()) != 0 && strcmp("Fake", inHistName[inHistNum].c_str()) != 0 && strcmp("RecoGenDR", inHistName[inHistNum].c_str()) != 0 && strcmp("RecoGenDPhi", inHistName[inHistNum].c_str()) != 0 && strcmp("RecoGenDEta", inHistName[inHistNum].c_str()) != 0) continue;

	    if(name2.Index("PERP") >= 0) continue;

	    //	    if(name2.Index("Reco5") >= 0) continue;
	    //	    if(name2.Index("Reco10") >= 0) continue;
	    //	    if(name2.Index("Reco15") >= 0) continue;
	    //	    if(name2.Index("Reco20") >= 0) continue;
	    if(name2.Index("Reco30") >= 0) continue;

	    if(name2.Index("Res") >= 0 && name2.Index("_Q_") >= 0) continue;
	    if(name2.Index("Res") >= 0 && name2.Index("_G_") >= 0) continue;

	    if(name2.Index("_Q_") >= 0) continue;
	    if(name2.Index("_G_") >= 0) continue;

	    if(!strcmp("Fake", inHistName[inHistNum].c_str()) && (name2.Index("_Q_") >= 0 || name2.Index("_G_") >= 0)) continue;
	    if(!strcmp("Eff", inHistName[inHistNum].c_str()) && (name2.Index("_Q_") >= 0 || name2.Index("_G_") >= 0)) continue;

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
  
  TCanvas* th1Canv_p[nDir];

  Int_t nXPanel = nCentBins;
  Int_t nYPanel = 1;
  //  Int_t nYPanel = 2;
  if(!isPbPb){
    nXPanel = 1;
    nYPanel = 1;
  }

  for(Int_t iter = 0; iter < nDir; iter++){
    th1Canv_p[iter] = new TCanvas(Form("%s_%s%s%sc", dirNames_p->at(iter).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str()), Form("%s_%s%s%sc", dirNames_p->at(iter).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str()), nXPanel*300, nYPanel*325);
    th1Canv_p[iter]->Divide(nXPanel, nYPanel, 0.0, 0.0);
  }

  Bool_t isDrawn[nDir][nXPanel*nYPanel];
  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t iter2 = 0; iter2 < nXPanel*nYPanel; iter2++){
      isDrawn[iter][iter2] = false;
    }
  }
  
  Bool_t legAdded[2][nQG] = {{false}, {false}};

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  Float_t dirMax[nDir];
  Float_t dirMin[nDir];

  for(Int_t iter = 0; iter < nDir; iter++){
    dirMax[iter] = -1;
    dirMin[iter] = 100;
  }

  Float_t th1XMin = 30;

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


  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    th1_p[iter] = (TH1F*)inFile_p->Get(th1Name.c_str());

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
      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    Float_t dirMinBound = 0.5;
    if(!strcmp("Fake", inHistName[inHistNum].c_str())) dirMinBound = 0.0;
    else if(!strcmp("Res", meanResStr[meanResNum].c_str())) dirMinBound = 0.0;
    else if(!strcmp("RecoGenDEta", inHistName[inHistNum].c_str())) dirMinBound = -10000;
    else if(!strcmp("RecoGenDPhi", inHistName[inHistNum].c_str())) dirMinBound = -10000;

    Float_t dirMaxBound = 1.5;

    for(Int_t binIter = xMinBin-1; binIter < th1_p[iter]->GetNbinsX(); binIter++){
      if(th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1) > dirMax[dirPos] && th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1) < dirMaxBound) dirMax[dirPos] = th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1);

      if(th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1) < dirMin[dirPos] && th1_p[iter]->GetBinContent(binIter+1) != 0 && th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1) > dirMinBound) dirMin[dirPos] = th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1);
    }
  }

  Float_t dirMaxMin = 1.15;
  if(!strcmp("Res", meanResStr[meanResNum].c_str())) dirMaxMin = 0;

  if(strcmp(inHistName[inHistNum].c_str(), "RecoGenDR") != 0 && strcmp(inHistName[inHistNum].c_str(), "RecoGenDEta") != 0 && strcmp(inHistName[inHistNum].c_str(), "RecoGenDPhi") != 0){
    for(Int_t iter = 0; iter < nDir; iter++){
      if(dirMax[iter] < dirMaxMin) dirMax[iter] = dirMaxMin;
      
      if(!strcmp(inHistName[inHistNum].c_str(), "Eff")) dirMax[iter] = 1. + (1-dirMin[iter])/5;
      else{
	if(iter == 2 && !strcmp("Res", meanResStr[meanResNum].c_str())) std::cout << "max: " << dirMax[iter] << std::endl;
	
	dirMax[iter] = setMaxMinNice(dirMax[iter], true);
	dirMin[iter] = setMaxMinNice(dirMin[iter], false);
	
	if(iter == 2 && !strcmp("Res", meanResStr[meanResNum].c_str())) std::cout << "max: " << dirMax[iter] << std::endl;
      }
    }
  }

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
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

    if(centPos == -1){
      //      std::cout << "NO CENT POS FOUND" << std::endl;
      centPos = 0;
    }

    th1Canv_p[dirPos]->cd(centPos+1);

    th1_p[iter]->SetMaximum(dirMax[dirPos]);
    th1_p[iter]->SetMinimum(dirMin[dirPos]);

    th1_p[iter]->GetXaxis()->CenterTitle();
    th1_p[iter]->GetXaxis()->SetTitleOffset(1.);
    th1_p[iter]->GetXaxis()->SetTitleFont(43);
    th1_p[iter]->GetXaxis()->SetTitleSize(26);
    th1_p[iter]->GetXaxis()->SetNdivisions(404);
    th1_p[iter]->GetYaxis()->CenterTitle();
    th1_p[iter]->GetYaxis()->SetTitleOffset(1.1);
    th1_p[iter]->GetYaxis()->SetNdivisions(505);    
    
    if(!strcmp(ptEtaStr2[ptEtaNum].c_str(), "Pt")) gPad->SetLogx();

    std::string legString;
    
    Int_t colPos = 0;
    if(strcmp("Eff", inHistName[inHistNum].c_str()) != 0){
      if(th1Name.find("_Q_") != std::string::npos){
	colPos = 1;
	legString = qgStr2[1];
      }
      else if(th1Name.find("_G_") != std::string::npos){
	colPos = 2;
	legString = qgStr2[2];
      }
      else legString = qgStr2[0];
    }
    else{
      if(th1Name.find("Reco10") != std::string::npos){
	colPos = 1;
	legString = ptCutStr[1];
      }
      else if(th1Name.find("Reco15") != std::string::npos){
	colPos = 2;
	legString = ptCutStr[2];
      }
      else if(th1Name.find("Reco20") != std::string::npos){
	colPos = 3;
	legString = ptCutStr[3];
      }
      else if(th1Name.find("Reco25") != std::string::npos){
	colPos = 4;
	legString = ptCutStr[4];
      }
      else legString = ptCutStr[0];
    }
    
    if(th1Name.find("Fit") != std::string::npos){
      th1_p[iter]->SetMarkerColor(ptCutCol[colPos]);
      th1_p[iter]->SetLineColor(ptCutCol[colPos]);
      if(!legAdded[0][colPos]){
	if(!strcmp(inHistName[inHistNum].c_str(), "Eff")) meanLeg_p->AddEntry(th1_p[iter], Form("%s", legString.c_str()), "P L");
	else meanLeg_p->AddEntry(th1_p[iter], Form("Fit %s %s", meanResStr2[meanResNum].c_str(), legString.c_str()), "P L");
	legAdded[0][colPos] = true;
      }
    }
    else{
      th1_p[iter]->SetMarkerColor(ptCutCol[colPos]);
      th1_p[iter]->SetMarkerStyle(24);
      th1_p[iter]->SetLineColor(ptCutCol[colPos]);

      if(!legAdded[1][colPos]){
	if(!strcmp(inHistName[inHistNum].c_str(), "Eff")){
	  if(isPbPb) meanLeg2_p->AddEntry(th1_p[iter], Form("%s", legString.c_str()), "P L");
	  else meanLeg_p->AddEntry(th1_p[iter], Form("%s", legString.c_str()), "P L");
	}
	else{
	  if(isPbPb) meanLeg2_p->AddEntry(th1_p[iter], Form("True %s %s", meanResStr2[meanResNum].c_str(), legString.c_str()), "P L");
	  else meanLeg_p->AddEntry(th1_p[iter], Form("True %s %s", meanResStr2[meanResNum].c_str(), legString.c_str()), "P L");
	}

	legAdded[1][colPos] = true;
      }
    }

    if(!isDrawn[dirPos][centPos]){
      th1_p[iter]->DrawCopy("E1 P");
      isDrawn[dirPos][centPos] = true;
    }
    else th1_p[iter]->DrawCopy("E1 P SAME");
    
    if(dirMin[dirPos] < 1 && dirMax[dirPos] > 1 && strcmp("Res", meanResStr[meanResNum].c_str()) != 0){
      TLine* oneLine_p = new TLine(th1XMin, 1, th1_p[iter]->GetXaxis()->GetXmax(), 1);
      oneLine_p->SetLineStyle(2);
      oneLine_p->DrawClone();
      delete oneLine_p;

      TLine* thirtyLine_p = new TLine(30, dirMin[dirPos], 30, dirMax[dirPos]);
      thirtyLine_p->SetLineStyle(2);
      thirtyLine_p->DrawClone();
      delete thirtyLine_p;

      if(!strcmp(inHistName[inHistNum].c_str(), "RecoOverGen")){
	TLine* oneLineUp_p = new TLine(th1XMin, 1.05, th1_p[iter]->GetXaxis()->GetXmax(), 1.05);
	oneLineUp_p->SetLineStyle(2);
	oneLineUp_p->DrawClone();
	delete oneLineUp_p;

	TLine* oneLineDown_p = new TLine(th1XMin, .95, th1_p[iter]->GetXaxis()->GetXmax(), .95);
	oneLineDown_p->SetLineStyle(2);
	oneLineDown_p->DrawClone();
	delete oneLineDown_p;
      }
    }    

    if(centPos == 0) label_p->DrawLatex(.30, .9, Form("#bf{#color[2]{%s}}", dirNames_p->at(dirPos).c_str()));

    if(centPos == 1) label_p->DrawLatex(.10, .9, Form("#bf{#color[2]{|#eta_{jet}|<2.0}}"));
    if(centPos == 2 && !strcmp(ptEtaStr2[ptEtaNum].c_str(), "Eta" )) label_p->DrawLatex(.10, .9, Form("#bf{#color[2]{Gen. p_{T}>30}}"));

    if(isPbPb){
      if(centPos%4 == 0 && centPos < 4) label_p->DrawLatex(.28, .25, centStrings2[centPos].c_str());
      else if(centPos%4 == 0 && centPos > 3) label_p->DrawLatex(.28, .25, centStrings2[centPos].c_str());
      else if(centPos%4 != 0 && centPos < 4) label_p->DrawLatex(.08, .25, centStrings2[centPos].c_str());
      else label_p->DrawLatex(.08, .25, centStrings2[centPos].c_str());
    }
    else{
      if(centPos%4 == 0) label_p->DrawLatex(.78, .9, "PP");
    }
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
    if(!isPbPb) th1Canv_p[iter]->cd();
    else th1Canv_p[iter]->cd(2);
    //    else th1Canv_p[iter]->cd(3);
    meanLeg_p->Draw("SAME");

    if(plotTrue || !strcmp(inHistName[inHistNum].c_str(), "Eff")){
      if(!isPbPb) th1Canv_p[iter]->cd();
      else if(!strcmp("Eff", inHistName[inHistNum].c_str())) th1Canv_p[iter]->cd(1);
      else th1Canv_p[iter]->cd(4);

      if(isPbPb) meanLeg2_p->Draw("SAME");      
    }

    th1Canv_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(th1Canv_p[iter], Form("pdfDir/%s", th1Canv_p[iter]->GetName()), "pdf");
    claverCanvasSaving(th1Canv_p[iter], Form("pdfDir/%s", th1Canv_p[iter]->GetName()), "C");
  }

  outFile_p->Close();
  delete outFile_p;

  delete meanLeg_p;

  delete label_p;

  for(Int_t iter = 0; iter < nDir; iter++){
    delete th1Canv_p[iter];
  }

  th1Names_p->clear();
  delete th1Names_p;

  dirNames_p->clear();
  delete dirNames_p;

  inFile_p->Close();
  delete inFile_p;

  return;
}


void makeJECPlotMeanRes_Cent(const std::string inFileName, const Int_t inHistNum, const Int_t meanResNum, const Bool_t isPbPb)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nContents = inFile_p->GetListOfKeys()->GetEntries();

  Int_t nTH1Temp = 0;
  std::vector<std::string>* th1Names_p = new std::vector<std::string>;

  Int_t nDirTemp = 0;
  std::vector<std::string>* dirNames_p = new std::vector<std::string>;

  std::vector<std::string>* ptNames_p = new std::vector<std::string>;
  std::vector<std::string>* etaNames_p = new std::vector<std::string>;

  for(Int_t iter = 0; iter < nContents; iter++){
    TString className = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetClassName();
    TString name = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetName();

    if(className.Index("TDirectoryFile") >= 0){
      nDirTemp++;
      dirNames_p->push_back(name.Data());

      TDirectoryFile* tempDir_p = (TDirectoryFile*)inFile_p->Get(name);

      const Int_t nDirContents = tempDir_p->GetListOfKeys()->GetEntries();
      
      for(Int_t dirIter = 0; dirIter < nDirContents; dirIter++){
	TString name2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetName();
	TString className2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetClassName();

	if(className2.Index("TH1") >= 0){
	  if(name2.Index(Form("%s_", meanResStr[meanResNum].c_str())) >= 0){
	    if(name2.Index("VCent") < 0) continue;
	    if(!plotTrue && name2.Index(Form("Fit")) < 0) continue;

	    if(name2.Index(inHistName[inHistNum].c_str()) >= 0){
	      nTH1Temp++;
	      th1Names_p->push_back(Form("%s/%s", name.Data(), name2.Data()));

	      std::string tempPtStr = th1Names_p->at(nTH1Temp-1).substr(th1Names_p->at(nTH1Temp-1).find("_Pt")+1, th1Names_p->at(nTH1Temp-1).length());
	      std::string tempEtaStr = th1Names_p->at(nTH1Temp-1).substr(th1Names_p->at(nTH1Temp-1).find("_Eta")+1, th1Names_p->at(nTH1Temp-1).length());

	      tempPtStr = tempPtStr.substr(0, tempPtStr.find("_"));
	      tempEtaStr = tempEtaStr.substr(0, tempEtaStr.find("_"));

	      Bool_t isPtVect = false;
	      Bool_t isEtaVect = false;

	      for(Int_t ptIter = 0; ptIter < (Int_t)ptNames_p->size(); ptIter++){
		if(!strcmp(ptNames_p->at(ptIter).c_str(), tempPtStr.c_str())){
		  isPtVect = true;
		  break;
		}
	      }

	      for(Int_t etaIter = 0; etaIter < (Int_t)etaNames_p->size(); etaIter++){
		if(!strcmp(etaNames_p->at(etaIter).c_str(), tempEtaStr.c_str())){
		  isEtaVect = true;
		  break;
		}
	      }

	      if(!isPtVect) ptNames_p->push_back(tempPtStr);
	      if(!isEtaVect) etaNames_p->push_back(tempEtaStr);
	    }
	  }
	}
      }
    }
  }

  const Int_t nTH1 = nTH1Temp;
  const Int_t nDir = nDirTemp;
  TH1F* th1_p[nTH1];
  
  const Int_t nPtCanv = (Int_t)ptNames_p->size();

  TCanvas* th1Canv_p[nDir][nPtCanv];

  //editing here

  Int_t nXPanel = (Int_t)etaNames_p->size();
  Int_t nYPanel = 1;

  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t ptIter = 0; ptIter < nPtCanv; ptIter++){
      th1Canv_p[iter][ptIter] = new TCanvas(Form("%s_%s%sVCent_%s_c", dirNames_p->at(iter).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str(), ptNames_p->at(ptIter).c_str()), Form("%s_%s%sVCent_%s_c", dirNames_p->at(iter).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str(), ptNames_p->at(ptIter).c_str()), nXPanel*300, nYPanel*325);
      th1Canv_p[iter][ptIter]->Divide(nXPanel, nYPanel, 0.0, 0.0);
    }
  }

  Bool_t isDrawn[nDir][nPtCanv][nXPanel*nYPanel];
  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t ptIter = 0; ptIter < nPtCanv; ptIter++){
      for(Int_t iter2 = 0; iter2 < nXPanel*nYPanel; iter2++){
	isDrawn[iter][ptIter][iter2] = false;
      }
    }
  }

  Bool_t legAdded[2][nQG] = {{false}, {false}};

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  Float_t dirMax[nDir][nPtCanv];
  Float_t dirMin[nDir][nPtCanv];

  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t ptIter = 0; ptIter < nPtCanv; ptIter++){
      dirMax[iter][ptIter] = -1;
      dirMin[iter][ptIter] = 100;
    }
  }

  TLegend* meanLeg_p;
  if(isPbPb) meanLeg_p = new TLegend(.50, .68, .90, .95);
  else meanLeg_p = new TLegend(.50, .58, .90, .85);
  meanLeg_p->SetBorderSize(0);
  meanLeg_p->SetFillColor(0);
  meanLeg_p->SetFillStyle(0);
  meanLeg_p->SetTextFont(43);
  meanLeg_p->SetTextSize(18);


  TLegend* meanLeg2_p = new TLegend(.50, .68, .90, .95);
  meanLeg2_p->SetBorderSize(0);
  meanLeg2_p->SetFillColor(0);
  meanLeg2_p->SetFillStyle(0);
  meanLeg2_p->SetTextFont(43);
  meanLeg2_p->SetTextSize(18);


  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);
    th1_p[iter] = (TH1F*)inFile_p->Get(th1Name.c_str());

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    Int_t ptPos = -1;
    for(Int_t ptIter = 0; ptIter < nPtCanv; ptIter++){
      std::size_t pos = th1Name.find(ptNames_p->at(ptIter));
      if(pos != std::string::npos){
	ptPos = ptIter;
	break;
      }
    }

    for(Int_t binIter = 0; binIter < th1_p[iter]->GetNbinsX(); binIter++){
      if(th1_p[iter]->GetBinContent(binIter+1) > dirMax[dirPos][ptPos]) dirMax[dirPos][ptPos] = th1_p[iter]->GetBinContent(binIter+1);

      if(th1_p[iter]->GetBinContent(binIter+1) < dirMin[dirPos][ptPos] && th1_p[iter]->GetBinContent(binIter+1) != 0) dirMin[dirPos][ptPos] = th1_p[iter]->GetBinContent(binIter+1);
    }
  }

  

  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t ptIter = 0; ptIter < nPtCanv; ptIter++){
      dirMax[iter][ptIter] += (dirMax[iter][ptIter] - dirMin[iter][ptIter])/10.;
      dirMin[iter][ptIter] -= (dirMax[iter][ptIter] - dirMin[iter][ptIter])/10.;

      if(dirMax[iter][ptIter] < 1.2) dirMax[iter][ptIter] = 1.2;
      if(dirMin[iter][ptIter] > .8) dirMin[iter][ptIter] = .8;

      dirMax[iter][ptIter] = setMaxMinNice(dirMax[iter][ptIter], true);
      dirMin[iter][ptIter] = setMaxMinNice(dirMin[iter][ptIter], false);
    }
  }

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    Int_t ptPos = -1;
    for(Int_t ptIter = 0; ptIter < nPtCanv; ptIter++){
      std::size_t pos = th1Name.find(ptNames_p->at(ptIter));
      if(pos != std::string::npos){
        ptPos = ptIter;
        break;
      }
    }

    Int_t etaPos = -1;
    for(Int_t etaIter = 0; etaIter < nXPanel; etaIter++){
      std::size_t pos = th1Name.find(etaNames_p->at(etaIter));
      if(pos != std::string::npos){
        etaPos = etaIter;
        break;
      }
    }

    th1Canv_p[dirPos][ptPos]->cd(etaPos+1);

    th1_p[iter]->SetMaximum(dirMax[dirPos][ptPos]);
    th1_p[iter]->SetMinimum(dirMin[dirPos][ptPos]);

    th1_p[iter]->GetXaxis()->CenterTitle();
    th1_p[iter]->GetXaxis()->SetTitleOffset(1.);
    th1_p[iter]->GetXaxis()->SetTitleFont(43);
    th1_p[iter]->GetXaxis()->SetTitleSize(26);
    th1_p[iter]->GetXaxis()->SetNdivisions(505);
    th1_p[iter]->GetYaxis()->CenterTitle();
    th1_p[iter]->GetYaxis()->SetTitleOffset(1.2);
    th1_p[iter]->GetYaxis()->SetNdivisions(505);

    Int_t colPos = 0;
    if(th1Name.find("_Q_") != std::string::npos) colPos = 1;
    else if(th1Name.find("_G_") != std::string::npos) colPos = 2;

    if(th1Name.find("Fit") != std::string::npos){
      th1_p[iter]->SetMarkerColor(qgCol[colPos]);
      th1_p[iter]->SetLineColor(qgCol[colPos]);
      if(!legAdded[0][colPos]){
	meanLeg_p->AddEntry(th1_p[iter], Form("Fit %s %s", meanResStr2[meanResNum].c_str(), qgStr2[colPos].c_str()), "P L");
	legAdded[0][colPos] = true;
      }
    }
    else{
      th1_p[iter]->SetMarkerColor(qgCol[colPos]);
      th1_p[iter]->SetMarkerStyle(24);
      th1_p[iter]->SetLineColor(qgCol[colPos]);

      if(!legAdded[1][colPos]){
	if(isPbPb) meanLeg2_p->AddEntry(th1_p[iter], Form("True %s %s", meanResStr2[meanResNum].c_str(), qgStr2[colPos].c_str()), "P L");
	else meanLeg_p->AddEntry(th1_p[iter], Form("True %s %s", meanResStr2[meanResNum].c_str(), qgStr2[colPos].c_str()), "P L");
	legAdded[1][colPos] = true;
      }
    }

    //    gStyle->SetErrorX(1);

    if(!isDrawn[dirPos][ptPos][etaPos]){
      th1_p[iter]->DrawCopy("E1 P ");
      isDrawn[dirPos][ptPos][etaPos] = true;
    }
    else th1_p[iter]->DrawCopy("E1 P SAME");
    
    if(dirMin[dirPos][ptPos] < 1 && dirMax[dirPos][ptPos] > 1 && strcmp("Res", meanResStr[meanResNum].c_str()) != 0){
      TLine* oneLine_p = new TLine(th1_p[iter]->GetXaxis()->GetXmin(), 1, th1_p[iter]->GetXaxis()->GetXmax(), 1);
      oneLine_p->SetLineStyle(2);
      oneLine_p->DrawClone();
      delete oneLine_p;

      if(!strcmp(inHistName[inHistNum].c_str(), "RecoOverGen")){
	TLine* oneLineUp_p = new TLine(th1_p[iter]->GetXaxis()->GetXmin(), 1.05, th1_p[iter]->GetXaxis()->GetXmax(), 1.05);
	oneLineUp_p->SetLineStyle(2);
	oneLineUp_p->DrawClone();
	delete oneLineUp_p;

	TLine* oneLineDown_p = new TLine(th1_p[iter]->GetXaxis()->GetXmin(), .95, th1_p[iter]->GetXaxis()->GetXmax(), .95);
	oneLineDown_p->SetLineStyle(2);
	oneLineDown_p->DrawClone();
	delete oneLineDown_p;
      }
    }    

    std::string tempPtStr = ptNames_p->at(ptPos);
    tempPtStr.replace(tempPtStr.find("Pt"), 2, "");
    tempPtStr.replace(tempPtStr.find("p"), 1, ".");
    tempPtStr.replace(tempPtStr.find("p"), 1, ".");
    tempPtStr.replace(tempPtStr.find("To"), 2, "<p_{T}<");

    if(etaPos == 0) label_p->DrawLatex(.30, .9, Form("#bf{#color[2]{%s}}", dirNames_p->at(dirPos).c_str()));
    if(etaPos == 1) label_p->DrawLatex(.10, .9, Form("#bf{#color[2]{%s}}", tempPtStr.c_str()));

    std::string tempEtaStr = etaNames_p->at(etaPos);
    tempEtaStr.replace(tempEtaStr.find("Eta"), 3, "");
    tempEtaStr.replace(tempEtaStr.find("p"), 1, ".");
    tempEtaStr.replace(tempEtaStr.find("p"), 1, ".");
    tempEtaStr.replace(tempEtaStr.find("To"), 2, "<#eta<");
    if(etaPos == 0) label_p->DrawLatex(.3, .3, Form("#bf{%s}", tempEtaStr.c_str()));
    else label_p->DrawLatex(.1, .3, Form("#bf{%s}", tempEtaStr.c_str()));
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
    for(Int_t ptIter = 0; ptIter < nPtCanv; ptIter++){
      if(!isPbPb) th1Canv_p[iter][ptIter]->cd();
      else th1Canv_p[iter][ptIter]->cd(3);
      meanLeg_p->Draw("SAME");

      if(plotTrue || !strcmp(inHistName[inHistNum].c_str(), "Eff")){
	if(!isPbPb) th1Canv_p[iter][ptIter]->cd();
	else th1Canv_p[iter][ptIter]->cd(2);
	
	if(isPbPb) meanLeg2_p->Draw("SAME");      
      }
      
      th1Canv_p[iter][ptIter]->Write("", TObject::kOverwrite);
      claverCanvasSaving(th1Canv_p[iter][ptIter], Form("pdfDir/%s", th1Canv_p[iter][ptIter]->GetName()), "pdf");
    }
  }

  outFile_p->Close();
  delete outFile_p;

  delete meanLeg_p;

  delete label_p;

  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t ptIter = 0; ptIter < nPtCanv; ptIter++){
      delete th1Canv_p[iter][ptIter];
    }
  }

  th1Names_p->clear();
  delete th1Names_p;

  dirNames_p->clear();
  delete dirNames_p;

  inFile_p->Close();
  delete inFile_p;

  return;
}








void makeJECPlotMeanRes_Rho(const std::string inFileName, const Int_t inHistNum, const Int_t meanResNum, const Bool_t isPbPb)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nContents = inFile_p->GetListOfKeys()->GetEntries();

  Int_t nTH1Temp = 0;
  std::vector<std::string>* th1Names_p = new std::vector<std::string>;

  Int_t nDirTemp = 0;
  std::vector<std::string>* dirNames_p = new std::vector<std::string>;

  std::vector<std::string>* rhoNames_p = new std::vector<std::string>;
  std::vector<std::string>* rhoLabels_p = new std::vector<std::string>;

  for(Int_t iter = 0; iter < nContents; iter++){
    TString className = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetClassName();
    TString name = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetName();

    if(className.Index("TDirectoryFile") >= 0){
      nDirTemp++;
      dirNames_p->push_back(name.Data());

      TDirectoryFile* tempDir_p = (TDirectoryFile*)inFile_p->Get(name);

      const Int_t nDirContents = tempDir_p->GetListOfKeys()->GetEntries();

      for(Int_t dirIter = 0; dirIter < nDirContents; dirIter++){
	TString name2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetName();
	TString className2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetClassName();

	if(className2.Index("TH1") >= 0){
	  if(name2.Index(Form("%s_", meanResStr[meanResNum].c_str())) >= 0){
	    if(name2.Index("RecoRho") < 0) continue;
	    if(!plotTrue && name2.Index(Form("Fit")) < 0) continue;

	    if(name2.Index("_Q_") >= 0) continue;
	    if(name2.Index("_G_") >= 0) continue;

	    if(name2.Index(inHistName[inHistNum].c_str()) >= 0){
	      nTH1Temp++;
	      th1Names_p->push_back(Form("%s/%s", name.Data(), name2.Data()));

	      std::string tempRhoStr = th1Names_p->at(nTH1Temp-1).substr(th1Names_p->at(nTH1Temp-1).find("RecoRho"), th1Names_p->at(nTH1Temp-1).length());

	      tempRhoStr = tempRhoStr.substr(0, tempRhoStr.find("Over"));

	      Bool_t isRhoVect = false;

	      for(Int_t rhoIter = 0; rhoIter < (Int_t)rhoNames_p->size(); rhoIter++){
		if(!strcmp(rhoNames_p->at(rhoIter).c_str(), tempRhoStr.c_str())){
		  isRhoVect = true;
		  break;
		}
	      }

	      if(!isRhoVect) rhoNames_p->push_back(tempRhoStr);
	      
	      tempRhoStr = th1Names_p->at(nTH1Temp-1).substr(th1Names_p->at(nTH1Temp-1).find("_Rho")+1, th1Names_p->at(nTH1Temp-1).length());

	      tempRhoStr = tempRhoStr.substr(0, tempRhoStr.find("_"));

	      isRhoVect = false;

	      for(Int_t rhoIter = 0; rhoIter < (Int_t)rhoLabels_p->size(); rhoIter++){
		if(!strcmp(rhoLabels_p->at(rhoIter).c_str(), tempRhoStr.c_str())){
		  isRhoVect = true;
		  break;
		}
	      }

	      if(!isRhoVect) rhoLabels_p->push_back(tempRhoStr);

	    }
	  }
	}
      }
    }
  }

  const Int_t nTH1 = nTH1Temp;
  const Int_t nDir = nDirTemp;
  TH1F* th1_p[nTH1];
  
  TCanvas* th1Canv_p[nDir];

  //editing here

  Int_t xBins = 0;
  Int_t yBins = 0;
  Float_t breakBool = false;

  const Int_t nRhoBins = (Int_t)rhoNames_p->size();

  for(Int_t iter = 0; iter < nRhoBins; iter++){
    if(breakBool) break;
    for(Int_t iter2 = 0; iter2 < 4; iter2++){
      xBins = iter + 4 - iter2;
      yBins = iter;

      if(xBins >= 6){
	if(xBins*yBins > nRhoBins && xBins*yBins - nRhoBins < 4){
	  breakBool = true;
	  break;
	}
      }
      else if(xBins*yBins > nRhoBins && xBins*yBins - nRhoBins < ((Float_t)xBins)/2.){
	breakBool = true;
	break;
      }
    }
  }


  Int_t nXPanel = xBins;
  Int_t nYPanel = yBins;

  for(Int_t iter = 0; iter < nDir; iter++){
    th1Canv_p[iter] = new TCanvas(Form("%s_%s%sVPt_RHO_c", dirNames_p->at(iter).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str()), Form("%s_%s%sVPt_RHO_c", dirNames_p->at(iter).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str()), nXPanel*300, nYPanel*325);
    th1Canv_p[iter]->Divide(nXPanel, nYPanel, 0.0, 0.0);
  }


  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  Float_t dirMax[nDir];
  Float_t dirMin[nDir];

  for(Int_t iter = 0; iter < nDir; iter++){
    dirMax[iter] = -1;
    dirMin[iter] = 100;
  }

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);
    th1_p[iter] = (TH1F*)inFile_p->Get(th1Name.c_str());

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    for(Int_t binIter = 0; binIter < th1_p[iter]->GetNbinsX(); binIter++){
      if(th1_p[iter]->GetBinContent(binIter+1) > dirMax[dirPos]) dirMax[dirPos] = th1_p[iter]->GetBinContent(binIter+1);

      if(th1_p[iter]->GetBinContent(binIter+1) < dirMin[dirPos] && th1_p[iter]->GetBinContent(binIter+1) != 0) dirMin[dirPos] = th1_p[iter]->GetBinContent(binIter+1);
    }
  }

  

  for(Int_t iter = 0; iter < nDir; iter++){
    dirMax[iter] += (dirMax[iter] - dirMin[iter])/10.;
    dirMin[iter] -= (dirMax[iter] - dirMin[iter])/10.;

    if(strcmp(meanResStr[meanResNum].c_str(), "Res") != 0){
      if(dirMax[iter] < 1.2) dirMax[iter] = 1.2;
      if(dirMin[iter] > .8) dirMin[iter] = .8;
    }

    dirMax[iter] = setMaxMinNice(dirMax[iter], true);
    dirMin[iter] = setMaxMinNice(dirMin[iter], false);
  }

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    Int_t rhoPos = -1;
    for(Int_t rhoIter = 0; rhoIter < nXPanel*nYPanel; rhoIter++){
      std::size_t pos = th1Name.find(rhoNames_p->at(rhoIter));
      if(pos != std::string::npos){
        rhoPos = rhoIter;
        break;
      }
    }

    th1Canv_p[dirPos]->cd(rhoPos+1);

    th1_p[iter]->SetMaximum(dirMax[dirPos]);
    th1_p[iter]->SetMinimum(dirMin[dirPos]);

    th1_p[iter]->GetXaxis()->CenterTitle();
    th1_p[iter]->GetXaxis()->SetTitleOffset(1.5);
    th1_p[iter]->GetXaxis()->SetTitleFont(43);
    th1_p[iter]->GetXaxis()->SetTitleSize(26);
    th1_p[iter]->GetXaxis()->SetNdivisions(505);
    th1_p[iter]->GetYaxis()->CenterTitle();
    th1_p[iter]->GetYaxis()->SetTitleOffset(2.0);
    th1_p[iter]->GetYaxis()->SetNdivisions(505);

    Int_t colPos = 0;
    if(th1Name.find("_Q_") != std::string::npos) colPos = 1;
    else if(th1Name.find("_G_") != std::string::npos) colPos = 2;

    if(th1Name.find("Fit") != std::string::npos){
      th1_p[iter]->SetMarkerColor(qgCol[colPos]);
      th1_p[iter]->SetLineColor(qgCol[colPos]);
    }
    else{
      th1_p[iter]->SetMarkerColor(qgCol[colPos]);
      th1_p[iter]->SetMarkerStyle(24);
      th1_p[iter]->SetLineColor(qgCol[colPos]);
    }

    th1_p[iter]->DrawCopy("E1 P ");
    
    gPad->SetLogx();

    if(dirMin[dirPos] < 1 && dirMax[dirPos] > 1 && strcmp("Res", meanResStr[meanResNum].c_str()) != 0){
      TLine* oneLine_p = new TLine(th1_p[iter]->GetXaxis()->GetXmin(), 1, th1_p[iter]->GetXaxis()->GetXmax(), 1);
      oneLine_p->SetLineStyle(2);
      oneLine_p->DrawClone();
      delete oneLine_p;

      if(!strcmp(inHistName[inHistNum].c_str(), "RecoOverGen")){
	TLine* oneLineUp_p = new TLine(th1_p[iter]->GetXaxis()->GetXmin(), 1.05, th1_p[iter]->GetXaxis()->GetXmax(), 1.05);
	oneLineUp_p->SetLineStyle(2);
	oneLineUp_p->DrawClone();
	delete oneLineUp_p;

	TLine* oneLineDown_p = new TLine(th1_p[iter]->GetXaxis()->GetXmin(), .95, th1_p[iter]->GetXaxis()->GetXmax(), .95);
	oneLineDown_p->SetLineStyle(2);
	oneLineDown_p->DrawClone();
	delete oneLineDown_p;
      }
    }    


    std::string tempRhoStr = rhoLabels_p->at(rhoPos);
    tempRhoStr.replace(tempRhoStr.find("Rho"), 3, "");
    tempRhoStr.replace(tempRhoStr.find("p"), 1, ".");
    tempRhoStr.replace(tempRhoStr.find("p"), 1, ".");
    tempRhoStr.replace(tempRhoStr.find("to"), 2, "<#rho<");

    if(rhoPos%nXPanel == 0){
      if(rhoPos == 0) label_p->DrawLatex(.30, .8, Form("#bf{#color[2]{%s}}", dirNames_p->at(dirPos).c_str()));
      label_p->DrawLatex(.30, .9, Form("#bf{#color[2]{%s}}", tempRhoStr.c_str()));
    }
    else{
      if(rhoPos == 1){
	if(isPbPb) label_p->DrawLatex(.10, .8, Form("#bf{#color[2]{PYTHIA+HYDJET}}"));
	else label_p->DrawLatex(.10, .8, Form("#bf{#color[2]{PYTHIA}}"));
      }
      label_p->DrawLatex(.10, .9, Form("#bf{#color[2]{%s}}", tempRhoStr.c_str()));
    }
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
    th1Canv_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(th1Canv_p[iter], Form("pdfDir/%s", th1Canv_p[iter]->GetName()), "pdf");
  }

  outFile_p->Close();
  delete outFile_p;
  delete label_p;

  for(Int_t iter = 0; iter < nDir; iter++){
    delete th1Canv_p[iter];
  }

  th1Names_p->clear();
  delete th1Names_p;

  dirNames_p->clear();
  delete dirNames_p;

  inFile_p->Close();
  delete inFile_p;

  return;
}






void makeJECPlotMeanRes_Comp(const std::string inFileName, const Int_t alg1, const Int_t alg2, const Int_t inHistNum, const Int_t meanResNum, const Int_t ptEtaNum, const Bool_t isPbPb)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nContents = inFile_p->GetListOfKeys()->GetEntries();

  Int_t nTH1Temp = 0;
  std::vector<std::string>* th1Names_p = new std::vector<std::string>;

  Int_t nDirTemp = 0;
  std::vector<std::string>* dirNames_p = new std::vector<std::string>;

  for(Int_t iter = 0; iter < nContents; iter++){
    TString className = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetClassName();
    TString name = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetName();

    if(name.Index(jetAlgo[alg1].c_str()) < 0 && name.Index(jetAlgo[alg2].c_str()) < 0) continue;

    if(className.Index("TDirectoryFile") >= 0){
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

	    if(!plotTrue && name2.Index(Form("Fit")) < 0 && strcmp("Eff", inHistName[inHistNum].c_str()) != 0 && strcmp("Fake", inHistName[inHistNum].c_str()) != 0 && strcmp("RecoGenDPhi", inHistName[inHistNum].c_str()) != 0) continue;

	    if(name2.Index("PERP") >= 0) continue;

	    if(name2.Index("_Inc_") < 0) continue;

 	    if(name2.Index("Reco5") >= 0) continue;
	    if(name2.Index("Reco10") >= 0) continue;
	    if(name2.Index("Reco15") >= 0) continue;
	    if(name2.Index("Reco20") >= 0) continue;
	    if(name2.Index("Reco30") >= 0) continue;

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
  
  TCanvas* th1Canv_p;

  Int_t nXPanel = 2;
  Int_t nYPanel = 1;
  if(!isPbPb){
    nXPanel = 1;
    nYPanel = 1;
  }

  th1Canv_p = new TCanvas(Form("%sV%s_%s%s%sc", dirNames_p->at(0).c_str(), dirNames_p->at(1).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str()), Form("%sV%s_%s%s%sc", dirNames_p->at(0).c_str(), dirNames_p->at(1).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str()), nXPanel*300, nYPanel*325);
  th1Canv_p->Divide(nXPanel, nYPanel, 0.0, 0.0);

  Bool_t isDrawn[nXPanel*nYPanel];
  for(Int_t iter2 = 0; iter2 < nXPanel*nYPanel; iter2++){
    isDrawn[iter2] = false;
  }
  
  Bool_t legAdded[2][3] = {{false}, {false}};

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  Float_t dirMax = -1;
  Float_t dirMin = 100;

  Float_t th1XMin = 30;

  TLegend* meanLeg_p;
  if(strcmp(inHistName[inHistNum].c_str(), "RecoGenDPhi") != 0){
    if(isPbPb) meanLeg_p = new TLegend(.30, .68, .90, .95);
    else meanLeg_p = new TLegend(.30, .58, .90, .85);
  }
  else{
    if(isPbPb) meanLeg_p = new TLegend(.30, .38, .90, .95);
    else meanLeg_p = new TLegend(.30, .38, .90, .85);
  }
  meanLeg_p->SetBorderSize(0);
  meanLeg_p->SetFillColor(0);
  meanLeg_p->SetFillStyle(0);
  meanLeg_p->SetTextFont(43);
  meanLeg_p->SetTextSize(18);

  TLegend* meanLeg2_p;
  if(strcmp(inHistName[inHistNum].c_str(), "RecoGenDPhi") != 0) meanLeg2_p = new TLegend(.30, .68, .90, .95);
  else meanLeg2_p = new TLegend(.30, .38, .90, .58); 
  meanLeg2_p->SetBorderSize(0);
  meanLeg2_p->SetFillColor(0);
  meanLeg2_p->SetFillStyle(0);
  meanLeg2_p->SetTextFont(43);
  meanLeg2_p->SetTextSize(18);

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    th1_p[iter] = (TH1F*)inFile_p->Get(th1Name.c_str());

    if(iter == 0 && strcmp("Pt", ptEtaStr2[ptEtaNum].c_str()) != 0) th1XMin = th1_p[iter]->GetBinCenter(1);

    Int_t xMinBin = th1_p[iter]->FindBin(th1XMin);

    if(iter == 0) th1XMin = th1_p[iter]->GetBinLowEdge(xMinBin);

    //insert new xmin
    th1XMin = th1_p[iter]->GetXaxis()->GetXmin();
    th1_p[iter]->SetAxisRange(th1XMin, th1_p[iter]->GetXaxis()->GetXmax(), "X");

    Float_t dirMinBound = 0.5;
    if(!strcmp("Fake", inHistName[inHistNum].c_str())) dirMinBound = 0.0;
    else if(!strcmp("Res", meanResStr[meanResNum].c_str())) dirMinBound = 0.0;

    Float_t dirMaxBound = 1.5;

    for(Int_t binIter = xMinBin-1; binIter < th1_p[iter]->GetNbinsX(); binIter++){
      if(th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1) > dirMax && th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1) < dirMaxBound) dirMax = th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1);

      if(th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1) < dirMin && th1_p[iter]->GetBinContent(binIter+1) != 0 && th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1) > dirMinBound) dirMin = th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1);
    }
  }

  Float_t dirMaxMin = 1.15;
  if(!strcmp("Res", meanResStr[meanResNum].c_str())) dirMaxMin = 0;

  if(dirMax < dirMaxMin) dirMax = dirMaxMin;

  if(!strcmp(inHistName[inHistNum].c_str(), "Eff")) dirMax = 1. + (1-dirMin)/5;
  else{
    dirMax = setMaxMinNice(dirMax, true);
    dirMin = setMaxMinNice(dirMin, false);    
  }
  
  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

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

    th1Canv_p->cd(centPos+1);

    th1_p[iter]->SetMaximum(dirMax);
    th1_p[iter]->SetMinimum(dirMin);

    th1_p[iter]->GetXaxis()->CenterTitle();
    th1_p[iter]->GetXaxis()->SetTitleOffset(th1_p[iter]->GetXaxis()->GetTitleOffset() - .2);
    th1_p[iter]->GetXaxis()->SetTitleFont(43);
    th1_p[iter]->GetXaxis()->SetTitleSize(26);
    th1_p[iter]->GetYaxis()->CenterTitle();
    th1_p[iter]->GetYaxis()->SetTitleOffset(th1_p[iter]->GetYaxis()->GetTitleOffset() - .8);
    
    if(!strcmp(ptEtaStr2[ptEtaNum].c_str(), "Pt")) gPad->SetLogx();
    
    Int_t colPos = 0;
    if(th1Name.find(jetAlgo[alg1].c_str()) != std::string::npos) colPos = 1;
    else if(th1Name.find(jetAlgo[alg2].c_str()) != std::string::npos) colPos = 2;

    Int_t colPos2[3] = {-1, alg1, alg2};

    if(th1Name.find("Fit") != std::string::npos){
      th1_p[iter]->SetMarkerColor(qgCol[colPos]);
      th1_p[iter]->SetLineColor(qgCol[colPos]);
      if(!legAdded[0][colPos]){
	meanLeg_p->AddEntry(th1_p[iter], Form("Fit %s %s", meanResStr2[meanResNum].c_str(), jetAlgo[colPos2[colPos]].c_str()), "P L");
	legAdded[0][colPos] = true;
      }
    }
    else{
      th1_p[iter]->SetMarkerColor(qgCol[colPos]);
      th1_p[iter]->SetMarkerStyle(24);
      th1_p[iter]->SetLineColor(qgCol[colPos]);

      if(!legAdded[1][colPos]){
	if(isPbPb) meanLeg2_p->AddEntry(th1_p[iter], Form("True %s %s", meanResStr2[meanResNum].c_str(), jetAlgo[colPos2[colPos]].c_str()), "P L");
	else meanLeg_p->AddEntry(th1_p[iter], Form("True %s %s", meanResStr2[meanResNum].c_str(), jetAlgo[colPos2[colPos]].c_str()), "P L");
	legAdded[1][colPos] = true;

	std::cout << Form("True %s %s", meanResStr2[meanResNum].c_str(), jetAlgo[colPos2[colPos]].c_str()) << std::endl;
      }
    }

    if((inHistNum == 3 || (inHistNum == 0 && meanResNum == 1) || (inHistNum == 8 && meanResNum == 1)) && ptEtaNum == 0){
      std::cout << "SHITS BROKE" << std::endl;
      std::cout << th1_p[iter]->GetName();
      gStyle->SetOptFit(0);
      th1_p[iter]->GetFunction("f1_p")->SetBit(TF1::kNotDraw);
      std::cout << "YOU GONNA FIX IT?" << std::endl;
    }

    th1_p[iter]->GetYaxis()->SetNdivisions(505);
    th1_p[iter]->GetXaxis()->SetNdivisions(505);

    if(!isDrawn[centPos]){
      th1_p[iter]->DrawCopy("E1 P");
      isDrawn[centPos] = true;
    }
    else th1_p[iter]->DrawCopy("E1 P SAME");
    
    if(dirMin < 1 && dirMax > 1 && strcmp("Res", meanResStr[meanResNum].c_str()) != 0){
      TLine* oneLine_p = new TLine(th1XMin, 1, th1_p[iter]->GetXaxis()->GetXmax(), 1);
      oneLine_p->SetLineStyle(2);
      oneLine_p->DrawClone();
      delete oneLine_p;

      TLine* thirtyLine_p = new TLine(30, dirMin, 30, dirMax);
      thirtyLine_p->SetLineStyle(2);
      thirtyLine_p->DrawClone();
      delete thirtyLine_p;

      if(!strcmp(inHistName[inHistNum].c_str(), "RecoOverGen")){
	TLine* oneLineUp_p = new TLine(th1XMin, 1.05, th1_p[iter]->GetXaxis()->GetXmax(), 1.05);
	oneLineUp_p->SetLineStyle(2);
	oneLineUp_p->DrawClone();
	delete oneLineUp_p;

	TLine* oneLineDown_p = new TLine(th1XMin, .95, th1_p[iter]->GetXaxis()->GetXmax(), .95);
	oneLineDown_p->SetLineStyle(2);
	oneLineDown_p->DrawClone();
	delete oneLineDown_p;
      }
    }    

    if(centPos == 1) label_p->DrawLatex(.10, .9, Form("#bf{#color[2]{|#eta_{jet}|<2.0}}"));
    if(centPos == 2 && !strcmp(ptEtaStr2[ptEtaNum].c_str(), "Eta" )) label_p->DrawLatex(.10, .9, Form("#bf{#color[2]{Gen. p_{T}>30}}"));

    if(isPbPb){
      if(centPos%4 == 0 && centPos < 4) label_p->DrawLatex(.28, .5, centStrings2[centPos].c_str());
      else if(centPos%4 == 0 && centPos > 3) label_p->DrawLatex(.28, .5, centStrings2[centPos].c_str());
      else if(centPos%4 != 0 && centPos < 4) label_p->DrawLatex(.08, .5, centStrings2[centPos].c_str());
      else label_p->DrawLatex(.08, .5, centStrings2[centPos].c_str());
    }
    else{
      if(centPos%4 == 0) label_p->DrawLatex(.78, .9, "PP");
    }
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

  if(!isPbPb) th1Canv_p->cd();
  else th1Canv_p->cd(3);
  meanLeg_p->Draw("SAME");
  
  if(plotTrue || !strcmp(inHistName[inHistNum].c_str(), "Eff") || !strcmp(inHistName[inHistNum].c_str(), "RecoGenDPhi")){
    if(!isPbPb) th1Canv_p->cd();
    else th1Canv_p->cd(4);
    
    if(isPbPb) meanLeg2_p->Draw("SAME");      
  }
  
  th1Canv_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(th1Canv_p, Form("pdfDir/%s", th1Canv_p->GetName()), "pdf");

  outFile_p->Close();
  delete outFile_p;

  delete meanLeg_p;

  delete label_p;

  delete th1Canv_p;

  th1Names_p->clear();
  delete th1Names_p;

  dirNames_p->clear();
  delete dirNames_p;

  inFile_p->Close();
  delete inFile_p;

  return;
}



void makeJECPlotMeanPts(const std::string inFileName, const Int_t inHistNum, const Int_t ptEtaNum, const Int_t qgNum, const Bool_t isPbPb)
{
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

    if(className.Index("TDirectoryFile") >= 0){
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

	  if(name2.Index("PERP") >= 0) continue;

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

  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      std::string centStr = centStrings[centIter];
      if(!isPbPb) centStr = "PP";

      th1Canv_p[iter][centIter] = new TCanvas(Form("%s_%s_MeanPts%s%s%s_c", dirNames_p->at(iter).c_str(), qgStr[qgNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str(), centStr.c_str()), Form("%s_%s_MeanPts%s%s%s_c", dirNames_p->at(iter).c_str(), qgStr[qgNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str(), centStr.c_str()), xBins*300, yBins*325);
      th1Canv_p[iter][centIter]->Divide(xBins, yBins, 0.0, 0.0);
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

    gStyle->SetOptStat(0);

    th1_p[iter]->SetMaximum(dirMax[dirPos][centPos]);
    th1_p[iter]->SetMinimum(dirMin[dirPos][centPos]);

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
    if(ptPos == 1) label_p->DrawLatex(.15, .78, Form("#bf{#color[2]{|#eta_{jet}|<2.0}}"));
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

  return;
}


void makeJECPlotScatter(const std::string inFileName, const Int_t inHistNum, const Int_t ptEtaNum, const Int_t qgNum, const Bool_t isPbPb)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nContents = inFile_p->GetListOfKeys()->GetEntries();

  Int_t nTH2Temp = 0;
  std::vector<std::string>* th2Names_p = new std::vector<std::string>;

  Int_t nDirTemp = 0;
  std::vector<std::string>* dirNames_p = new std::vector<std::string>;

  for(Int_t iter = 0; iter < nContents; iter++){
    TString className = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetClassName();
    TString name = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetName();

    if(className.Index("TDirectoryFile") >= 0){
      nDirTemp++;
      dirNames_p->push_back(name.Data());

      TDirectoryFile* tempDir_p = (TDirectoryFile*)inFile_p->Get(name);

      const Int_t nDirContents = tempDir_p->GetListOfKeys()->GetEntries();
      for(Int_t dirIter = 0; dirIter < nDirContents; dirIter++){
        TString name2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetName();
        TString className2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetClassName();

        if(className2.Index("TH2") >= 0){
	  if(name2.Index(ptEtaStr[ptEtaNum].c_str()) < 0 && name2.Index("RecoVGen") < 0) continue;
	  if(name2.Index(Form("_%s_", qgStr[qgNum].c_str())) < 0) continue;

          if(name2.Index(inHistName[inHistNum].c_str()) >= 0){
            nTH2Temp++;
            th2Names_p->push_back(Form("%s/%s", name.Data(), name2.Data()));
          }
        }
      }
    }    
  }
  
  const Int_t nTH2 = nTH2Temp;
  TH2F* th2_p[nTH2];
  TCanvas* th2Canv_p[nTH2];

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  for(Int_t iter = 0; iter < nTH2; iter++){
    std::string th2Name = th2Names_p->at(iter);

    th2_p[iter] = (TH2F*)inFile_p->Get(th2Name.c_str());

    std::size_t pos = th2Name.find("_h");
    if(pos != std::string::npos) th2Name.replace(pos, pos+2, "_c");

    pos = th2Name.find("/");
    if(pos != std::string::npos) th2Name.replace(0, pos+1, "");

    th2Canv_p[iter] = new TCanvas(th2Name.c_str(), th2Name.c_str(), 700, 700);

    gPad->SetLeftMargin(gPad->GetLeftMargin() - .06);
    gPad->SetBottomMargin(gPad->GetBottomMargin() - .06);
    gPad->SetRightMargin(gPad->GetRightMargin() + .11);

    if(!strcmp(ptEtaStr2[ptEtaNum].c_str(), "Pt")) gPad->SetLogx();
    if(strcmp(ptEtaStr2[ptEtaNum].c_str(), "PtEta") != 0) gPad->SetLogz();
    if(!strcmp(ptEtaStr2[ptEtaNum].c_str(), "PtEta")) gPad->SetLogy();

    if(!strcmp(inHistName[inHistNum].c_str(), "RecoVGen")) gPad->SetLogy();

    th2_p[iter]->GetXaxis()->CenterTitle();
    th2_p[iter]->GetXaxis()->SetTitleOffset(th2_p[iter]->GetXaxis()->GetTitleOffset() - .20);
    th2_p[iter]->GetYaxis()->CenterTitle();
    th2_p[iter]->GetYaxis()->SetTitleOffset(th2_p[iter]->GetYaxis()->GetTitleOffset() - .20);

    th2_p[iter]->DrawCopy("COLZ");

    if(isPbPb){
      pos = th2Name.find("cent");
      if(pos != std::string::npos){
	std::string centStr = th2Name.substr(pos);
	pos = centStr.find("_");
	if(pos != std::string::npos) centStr = centStr.substr(4, pos-4);
	
	pos = centStr.find("to");
	
	if(pos != std::string::npos){
	  if(strcmp(ptEtaStr2[ptEtaNum].c_str(), "PtEta") != 0) label_p->DrawLatex(.65, .85, Form("#bf{#color[2]{%s-%s%%}}", centStr.substr(0, pos).c_str(), centStr.substr(pos+2).c_str()));
	  else label_p->DrawLatex(.65, .85, Form("#bf{%s-%s%%}", centStr.substr(0, pos).c_str(), centStr.substr(pos+2).c_str()));
	}
      }
      label_p->DrawLatex(.65, .70, qgStr2[qgNum].c_str());
    }
    else{
      if(strcmp(ptEtaStr2[ptEtaNum].c_str(), "PtEta") != 0) label_p->DrawLatex(.65, .85, Form("#bf{#color[2]{PP}}"));
      else label_p->DrawLatex(.65, .85, Form("#bf{PP}"));

      label_p->DrawLatex(.65, .70, qgStr2[qgNum].c_str());
    }      

    pos = th2Name.find("ak");
    if(pos != std::string::npos){
      std::string centStr = th2Name.substr(pos);
      pos = centStr.find("_");
      if(pos != std::string::npos) centStr = centStr.substr(0, pos);
      if(pos != std::string::npos) label_p->DrawLatex(.65, .92, Form("#bf{#color[2]{%s}}", centStr.c_str()));
    }

    label_p->DrawLatex(.65, .78, "#bf{#color[2]{|#eta_{jet}|<2.0}}");
    if(!strcmp(ptEtaStr2[ptEtaNum].c_str(), "Eta" )) label_p->DrawLatex(.65, .71, Form("#bf{#color[2]{Gen. p_{T}>30}}"));
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

  for(Int_t iter = 0; iter < nTH2; iter++){
    th2Canv_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(th2Canv_p[iter], Form("pdfDir/%s", th2Canv_p[iter]->GetName()), "pdf");
  }

  outFile_p->Close();
  delete outFile_p;

  delete label_p;

  for(Int_t iter = 0; iter < nTH2; iter++){
    delete th2Canv_p[iter];
  }

  dirNames_p->clear();
  delete dirNames_p;

  th2Names_p->clear();
  delete th2Names_p;

  inFile_p->Close();
  delete inFile_p;

  return;
}


void makeCentVRhoPlot(const std::string inFileName, Bool_t isPbPb)
{
  if(!isPbPb) return;

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TH2F* centVRho_p = (TH2F*)inFile_p->Get("centVRho_h");

  TCanvas* centVRhoCanv_p = new TCanvas("centVRho_c", "centVRho_c", 350, 300);
  centVRho_p->GetXaxis()->SetNdivisions(404);
  centVRho_p->GetYaxis()->SetNdivisions(404);

  centVRho_p->GetXaxis()->SetTitleFont(43);
  centVRho_p->GetXaxis()->SetTitleSize(24);
  centVRho_p->GetXaxis()->SetLabelFont(43);
  centVRho_p->GetXaxis()->SetLabelSize(24);
  centVRho_p->GetXaxis()->CenterTitle();
  centVRho_p->GetYaxis()->SetTitleFont(43);
  centVRho_p->GetYaxis()->SetTitleSize(24);
  centVRho_p->GetYaxis()->SetLabelFont(43);
  centVRho_p->GetYaxis()->SetLabelSize(24);
  centVRho_p->GetYaxis()->SetTitleOffset(1);
  centVRho_p->GetYaxis()->CenterTitle();

  centVRho_p->DrawCopy("COLZ");
  gPad->SetLeftMargin(.20);
  gPad->SetRightMargin(.20);
  gPad->SetBottomMargin(.20);

  gPad->SetLogz();
  gPad->SetLogy();

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(20);

  label_p->DrawLatex(.30, .9, Form("#bf{PYTHIA+HYDJET}"));

  claverCanvasSaving(centVRhoCanv_p, Form("pdfDir/%s", centVRhoCanv_p->GetName()), "pdf");
  delete centVRhoCanv_p;
  inFile_p->Close();
  delete inFile_p;

  return;
}


void makeJECPlot(const std::string inFileName, const Bool_t isPbPb)
{
  //  makeCentVRhoPlot(inFileName, isPbPb);

  //  makeJECPlotMeanRes_Rho(inFileName, nHistName-1+4+10+1, 1, isPbPb);

  //  return;

 for(Int_t iter = 0; iter < nHistName; iter++){
    //    if(strcmp(inHistName[iter].c_str(), "Eff") != 0) continue;

   //   if(iter > 0) continue;

    for(Int_t ptEtaIter = 0; ptEtaIter < nPtEta; ptEtaIter++){
      //      if(ptEtaIter > 0) continue;

      for(Int_t iter2 = 0; iter2 < nMeanRes; iter2++){
	//	if(iter2 == 0) continue;

	if((!strcmp(inHistName[iter].c_str(), "Eff") || !strcmp(inHistName[iter].c_str(), "Fake")) && !strcmp(meanResStr[iter2].c_str(), "Res")) continue;

	if(ptEtaIter < 2) makeJECPlotMeanRes(inFileName, iter, iter2, ptEtaIter, isPbPb);
     }
      
      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	if(!strcmp(inHistName[iter].c_str(), "Fake") && qgIter > 0) continue;

	if(ptEtaIter < 2 && iter != 3) makeJECPlotMeanPts(inFileName, iter, ptEtaIter, qgIter, isPbPb);
	//	makeJECPlotScatter(inFileName, iter, ptEtaIter, qgIter, isPbPb);
      }
      
    }
  }

 
  for(Int_t ptEtaIter = 0; ptEtaIter < nPtEta-1; ptEtaIter++){
    makeJECPlotMeanRes(inFileName, nHistName, 0, ptEtaIter, isPbPb);
    makeJECPlotMeanRes(inFileName, nHistName+2, 0, ptEtaIter, isPbPb);
    makeJECPlotMeanRes(inFileName, nHistName+3, 0, ptEtaIter, isPbPb);

    makeJECPlotMeanRes(inFileName, nHistName, 1, ptEtaIter, isPbPb);
    makeJECPlotMeanRes(inFileName, nHistName+2, 1, ptEtaIter, isPbPb);
    makeJECPlotMeanRes(inFileName, nHistName+3, 1, ptEtaIter, isPbPb);
  }
 
  /*
  makeJECPlotMeanRes(inFileName, nHistName+4, 1, 0, isPbPb);
  makeJECPlotMeanRes(inFileName, nHistName+5, 1, 0, isPbPb);
  makeJECPlotMeanRes(inFileName, nHistName+6, 1, 0, isPbPb);
  makeJECPlotMeanRes(inFileName, nHistName+7, 1, 0, isPbPb);
  makeJECPlotMeanRes(inFileName, nHistName+8, 1, 0, isPbPb);
  makeJECPlotMeanRes(inFileName, nHistName+9, 1, 0, isPbPb);
  makeJECPlotMeanRes(inFileName, nHistName+10, 1, 0, isPbPb);
  makeJECPlotMeanRes(inFileName, nHistName+11, 1, 0, isPbPb);
  makeJECPlotMeanRes(inFileName, nHistName+12, 1, 0, isPbPb);
  makeJECPlotMeanRes(inFileName, nHistName+13, 1, 0, isPbPb);
  */
  for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
    //    makeJECPlotScatter(inFileName, nHistName, 0, qgIter, isPbPb);
    //    makeJECPlotScatter(inFileName, nHistName+1, 0, qgIter, isPbPb);
  }
 
  /*
  for(Int_t mIter = 0; mIter < nMeanRes; mIter++){
    makeJECPlotMeanRes_Cent(inFileName, 0, mIter, isPbPb);
  }
  */

  
  for(Int_t iter = 0; iter < nJetAlgo-1; iter++){
    for(Int_t iter2 = iter+1; iter2 < nJetAlgo; iter2++){
      
      std::cout << "iters: " << iter << ", " << iter2 << std::endl;
      
      for(Int_t ptEtaIter = 0; ptEtaIter < nPtEta-1; ptEtaIter++){
	//makeJECPlotMeanRes_Comp(inFileName, iter, iter2, 3, 0, ptEtaIter, isPbPb);
	makeJECPlotMeanRes_Comp(inFileName, iter, iter2, 0, 0, ptEtaIter, isPbPb);
	makeJECPlotMeanRes_Comp(inFileName, iter, iter2, 0, 1, ptEtaIter, isPbPb);
	
	makeJECPlotMeanRes_Comp(inFileName, iter, iter2, 3, 0, ptEtaIter, isPbPb);
	
	makeJECPlotMeanRes_Comp(inFileName, iter, iter2, 8, 1, ptEtaIter, isPbPb);
      }
    }
  }
  
  return;
}
