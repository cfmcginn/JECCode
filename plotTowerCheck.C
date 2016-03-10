#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TKey.h"
#include "TLatex.h"
#include "TDatime.h"
#include "TLegend.h"

#include <iostream>

void claverCanvasSaving(TCanvas* c, TString s,TString format="gif"){
  TDatime* date = new TDatime();
  c->SaveAs(Form("%s_%d.%s",s.Data(),date->GetDate(), format.Data()));
  return;
}

void plotTowerCheck_Mean(const std::string inFileName, Bool_t isPbPb)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nContents = inFile_p->GetListOfKeys()->GetEntries();

  Int_t nTH1Temp = 0;
  std::vector<std::string>* th1Names_p = new std::vector<std::string>;

  for(Int_t iter = 0; iter < nContents; iter++){
    TString name = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetName();
    TString className = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetClassName();

    if(className.Index("TH1") < 0) continue;
    //    if(name.Index("_Mean_") >= 0) std::cout << name << std::endl;
    if(name.Index("_Mean_") < 0) continue;

    nTH1Temp++;
    th1Names_p->push_back(name.Data());
  }

  const Int_t nTH1 = nTH1Temp;
  TH1F* th1_p[nTH1];
  Int_t centBound[nTH1];
  Int_t centBound2[nTH1];
  Int_t centPos[nTH1];

  Float_t plotMax = -1;
  Float_t plotMin = 100;

  for(Int_t iter = 0; iter < nTH1; iter++){
    th1_p[iter] = (TH1F*)inFile_p->Get(th1Names_p->at(iter).c_str());

    th1_p[iter]->GetXaxis()->SetTitleFont(43);
    th1_p[iter]->GetXaxis()->SetTitleSize(24);
    th1_p[iter]->GetXaxis()->SetTitleOffset(1);
    th1_p[iter]->GetXaxis()->CenterTitle();

    th1_p[iter]->GetYaxis()->SetTitleFont(43);
    th1_p[iter]->GetYaxis()->SetTitleSize(24);
    th1_p[iter]->GetYaxis()->SetTitleOffset(1.4);
    th1_p[iter]->GetYaxis()->CenterTitle();

    if(th1Names_p->at(iter).find("Vs") != std::string::npos){
      th1_p[iter]->SetMarkerColor(kRed);
      th1_p[iter]->SetLineColor(kRed);
    }
    else{
      th1_p[iter]->SetMarkerColor(kBlue);
      th1_p[iter]->SetLineColor(kBlue);
    }

    for(Int_t binIter = 0; binIter < th1_p[iter]->GetNbinsX(); binIter++){
      if(th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1) > plotMax) plotMax = th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1);

      if(th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1) < plotMin && th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1) > 0) plotMin = th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1);
    }

    centBound[iter] =  std::stoi(th1Names_p->at(iter).substr(th1Names_p->at(iter).find("cent")+4, th1Names_p->at(iter).find("to") - (th1Names_p->at(iter).find("cent")+4)));
    centBound2[iter] = std::stoi(th1Names_p->at(iter).substr(th1Names_p->at(iter).find("to")+2, th1Names_p->at(iter).substr(th1Names_p->at(iter).find("to")+2, th1Names_p->at(iter).length()).find("_") - (th1Names_p->at(iter).find("to")+2)));
    centPos[iter] = iter;

    for(Int_t iter2 = 0; iter2 < iter; iter2++){
      if(centBound[iter] > centBound[iter2]){
	Int_t tempCentPos = centPos[iter];
	centPos[iter] = centPos[iter2];
	centPos[iter2] = tempCentPos;
      }
    }
  }


  for(Int_t iter = 0; iter < nTH1; iter++){
    std::cout << "Order: " << centBound[iter] << ", " << centBound2[iter] << ", " << centPos[iter] << std::endl;
  }

  TCanvas* outCanv_p = new TCanvas(Form("jtTowerOverRawVRawPt_Mean_c"), Form("jtTowerOverRawVRawPt_Mean_c"), 4*300, 1*325);
  outCanv_p->Divide(4, 1, 0.0, 0.0);

  Bool_t isDrawn[4] = {false, false, false, false};

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  TLegend* meanLeg_p = new TLegend(.50, .68, .90, .95);
  meanLeg_p->SetBorderSize(0);
  meanLeg_p->SetFillColor(0);
  meanLeg_p->SetFillStyle(0);
  meanLeg_p->SetTextFont(43);
  meanLeg_p->SetTextSize(18);

  
  
  for(Int_t iter = 0; iter < nTH1; iter++){
    outCanv_p->cd(centPos[iter]/2+1);
    th1_p[iter]->SetMaximum(plotMax);
    th1_p[iter]->SetMinimum(plotMin);
    if(isDrawn[centPos[iter]/2]) th1_p[iter]->Draw("E1 P SAME");
    else th1_p[iter]->Draw("E1 P");

    if(centPos[iter]/2 == 0){
      if(th1Names_p->at(iter).find("Vs") == std::string::npos) meanLeg_p->AddEntry(th1_p[iter], "Unsub. Towers", "P L");
      else meanLeg_p->AddEntry(th1_p[iter], "Vs sub. Towers", "P L");
    }

    isDrawn[centPos[iter]/2] = true;

    gPad->SetLogx();

    if(centPos[iter]/2 == 0) label_p->DrawLatex(.6, .65, Form("%d-%d%%", centBound[iter], centBound2[iter]));
    else label_p->DrawLatex(.6, .65, Form("%d-%d%%", centBound[iter], centBound2[iter]));
  }

  delete label_p;

  std::string outName = inFileName;
  const std::string inString = "_HIST.root";
  const std::string outString = "_PLOT.root";
  std::size_t strIndex = 0;

  strIndex = outName.find(inString);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, inString.length(), outString);
  }
  TFile* outFile_p = new TFile(outName.c_str(), "UPDATE");
  outCanv_p->cd(1);
  meanLeg_p->Draw("SAME");
  outCanv_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(outCanv_p, Form("pdfDir/%s", outCanv_p->GetName()), "pdf");
  delete meanLeg_p;
  outFile_p->Close();
  delete outFile_p;

  delete outCanv_p;

  th1Names_p->clear();
  delete th1Names_p;

  inFile_p->Close();
  delete inFile_p;

  return;
}



void plotTowerCheck_TH2(const std::string inFileName, Bool_t isPbPb)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nContents = inFile_p->GetListOfKeys()->GetEntries();

  Int_t nTH2Temp = 0;
  std::vector<std::string>* th2Names_p = new std::vector<std::string>;

  for(Int_t iter = 0; iter < nContents; iter++){
    TString name = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetName();
    TString className = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetClassName();

    if(className.Index("TH2") < 0) continue;

    nTH2Temp++;
    th2Names_p->push_back(name.Data());
  }

  const Int_t nTH2 = nTH2Temp;
  TH2F* th2_p[nTH2];
  Int_t centBound[nTH2];
  Int_t centBound2[nTH2];
  Int_t centPos[nTH2];

  TCanvas* outCanv_p[nTH2];

  for(Int_t iter = 0; iter < nTH2; iter++){
    th2_p[iter] = (TH2F*)inFile_p->Get(th2Names_p->at(iter).c_str());

    th2_p[iter]->GetXaxis()->SetTitleFont(43);
    th2_p[iter]->GetXaxis()->SetTitleSize(24);
    th2_p[iter]->GetXaxis()->SetTitleOffset(1);
    th2_p[iter]->GetXaxis()->CenterTitle();

    th2_p[iter]->GetYaxis()->SetTitleFont(43);
    th2_p[iter]->GetYaxis()->SetTitleSize(24);
    th2_p[iter]->GetYaxis()->SetTitleOffset(1.4);
    th2_p[iter]->GetYaxis()->CenterTitle();

    outCanv_p[iter] = new TCanvas(Form("%s_c", th2Names_p->at(iter).c_str()), Form("%s_c", th2Names_p->at(iter).c_str()), 600, 650);
  }

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);
  
  for(Int_t iter = 0; iter < nTH2; iter++){
    outCanv_p[iter]->cd();
    th2_p[iter]->Draw("COLZ");

    gPad->SetLogx();
    gPad->SetLogz();

    gPad->SetLeftMargin(.14);
    gPad->SetRightMargin(.14);
    gPad->SetBottomMargin(.14);

    std::string labelStr = th2Names_p->at(iter).substr(th2Names_p->at(iter).find("cent")+4, th2Names_p->at(iter).length());
    labelStr = labelStr.substr(0, labelStr.find("_"));
    labelStr.replace(labelStr.find("to"), 2, "-");

    label_p->DrawLatex(.6, .65, Form("%s%%", labelStr.c_str()));
  }

  delete label_p;

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
    outCanv_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(outCanv_p[iter], Form("pdfDir/%s", outCanv_p[iter]->GetName()), "pdf");
  }
  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < nTH2; iter++){
    delete outCanv_p[iter];
  }

  th2Names_p->clear();
  delete th2Names_p;

  inFile_p->Close();
  delete inFile_p;

  return;
}
