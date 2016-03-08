#ifndef getBkgEstimate_h
#define getBkgEstimate_h

#include "etaPhiFunc.h"
#include "TMath.h"

const Float_t jtPtCut = 30;

const Float_t etaBound = 2.0;
const Float_t r4Rad = 0.4;
const Float_t r4JetArea = TMath::Pi()*r4Rad*r4Rad;

void getBkgEstimate(Int_t algR, Int_t nPart, Float_t pt[], Float_t phi[], Float_t eta[], Int_t nJt, Float_t jtPt[], Float_t jtPhi[], Float_t jtEta[], Float_t &bkgPerJetCone, Float_t &area)
{
  Int_t jtMult = 0;
  Float_t ptSum = 0;

  const Float_t r4Rad = algR/10.;
  const Float_t r4JetArea = TMath::Pi()*r4Rad*r4Rad;

  for(Int_t iter = 0; iter < nJt; iter++){
    if(jtPt[iter] < jtPtCut) continue;
    if(TMath::Abs(jtEta[iter]) > etaBound) continue;

    jtMult++;
  }

  for(Int_t iter = 0; iter < nPart; iter++){
    if(TMath::Abs(eta[iter]) > etaBound) continue;
    Bool_t skipBool = false;

    for(Int_t jtIter = 0; jtIter < nJt; jtIter++){
      if(jtPt[jtIter] < jtPtCut) continue;
      if(TMath::Abs(jtEta[jtIter]) > etaBound) continue;
      if(jtIter > 1 && jtPt[jtIter] < jtPt[0]*.5 && jtPt[jtIter] < jtPt[1]*.5) continue;

      Float_t dR = getDR(jtEta[jtIter], jtPhi[jtIter], eta[iter], phi[iter]);

      if(dR <= r4Rad){
	skipBool = true;
	break;
      }
    }
    
    if(skipBool) continue;
    
    ptSum += pt[iter];
  }

  area = 2*etaBound*2*TMath::Pi() - jtMult*r4JetArea;

  bkgPerJetCone = ptSum/area;

  return;
}

void getBkgEstimate(Int_t algR, std::vector<float>* pt_p, std::vector<float>* phi_p, std::vector<float>* eta_p, Int_t nJt, Float_t jtPt[], Float_t jtPhi[], Float_t jtEta[], Float_t &bkgPerJetCone, Float_t &area, std::vector<Bool_t>* isBkg_p)
{
  Int_t jtMult = 0;
  Float_t ptSum = 0;

  const Float_t r4Rad = algR/10.;
  const Float_t r4JetArea = TMath::Pi()*r4Rad*r4Rad;

  for(Int_t iter = 0; iter < nJt; iter++){
    if(jtPt[iter] < jtPtCut) continue;
    if(TMath::Abs(jtEta[iter]) > etaBound) continue;
    //    std::cout << "Using jet: " << jtPt[iter] << ", " << jtPhi[iter] << ", " << jtEta[iter] << std::endl;

    jtMult++;
  }

  for(Int_t iter = 0; iter < (Int_t)pt_p->size(); iter++){
    if(TMath::Abs(eta_p->at(iter)) > etaBound) continue;
    Bool_t skipBool = false;

    for(Int_t jtIter = 0; jtIter < nJt; jtIter++){
      if(jtPt[jtIter] < jtPtCut) continue;
      if(TMath::Abs(jtEta[jtIter]) > etaBound) continue;

      Float_t dR = getDR(jtEta[jtIter], jtPhi[jtIter], eta_p->at(iter), phi_p->at(iter));

      if(dR <= r4Rad){
	skipBool = true;
	break;
      }
    }
    
    if(skipBool){
      isBkg_p->push_back(false);
      continue;
    }
    
    ptSum += pt_p->at(iter);
    isBkg_p->push_back(true);
  }

  area = 2*etaBound*2*TMath::Pi() - jtMult*r4JetArea;

  bkgPerJetCone = ptSum/area;

  return;
}

#endif
