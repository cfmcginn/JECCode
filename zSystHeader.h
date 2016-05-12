#ifndef KAYASYSTHEADER_h
#define KAYASYSTHEADER_h

#include "TMath.h"

Float_t returnXJZSystErr(Float_t xJZ, Float_t yVal, Bool_t isPbPb)
{
  Float_t retVal = 0;

  if(isPbPb){
    if(xJZ < 0.4){
      retVal += (.01*yVal)*(.01*yVal);
      retVal += (.001*yVal)*(.001*yVal);
      retVal += (.01*yVal)*(.01*yVal);
      retVal += (.2*yVal)*(.2*yVal);
      retVal += (.03*yVal)*(.03*yVal);
      retVal += (.1*yVal)*(.1*yVal);
      
      retVal = TMath::Sqrt(retVal);
    }
    else if(xJZ < 0.6){
      retVal += (.025*yVal)*(.025*yVal);
      retVal += (.001*yVal)*(.001*yVal);
      retVal += (.03*yVal)*(.03*yVal);
      retVal += (.1*yVal)*(.1*yVal);
      retVal += (.03*yVal)*(.03*yVal);
      retVal += (.05*yVal)*(.05*yVal);
      
      retVal = TMath::Sqrt(retVal);
    }
    else if(xJZ < 1){
      retVal += (.03*yVal)*(.03*yVal);
      retVal += (.01*yVal)*(.01*yVal);
      retVal += (.01*yVal)*(.01*yVal);
      retVal += (.15*yVal)*(.15*yVal);
      retVal += (.01*yVal)*(.01*yVal);
      retVal += (.03*yVal)*(.03*yVal);
      
      retVal = TMath::Sqrt(retVal);
    }
    else if(xJZ < 1.5){
      retVal += (.01*yVal)*(.01*yVal);
      retVal += (.01*yVal)*(.01*yVal);
      retVal += (.01*yVal)*(.01*yVal);
      retVal += (.1*yVal)*(.1*yVal);
      retVal += (.01*yVal)*(.01*yVal);
      retVal += (.03*yVal)*(.03*yVal);
      
      retVal = TMath::Sqrt(retVal);
    }
    else{
      retVal += (.001*yVal)*(.001*yVal);
      retVal += (.03*yVal)*(.03*yVal);
      retVal += (.001*yVal)*(.001*yVal);
      retVal += (.08*yVal)*(.08*yVal);
      retVal += (.001*yVal)*(.001*yVal);
      retVal += (.03*yVal)*(.03*yVal);
      
      retVal = TMath::Sqrt(retVal);
    }
  }
  else{
    if(xJZ < 0.4){
      retVal += (.01*yVal)*(.01*yVal);
      retVal += (.01*yVal)*(.01*yVal);
      retVal += (.01*yVal)*(.01*yVal);
      
      retVal = TMath::Sqrt(retVal);
    }
    else if(xJZ < 0.6){
      retVal += (.01*yVal)*(.01*yVal);
      retVal += (.06*yVal)*(.06*yVal);
      retVal += (.01*yVal)*(.01*yVal);
      
      retVal = TMath::Sqrt(retVal);
    }
    else if(xJZ < 1){
      retVal += (.001*yVal)*(.001*yVal);
      retVal += (.05*yVal)*(.05*yVal);
      retVal += (.01*yVal)*(.01*yVal);
      
      retVal = TMath::Sqrt(retVal);
    }
    else if(xJZ < 1.5){
      retVal += (.001*yVal)*(.001*yVal);
      retVal += (.08*yVal)*(.08*yVal);
      retVal += (.02*yVal)*(.02*yVal);
      
      retVal = TMath::Sqrt(retVal);
    }
    else{
      retVal += (.001*yVal)*(.001*yVal);
      retVal += (.01*yVal)*(.01*yVal);
      retVal += (.01*yVal)*(.01*yVal);
      
      retVal = TMath::Sqrt(retVal);
    }

  }

  return retVal;  
}
  
#endif
