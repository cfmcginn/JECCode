#ifndef JECCONFIGPARSER_H
#define JECCONFIGPARSER_H

#include <iostream> 
#include <string>
#include <vector>
#include <fstream>

#include "include/checkMakeDir.h"

//root headers

class jecConfigParser{
 private:
  const static unsigned int nMaxPtHat = 10;
  const std::string txtStr = ".txt";
  const std::string numStr = "0123456789";
  const std::string commaStr = ",";

  const static unsigned int nValidConfigVals = 4;
  const std::string validConfigVals[nValidConfigVals] = {"NPTHAT", "PTHAT", "INPUT", "ISPBPB"};

  unsigned int nPthats = 0;
  std::vector<unsigned int> pthats;
  std::vector<std::string> inputStrings;
  bool isPbPb = false;

 public:
  jecConfigParser();
  jecConfigParser(const std::string);
  bool SetConfigParser(const std::string);
  unsigned int GetNPthats();
  void PrintPthats();
};

jecConfigParser::jecConfigParser()
{
  nPthats = 0;
  return;
}

jecConfigParser::jecConfigParser(const std::string inConfigFile)
{
  if(SetConfigParser(inConfigFile)) return;
  else{
    std::cout << "Initializing jecConfigParser failed. nPthats set to 0." << std::endl;
    nPthats = 0;
  }

  return;
}

bool jecConfigParser::SetConfigParser(const std::string inConfigFile)
{
  std::vector<std::string> lines;
  if(!checkFile(inConfigFile)){
    std::cout << "Input jecConfig txt file, \'" << inConfigFile << "\', is not valid file. return false" << std::endl;
    return false;
  }  
  else if(inConfigFile.size() < txtStr.size()){
    std::cout << "Input jecConfig txt file, \'" << inConfigFile << "\', doesn't end in \'.txt\'. return false" << std::endl;
    return false;
  }
  else if(inConfigFile.substr(inConfigFile.size()-txtStr.size(), txtStr.size()).find(txtStr) == std::string::npos){
    std::cout << "Input jecConfig txt file, \'" << inConfigFile << "\', doesn't end in \'.txt\'. return false" << std::endl;
    return false;
  }

  std::ifstream inFile(inConfigFile.c_str());
  std::string tempStr;
  while(std::getline(inFile, tempStr)){
    if(tempStr.size() == 0) continue; //skip the empties
    if(tempStr.substr(0, 1).find("#") != std::string::npos) continue; //comment handling

    while(tempStr.find(" ") != std::string::npos){
      tempStr.replace(tempStr.find(" "), 1, "");
    }

    //skip invalid inputs and print skip
    bool isValid = false;
    for(unsigned int iter = 0; iter < nValidConfigVals; iter++){
      if(tempStr.substr(0, validConfigVals[iter].size()).find(validConfigVals[iter]) != std::string::npos && tempStr.find("=") != std::string::npos){
	isValid = true;
	break;
      }
    }
    

    if(!isValid){
      std::cout << "Config line \'" << tempStr << "\' is invalid, will be skipped. Please use valid inputs: ";
      for(unsigned int iter = 0; iter < nValidConfigVals; iter++){
	if(iter < nValidConfigVals-1) std::cout << validConfigVals[iter] << ", ";
	else std::cout << validConfigVals[iter] << ", and \'=\'." << std::endl;
      }
      continue;
    }

    lines.push_back(tempStr);
  }

  //begin setting class params
  const unsigned int nLines = lines.size();  
  for(unsigned int iter = 0; iter < nLines; iter++){
    tempStr = lines.at(iter);
    std::string valStr = tempStr.substr(tempStr.find("=")+1, tempStr.size() - tempStr.find("=")+1);
    tempStr.replace(tempStr.find("="), tempStr.size() - tempStr.find("="), "");

    if(tempStr.substr(0, validConfigVals[0].size()).find(validConfigVals[0]) != std::string::npos){
      if(valStr.find("-") != std::string::npos){ // check if negative
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be non-negative. Setting to 0" << std::endl;
	nPthats = 0;
	continue;
      }
      else if(valStr.find(".") != std::string::npos){ // check if integer
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be integer (no \'.\'). Setting to 0" << std::endl;
        nPthats = 0;
	continue;
      }
      
      bool isNum = true; // check if number
      for(unsigned int iter2 = 0; iter2 < valStr.size(); iter2++){
	if(numStr.find(valStr.at(iter2)) == std::string::npos){
	  isNum = false;
	  break;
	}
      }
      if(!isNum){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be pure number from \'" << numStr << "\'. Setting to 0" << std::endl;
	nPthats = 0;
	continue;
      }
      // setting
      nPthats = (unsigned int)std::stoi(valStr);
    
      // finally check its below class cap
      if(nPthats > nMaxPtHat){
	std::cout << tempStr << " value \'" << valStr << "\' is greater than class cap \'" << nMaxPtHat << "\'. Consider raising cap or using fewer pthats. Setting to 0." << std::endl;
	nPthats = 0;
      }
    }
    if(tempStr.substr(0, validConfigVals[1].size()).find(validConfigVals[1]) != std::string::npos){
      while(valStr.find(",,") != std::string::npos){
	valStr.replace(valStr.find(",,"), 2, ",");
      }

      if(valStr.find("-") != std::string::npos){ // check if negative
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be non-negative. Return empty." << std::endl;
	continue;
      }
      else if(valStr.find(".") != std::string::npos){ // check if integer
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be integer (no \'.\'). Return empty" << std::endl;
	continue;
      }
      
      bool isNum = true; // check if number
      for(unsigned int iter2 = 0; iter2 < valStr.size(); iter2++){
	if(numStr.find(valStr.at(iter2)) == std::string::npos && commaStr.find(valStr.at(iter2)) == std::string::npos){
	  isNum = false;
	  break;
	}
      }
      if(!isNum){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be pure number from \'" << numStr << "\' or \',\'. Return empty" << std::endl;
	continue;
      }
      // setting
      while(valStr.find(",") != std::string::npos){
	pthats.push_back((unsigned int)std::stoi(valStr.substr(0, valStr.find(","))));
	valStr.replace(0, valStr.find(",")+1, "");
      }
      if(valStr.size() != 0) pthats.push_back((unsigned int)std::stoi(valStr));
    } 
  }

  if(nPthats == 0){
    pthats.clear();
    inputStrings.clear();
    isPbPb = false;
    return false;
  }
  if(pthats.size() != nPthats){
    std::cout << "NPTHAT, \'" << nPthats << "\', does not match PTHAT size, \'" << pthats.size() << "\'. Return false." << std::endl;
    nPthats = 0;
    pthats.clear();
    inputStrings.clear();
    isPbPb = false;
    return false;
  }
  /*
  if(inputStrings.size() != nPthats){
    nPthats = 0;
    pthats.clear();
    inputStrings.clear();
    isPbPb = false;
    return false;
  }
  */

  return true;
}

unsigned int jecConfigParser::GetNPthats(){return nPthats;}

void jecConfigParser::PrintPthats()
{
  std::cout << "Pthats: ";

  for(unsigned int iter = 0; iter < nPthats; iter++){
    if(iter < nPthats-1) std::cout << pthats.at(iter) << ", ";
    else std::cout << pthats.at(nPthats-1) << ".";
  }

  std::cout << std::endl;

  return;
}

#endif
