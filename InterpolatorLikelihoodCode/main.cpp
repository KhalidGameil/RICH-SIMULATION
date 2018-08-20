//# CPUPROFILE=./tmp/main.cpp ./main.cpp
// C++
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <numeric>
#include <string>
#include <sys/stat.h>
#include <vector>

// Boost
// include "/usr/local/include/boost/algorithm/string.hpp"

// Root
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGFrame.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TRandom2.h"
#include "TString.h"
#include "TStyle.h"
#include "TVirtualPad.h"

// Classes
#include "interpolateLikelihoodClass.h"

std::vector<double> testFunction(std::vector<double> v, int pixels = 256) {
  std::vector<double> mean;
  for (int i = 0; i < pixels; i++) {
    mean.push_back(std::pow(3 * v.at(0), 2) * i + std::pow(i * v.at(1), 2) +
                   15 * v.at(0) * v.at(1) + v.at(2) + v.at(3) + v.at(4) +
                   std::pow(v.at(3), 2) + std::sqrt(v.at(4)) +
                   v.at(3) * v.at(4) * v.at(2));
    // mean.push_back(15*v.at(0) + 12*v.at(1) +std::pow(v.at(0),2) +
    // v.at(0)*v.at(1));
  }
  return mean;
}

int main(int argc, char **argv) {

  int interpolateMeanLinearly = 0;
  int interpolationTriCubicXYTheta = 1;
  int interpolateMeansLinearlyWithSpline = 2;

  std::cout << argc << std::endl;
  std::clock_t begin = std::clock();

  std::string filename = "MeanInterpolationDataBank_database_2018Jan26.txt";
  std::shared_ptr<InterpolateLikelihood> linearInter(
      std::make_shared<InterpolateLikelihood>(filename));
  //linearInter->readDataBaseFileMMapped(filename);
  std::clock_t end = std::clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout << "TIME THE CODE TOOK:" << elapsed_secs << std::endl;
  //linearInter->printVariablLists();

  // std::shared_ptr<InterpolateLikelihood>
  // linearInter(std::make_shared<InterpolateLikelihood>());
  
  std::vector<double> variables;
  variables.reserve(5);
  variables.push_back(0.098516);  // x
  variables.push_back(0.098516);  // y
  variables.push_back(0.193103);    // theta
  variables.push_back(0.5);    // phi
  variables.push_back(0.989); // beta

  std::vector<double> meanInterpBetaSplinePion;
  std::vector<bool> variableSwitchSpline;
  variableSwitchSpline.push_back(false);
  variableSwitchSpline.push_back(false);
  variableSwitchSpline.push_back(false);
  variableSwitchSpline.push_back(false);
  variableSwitchSpline.push_back(true);

  linearInter->getSurroundingPointsIndex(variables);

  
  meanInterpBetaSplinePion = linearInter->interpolateMeansLinearlyWithPolynomialApproximation(variables, variableSwitchSpline);
  std::cout << "MEAN BETA: " << std::endl;
  for (int i = 0; i < meanInterpBetaSplinePion.size(); i++) {
    std::cout << i << "   " << meanInterpBetaSplinePion[i] << std::endl;
  }
}
