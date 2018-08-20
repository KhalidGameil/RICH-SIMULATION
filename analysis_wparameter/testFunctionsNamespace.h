//# CPUPROFILE=./tmp/main.cpp ./main.cpp
//C++
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <numeric>
#include <stdio.h>
#include <stdlib.h> /* srand, rand */
#include <string>
#include <sys/stat.h>
#include <time.h> /* time */
#include <vector>

//Boost
//include "/usr/local/include/boost/algorithm/string.hpp"

//Root
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGFrame.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TPDF.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TRandom2.h"
#include "TString.h"
#include "TStyle.h"
#include "TVirtualPad.h"

//Classes
#include "aerogelClass.h"
#include "detectedEventClass.h"
 
#include "detectorClass.h"
#include "directionClass.h"
#include "interpolateLikelihoodClass.h"
#include "pmtClass.h"
#include "positionClass.h"
#include "randomEventClass.h"

//NameSpaces
#include "graphFunctionsNamespace.cpp"

double calcLogOfDetectedParticle(std::vector<double> mean, double pmtDN) {
  double ll_event = 0.0;
  for (unsigned long i = 0; i < mean.size(); i++) {
    mean.at(i) = mean.at(i) + pmtDN;
    ll_event = ll_event + -2 * log(1 - TMath::Poisson(0.0, mean.at(i)));
  }
  return ll_event;
}

double calcProbOfDetectedParticle(std::vector<double> mean, double pmtDN) {
  double ll_event = 0.0;
  for (unsigned long i = 0; i < mean.size(); i++) {
    mean.at(i) = mean.at(i) + pmtDN;
    ll_event = ll_event + (1 - TMath::Poisson(0.0, mean.at(i)));
  }
  return ll_event;
}

double calcLogOfNotDetectedParticle(std::vector<double> mean, double pmtDN) {
  double ll_event = 0.0;
  for (unsigned long i = 0; i < mean.size(); i++) {
    mean.at(i) = mean.at(i) + pmtDN;
    ll_event = ll_event + -2 * log(TMath::Poisson(0.0, mean.at(i)));
  }
  return ll_event;
}

double calcProbOfNotDetectedParticle(std::vector<double> mean, double pmtDN) {
  double ll_event = 0.0;
  for (unsigned long i = 0; i < mean.size(); i++) {
    mean.at(i) = mean.at(i) + pmtDN;
    ll_event = ll_event + TMath::Poisson(0.0, mean.at(i));
  }
  return ll_event;
}

/*
 * findVariableGivenPH
 */
std::vector<double> findVariableGivenPH(int pixel, std::vector<std::vector<double>> pH,
                                        double x, double y, double theta, double phi, double beta,
                                        std::vector<bool> variableSwitch,
                                        std::shared_ptr<InterpolateLikelihood> linearInter)

{
  std::vector<double> singleVariableList;
  for (unsigned long i = 0; i < pH.size(); i++) {
    linearInter->findVariableGivenMean(pixel, pH[i][pixel],
                                       x, y, theta, phi, beta,
                                       variableSwitch);
    if (variableSwitch.at(0)) {
      singleVariableList.push_back(x);
    } else if (variableSwitch.at(1)) {
      singleVariableList.push_back(y);
    } else if (variableSwitch.at(2)) {
      singleVariableList.push_back(theta);
    } else if (variableSwitch.at(3)) {
      singleVariableList.push_back(phi);
    } else if (variableSwitch.at(4)) {
      singleVariableList.push_back(beta);
    }
  }
  return singleVariableList;
}

void saveDiffScatterPlot(std::vector<double> ll1, std::vector<double> ll2,
                         std::string particleName, std::string variableNames,
                         bool saveAll = false, std::string xAxisName = "",
                         std::string yAxisName = "") {
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  std::string pathName = graphingAndAnalysisCherenkovPID::saveAllInFolder("logLikelihood", variableNames, saveAll);
  std::string fileName = pathName + "DifferenceScatter" + particleName;
  fileName = fileName + variableNames;

  std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>("Difference Scatter Plot"));
  std::shared_ptr<TGraphErrors> g1(std::make_shared<TGraphErrors>(ll1.size()));

  std::shared_ptr<TCanvas> c2(std::make_shared<TCanvas>("Difference Scatter Histogram"));
  std::shared_ptr<TH1F> h1(std::make_shared<TH1F>("Difference Scatter Histogram", "Difference Scatter Histogram", 100, -5, 5));

  double mean = 0;
  std::vector<double> deltas;
  for (unsigned long i = 0; i < ll1.size(); i++) {
    double delta = ll1.at(i) - ll2.at(i);
    g1->SetPoint(i, i + 1, delta);
    deltas.push_back(delta);
    mean = mean + delta;
  }
  mean = mean / ll1.size();
  c1->cd();
  g1->Draw("AC");
  c1->SetGrid();
  g1->GetXaxis()->SetTitle((xAxisName).c_str());
  g1->GetYaxis()->SetTitle((yAxisName).c_str());
  g1->SetTitle("Difference Scatter Plot");

  TPaveText *ps = new TPaveText(.75, .7, .95, .9, "NDC");
  ;
  ps->SetName("mystats");
  std::string text = "Mean: ";
  text = text + std::to_string(mean) + "s";
  ps->AddText(text.c_str());
  ps->Draw();

  std::string namePng1 = fileName + ".png";
  std::string nameRoot1 = fileName + ".root";
  std::string namePdf1 = fileName + ".pdf";

  c1->Modified();
  c1->Update();
  //gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  c1->SaveAs(namePng1.c_str());
  c1->SaveAs(nameRoot1.c_str());
  c1->SaveAs(namePdf1.c_str());
}
/*
 * saveSumMeanInterpScatterPlot = stupid code to take sum of mean
 * garbage
 */
void saveSumMeanInterpScatterPlot(std::vector<double> variableList,
                                  std::vector<std::vector<double>> pionMeanVariation,
                                  std::vector<bool> variableSwitch,
                                  std::vector<double> variables, bool betaBool,
                                  bool saveAllInSingleFolder = false) {
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  std::string variableName;
  double value = 0;
  graphingAndAnalysisCherenkovPID::getVariableFromSwitch(variables, variableSwitch, value, variableName);

  std::string variableNames = graphingAndAnalysisCherenkovPID::getFileName(betaBool, variables);
  std::string pathName = graphingAndAnalysisCherenkovPID::saveAllInFolder("ScatterPlotSumMean", variableNames, saveAllInSingleFolder);
  std::string fileName = pathName + "DifferenceScatterSumMean" + variableName;
  fileName = fileName + graphingAndAnalysisCherenkovPID::getFileName(betaBool, variables);

  std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>("Difference Scatter Plot"));
  std::shared_ptr<TGraphErrors> gs(std::make_shared<TGraphErrors>(pionMeanVariation.size()));
  for (unsigned long i = 0; i < pionMeanVariation.size(); i++) {
    double tempMeanSum = 0;
    for (unsigned long j = 0; j < pionMeanVariation.at(i).size(); j++) {
      tempMeanSum = tempMeanSum + pionMeanVariation[i][j];
    }
    gs->SetPoint(i, variableList.at(i), tempMeanSum);
  }

  std::string namePng = fileName + ".png";
  std::string nameRoot = fileName + ".root";
  std::string namePdf = fileName + ".pdf";

  gs->Draw("APL");
  gs->GetXaxis()->SetTitle(variableName.c_str());
  gs->GetYaxis()->SetTitle("Sum of Mean Value");
  gs->SetTitle("Difference Sum Scatter Plot");

  c1->Modified();
  c1->Update();
  c1->SaveAs(namePng.c_str());
  c1->SaveAs(nameRoot.c_str());
  c1->SaveAs(namePdf.c_str());
}
