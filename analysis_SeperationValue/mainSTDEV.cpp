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
#include "graphFunctionsNamespace.h"

//Geometric ToolBox
#include <Mathematics/GteIntpTricubic3.h>

double getSTDEV(std::vector<double> v) {
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  double mean = sum / v.size();

  double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
  double stdev = std::sqrt(sq_sum / v.size() - mean * mean);
  return stdev;
}

int main(int argc, char **argv) {

  int interpolateMeanLinearly = 0;
  int interpolationTriCubicXYTheta = 1;
  int interpolateMeansLinearlyWithSpline = 2;

  //Mass
  double pionMass = 137.2735; //MeV/c^2
  double kaonMass = 495.665;  //MeV/c^2

  //Position
  double aeroWidth = .12;     //m
  double aeroLength = .12;    //m
  double aeroThickness = .02; //m

  //Refractive Index
  double n1 = 1.010;
  double n2 = 1.012;

  /*Aerogel Parameters*/
  std::shared_ptr<Aerogel> aero1(std::make_shared<Aerogel>(aeroWidth, aeroLength, aeroThickness, n1));
  std::shared_ptr<Aerogel> aero2(std::make_shared<Aerogel>(aeroWidth, aeroLength, aeroThickness, n2));
  std::vector<std::shared_ptr<Aerogel>> aerogels;
  aerogels.push_back(aero1);
  aerogels.push_back(aero2);

  /*PMT Parameters*/
  double multiplePMTs = 4;                       // THIS IS A HACK IDEALLY IT SHOULD BE THAT YOU MAKE MULTIPLE PMT Objects //#CHANGE=1
  double pmtWidth = multiplePMTs * 48.5 / 1000;  //m
  double pmtLength = multiplePMTs * 48.5 / 1000; //m
  int pixelN = std::pow(multiplePMTs * 8, 2);
  std::shared_ptr<PMT> pmt(std::make_shared<PMT>(pmtWidth, pmtLength, pixelN));

  std::string filename = "MeanInterpolationDataBank_10Juluy2017_p1024.txt";
  std::shared_ptr<InterpolateLikelihood> linearInter(std::make_shared<InterpolateLikelihood>(filename));
  //linearInter->readDataBaseFileFast(filename);

  std::shared_ptr<Detector> detector(std::make_shared<Detector>(aerogels, pmt, .09));
  int eventCount = 10000;

  //Momentum
  double startMomentum = 1500;
  double endMomentum = 10000;
  double momentumCount = 1;
  double stepMomentum = graphingAndAnalysisCherenkovPID::calcStepSize(startMomentum, endMomentum, momentumCount);

  //Beta
  double startBeta = .9995; //FORCES TO BE SHOWN
  double endBeta = 1.999;
  double betaCount = 1;
  double stepBeta = graphingAndAnalysisCherenkovPID::calcStepSize(startBeta, endBeta, betaCount);
  bool betaBool = false;
  std::cout << "BETA:" << stepBeta << std::endl;

  //Phi
  double startPhi = 0;
  double endPhi = TMath::Pi() / 4; //2*TMath::Pi()/4;//radians => max = 2*PI
  double phiCount = 1;
  double stepPhi = graphingAndAnalysisCherenkovPID::calcStepSize(startPhi, endPhi, phiCount);
  std::cout << "PHI:" << stepPhi << std::endl;

  //Theta
  double startTheta = 0; //-TMath::Pi()/4;//radians
  double endTheta = 0.4; //TMath::Pi()/4;//radians
  double thetaCount = 1;
  double stepTheta = graphingAndAnalysisCherenkovPID::calcStepSize(startTheta, endTheta, thetaCount);
  std::cout << "THETA " << stepTheta << std::endl;

  //Position
  double beamx = 0.097;               //m //#CHANGE=.06
  double beamy = 0.097;               //m //#CHANGE=.06
  double pixelDist = 48.5 / (8000.0); //dimensions of a pixel
  int beamCountX = 1;
  int beamCountY = 1;
  double beamStepX = graphingAndAnalysisCherenkovPID::calcStepSize(beamx, beamx + pixelDist, beamCountX);
  double beamStepY = graphingAndAnalysisCherenkovPID::calcStepSize(beamy, beamy + pixelDist, beamCountY);
  double beamz = 0.00; //m //#CHANGE=.00

  double darknoise = pmt->getDarkNoise();

  bool saveAllInSingleFolder = true;
  bool completed = false;

  double startI = startMomentum;
  double endI = endMomentum;
  double stepI = stepMomentum;

  if (betaBool) {
    startI = startBeta;
    endI = endBeta;
    stepI = stepBeta;
  }
  gROOT->SetBatch(kTRUE);
  std::vector<double> xVariables;
  std::vector<double> yVariables;
  std::vector<double> thetaVariables;
  std::vector<double> phiVariables;
  std::vector<double> betaMomentumVariables;

  xVariables.push_back(beamx);
  xVariables.push_back(beamx + pixelDist);
  xVariables.push_back(beamStepX);

  yVariables.push_back(beamy);
  yVariables.push_back(beamy + pixelDist);
  yVariables.push_back(beamStepY);

  thetaVariables.push_back(startTheta);
  thetaVariables.push_back(endTheta);
  thetaVariables.push_back(stepTheta);

  phiVariables.push_back(startPhi);
  phiVariables.push_back(endPhi);
  phiVariables.push_back(stepPhi);

  if (betaBool) {
    betaMomentumVariables.push_back(startBeta);
    betaMomentumVariables.push_back(endBeta);
    betaMomentumVariables.push_back(stepBeta);
  } else {
    betaMomentumVariables.push_back(startMomentum);
    betaMomentumVariables.push_back(endMomentum);
    betaMomentumVariables.push_back(stepMomentum);
  }

  if (argc > 1) {
    //std::cout << argv[1] << std::endl;
    std::string filename(argv[1]);
    graphingAndAnalysisCherenkovPID::setVariablesFromJobFiles(filename,
                                                              &xVariables,
                                                              &yVariables,
                                                              &thetaVariables,
                                                              &phiVariables,
                                                              &betaMomentumVariables);
  }

  for (double x = xVariables.at(0); x <= xVariables.at(1); x = x + xVariables.at(2)) {
    for (double y = yVariables.at(0); y <= yVariables.at(1); y = y + yVariables.at(2)) {
      for (double i = betaMomentumVariables.at(0); i <= betaMomentumVariables.at(1); i = i + betaMomentumVariables.at(2)) //I index can either be momentum or beta
      {
        for (double t = thetaVariables.at(0); t <= thetaVariables.at(1); t = t + thetaVariables.at(2)) {
          for (double p = phiVariables.at(0); p <= phiVariables.at(1); p = p + phiVariables.at(2)) {

            //std::cout << "Test Variable [beta/x/y/theta/phi/none]?:";
            std::string variableLikelihoodName = "none";
            //std::cin >> variableLikelihoodName;

            if (betaBool == true) {
              i = startBeta;
            } else {
              i = startMomentum;
            }

            std::clock_t begin = std::clock();

            std::shared_ptr<Position> pos(std::make_shared<Position>(x, y, beamz));
            std::shared_ptr<Direction> dir(std::make_shared<Direction>(t, p));
            std::shared_ptr<Beam> beamPion;
            std::shared_ptr<Beam> beamKaon;

            if (betaBool == true) {
              beamPion = std::make_shared<Beam>(i, dir, pos);
              beamKaon = std::make_shared<Beam>(i, dir, pos);
              beamPion->setParticleMass(pionMass);
              beamKaon->setParticleMass(kaonMass);
            } else {
              beamPion = std::make_shared<Beam>(i, pionMass, dir, pos);
              beamKaon = std::make_shared<Beam>(i, kaonMass, dir, pos);
            }

            std::vector<double> variablesPion;
            variablesPion.reserve(5);
            variablesPion.push_back(x);                   //x
            variablesPion.push_back(y);                   //y
            variablesPion.push_back(t);                   //theta
            variablesPion.push_back(p);                   //phi
            variablesPion.push_back(beamPion->getBeta()); //beta

            std::vector<double> variablesKaon;
            variablesKaon.reserve(5);
            variablesKaon.push_back(x);                   //x
            variablesKaon.push_back(y);                   //y
            variablesKaon.push_back(t);                   //theta
            variablesKaon.push_back(p);                   //phi
            variablesKaon.push_back(beamKaon->getBeta()); //beta

            std::vector<std::vector<std::shared_ptr<DetectedEvent>>> dEPion = detector->getDetectedEvents(beamPion, eventCount);
            std::vector<std::vector<double>> pixelHitsPion = detector->getPixelHitsPerEvents(dEPion, false);
            std::vector<double> meanPion = detector->getMeanPerPixel(pixelHitsPion);

            std::vector<std::vector<std::shared_ptr<DetectedEvent>>> dEKaon = detector->getDetectedEvents(beamKaon, eventCount);
            std::vector<std::vector<double>> pixelHitsKaon = detector->getPixelHitsPerEvents(dEKaon, false);
            std::vector<double> meanKaon = detector->getMeanPerPixel(pixelHitsKaon);

            std::string variableNames = graphingAndAnalysisCherenkovPID::getFileName(betaBool, beamPion); //Should return the same for beamPion or beamKaon

            std::vector<bool> variableSwitch = std::vector<bool>();

            //Arrange Means and PixelHits to be the same as the interpolation
            //
            meanPion = graphingAndAnalysisCherenkovPID::organizeMeans(meanPion, pmt);
            meanKaon = graphingAndAnalysisCherenkovPID::organizeMeans(meanKaon, pmt);

            pixelHitsPion = graphingAndAnalysisCherenkovPID::quickOrg(pixelHitsPion);
            pixelHitsKaon = graphingAndAnalysisCherenkovPID::quickOrg(pixelHitsKaon);

            //double generatedSep = graphingAndAnalysisCherenkovPID::getSeperationValueFromPixelHitsAndMean(detector, pixelHitsPion, pixelHitsKaon, meanPion, meanKaon, variableNames);
            //double generatedSepNuissance = graphingAndAnalysisCherenkovPID::getSeperationValueFromPixelHitsAndMean(detector, pixelHitsPion, pixelHitsKaon, meanPion, meanKaon, variableNames, true); //turn on nuissance parameter

            //Interpolat Mean
            std::vector<double> meanInterpPion = linearInter->interpolationSwitch(variablesPion, interpolateMeanLinearly, variableSwitch,
                                                                                  detector->getPMT()->getWidthPixelCount(), detector->getTotalDistance(),
                                                                                  detector->getMaxWidth(), detector->getMaxLength(),
                                                                                  detector->getPMT()->getWidthPixel(), detector->getPMT()->getLengthPixel());

            std::vector<double> meanInterpKaon = linearInter->interpolationSwitch(variablesKaon, interpolateMeanLinearly, variableSwitch,
                                                                                  detector->getPMT()->getWidthPixelCount(), detector->getTotalDistance(),
                                                                                  detector->getMaxWidth(), detector->getMaxLength(),
                                                                                  detector->getPMT()->getWidthPixel(), detector->getPMT()->getLengthPixel());

            std::vector<double> meanPionDelta;
            std::vector<double> meanKaonDelta;
            for (int pixel = 0; pixel < meanPion.size(); pixel++) {
              meanPionDelta.push_back(meanPion[pixel] - meanInterpPion[pixel]);
              meanKaonDelta.push_back(meanKaon[pixel] - meanInterpKaon[pixel]);
            }
            std::cout << getSTDEV(meanPionDelta) << std::endl;
            std::cout << getSTDEV(meanKaonDelta) << std::endl;
            //double interpolatedSep = graphingAndAnalysisCherenkovPID::getSeperationValueFromPixelHitsAndMean(detector, pixelHitsPion, pixelHitsKaon, meanInterpPion, meanInterpKaon, variableNames);
            //double interpolatedSepNuissance = graphingAndAnalysisCherenkovPID::getSeperationValueFromPixelHitsAndMean(detector, pixelHitsPion, pixelHitsKaon, meanInterpPion, meanInterpKaon, variableNames, true); //turn on nuissance parameter

            //std::cout << "x: " << x << "    ";
            //std::cout << "y: " << y << "    ";
            //std::cout << "pionP: " << beamPion->getMomentum() << "     ";
            //std::cout << "pionB: " << beamPion->getBeta() << "     ";
            //std::cout << "kaonP: " << beamKaon->getMomentum() << "     ";
            //std::cout << "kaonB: " << beamKaon->getBeta() << "     ";
            //std::cout << "theta: " << t << "    ";
            //std::cout << "phi: " << p << "    ";
            //std::cout << "gS: " << generatedSep << "     ";
            //std::cout << "gSN: " << generatedSepNuissance << "     ";
            //std::cout << "iS: " << interpolatedSep << "     ";
            //std::cout << "iSN: " << interpolatedSepNuissance << "     ";
            //std::cout << std::endl;

            //graphingAndAnalysisCherenkovPID::saveSeperationParameter("seperation_Momentum", betaBool, beamPion,
            //generatedSep, generatedSepNuissance,
            //interpolatedSep, interpolatedSepNuissance);
          }
        }
      }
    }
  }
}
