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
#include "TSpline.h"
#include "TString.h"
#include "TStyle.h"
#include "TVirtualPad.h"

//Classes
#include "aerogelClass.h"
#include "detectedEventClass.h"
#include "detectorClass.h"
#include "directionClass.h"
#include "interpolatorClass.h"
#include "pmtClass.h"
#include "positionClass.h"
#include "randomEventClass.h"

//NameSpaces
#include "graphFunctionsNamespace.h"

double calcStepSize(double start, double end, double step) {
  return (end - start) / (step - 1);
}
int main(int argc, char **argv) {
   graphingAndAnalysisCherenkovPID::printArgumentHelperList();

	if (argc !=11){
                std::cout << "Number of Arguments:" << argc << std::endl;
		std::cout << "ARGUMENT ERROR" << std::endl;
		return 0;
	}
  
  graphingAndAnalysisCherenkovPID::printArgList(argc,argv);
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

  //File Name Append
  std::string fileNameAppend = argv[1];

  //Event
  int eventCount = std::atoi(argv[2]);
  int interpolatorSwitch = std::atoi(argv[3]);
  bool wSwitch = std::atoi(argv[4]);
  std::string filename = "MeanInterpolationDataBank_database_2018Jan26.txt";
  std::shared_ptr<InterpolateLikelihood> linearInter(std::make_shared<InterpolateLikelihood>(filename));
  std::shared_ptr<Detector> detector(std::make_shared<Detector>(aerogels, pmt, .09)); //.06)); //CHANGE=.09

  //Parameter Arguments
  double x = std::stod(argv[5]);                //m //#CHANGE=.06
  double y = std::stod(argv[6]);                //m //#CHANGE=.06
  double z = 0;
  double t = std::stod(argv[7]);
  double p = std::stod(argv[8]);
  bool betaBool = bool(std::atoi(argv[9]));
  double i = std::stod(argv[10]);

    std::cout << "x: " << x << std::endl;
    std::cout << "y: " << y << std::endl;
    std::cout << "theta: " << t << std::endl;
    std::cout << "phi: " << p << std::endl;
    std::cout << "beta: " << i << std::endl;

    std::clock_t begin = std::clock();

    std::shared_ptr<Position> pos(std::make_shared<Position>(x, y, z));
    std::shared_ptr<Direction> dir(std::make_shared<Direction>(t, p));
    std::shared_ptr<Beam> beamPion;

    if (betaBool == true) {
      beamPion = std::make_shared<Beam>(i, dir, pos);
      beamPion->setParticleMass(pionMass);
    } else {
      beamPion = std::make_shared<Beam>(i, pionMass, dir, pos);
    }

    std::vector<double> variablesPion;
    variablesPion.reserve(5);
    variablesPion.push_back(x);                   //x
    variablesPion.push_back(y);                   //y
    variablesPion.push_back(t);                   //theta
    variablesPion.push_back(p);                   //phi
    variablesPion.push_back(beamPion->getBeta()); //beta

    std::vector<std::vector<std::shared_ptr<DetectedEvent>>> dEPion = detector->getDetectedEvents(beamPion, eventCount);
    std::vector<std::vector<double>> pixelHitsPion = detector->getPixelHitsPerEvents(dEPion, false);
    std::vector<double> meanPion = detector->getMeanPerPixel(pixelHitsPion);

    std::string variableNames = graphingAndAnalysisCherenkovPID::getFileName(betaBool, beamPion); //Should return the same for beamPion or beamKaon

    std::cout << "PIXEL COUNT:" << meanPion.size() << "   ";
    std::vector<double> meanInterpPion;
    std::vector<double> meanInterpKaon;
    std::vector<bool> variableSwitch = std::vector<bool>();
    variableSwitch.push_back(false);//x
    variableSwitch.push_back(false);//y
    variableSwitch.push_back(false);//theta
    variableSwitch.push_back(false);//phi
    variableSwitch.push_back(true);//beta

    meanInterpPion = linearInter->interpolationSwitch(variablesPion, interpolatorSwitch, variableSwitch,
                                                      detector->getPMT()->getWidthPixelCount(), detector->getTotalDistance(),
                                                      detector->getMaxWidth(), detector->getMaxLength(),
                                                      detector->getPMT()->getWidthPixel(), detector->getPMT()->getLengthPixel());
    
    std::cout << "x: " << x << "    ";
    std::cout << "y: " << y << "    ";
    std::cout << "pionP: " << beamPion->getMomentum() << "     ";
    std::cout << "pionB: " << beamPion->getBeta() << "     ";
    std::cout << "theta: " << t << "    ";
    std::cout << "phi: " << p << "    ";
    std::cout << std::endl;

    std::vector<std::vector<double>> v;
    std::vector<std::vector<double>> v_err;
    std::vector<std::vector<double>> v_war;
    std::vector<double> ll;
    detector->getLogLikelihoodForVariableAndWparameterSingleParticle(pixelHitsPion, meanInterpPion, 
									v, v_err, v_war, 
									ll, beamPion, linearInter,
									interpolatorSwitch, wSwitch,
									"PION");


    std::string fileName = variableNames + "wParaMeter_test_Pion"+fileNameAppend;
    graphingAndAnalysisCherenkovPID::saveParameterAndLikelihood(fileName,
                                                                betaBool, beamPion,
                                                                v,
                                                                v_err,
                                                                v_war,
                                                                ll);

    std::clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "TIME THE CODE TOOK:" << elapsed_secs << std::endl;
    return 0;
}
