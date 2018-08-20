#ifndef INTERLIKE_CLASS_INCLUDE
#define INTERLIKE_CLASS_INCLUDE

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <numeric>
#include <sstream> // std::istringstream
#include <string>
#include <sys/stat.h>
#include <vector>

#include "TCanvas.h"
#include "TCollection.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TSpline.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystemDirectory.h"

//Geometric Toolbox
//#include <Mathematics/GteIntpTricubic3.h>

//Boost
//#include <boost/filesystem.hpp>
//#include <boost/iostreams/device/mapped_file.hpp> // for mmap
#include </usr/local/Cellar/boost/1.64.0_1/include/boost/iostreams/device/mapped_file.hpp>

struct measuredBatch {
  double p;
  double beta;
  double x;
  double y;
  double theta;
  double phi;
  std::vector<std::vector<double>> mean;
  std::vector<std::string> particle;
  std::vector<double> mass;
};

struct NColumnOnlyCmp {
  NColumnOnlyCmp(int columnN) { this->N = columnN; }
  bool operator()(const std::vector<double> &lhs,
                  const std::vector<double> &rhs) const {
    return lhs[N] < rhs[N];
  }
  int N;
};

class InterpolateLikelihood {
private:
  std::string folderLocation;
  std::vector<measuredBatch> measurements;
  void readTextFile(std::string);
  double calcBeta(double, double);
  double calcMomentum(double, double);
  double pionMass = 137.2735;
  double kaonMass = 495.665;
  double protonMass = 0; //not programed yet
  std::vector<double> betaList;
  std::vector<double> xList;
  std::vector<double> yList;
  std::vector<double> thetaList;
  std::vector<double> phiList;
  void readDataBaseFile(std::string);
  void readDataBaseFileFastC(std::string);
  void getVariableBoundariesList(std::vector<double>, int &, int &);
  std::vector<std::shared_ptr<TSpline3>> getSplines(std::vector<double>, std::vector<int>, std::string);
  void getListOfVariableIndex(std::vector<bool>, std::vector<double>, std::vector<int> &, std::vector<double> &);

public:
  ////InterpolateLikelihood(InterpolateLikelihood *);
  InterpolateLikelihood();
  InterpolateLikelihood(std::string);
  int countLines(boost::iostreams::mapped_file mmap, const char *f);
  void readDataBaseFileMMapped(std::string);
  void readDataBaseFileFast(std::string);
  //InterpolateLikelihood(std::string,
  //int, int,
  //int, int,
  //int);

  std::string getFolderLocation();
  void setFolderLocation(std::string);
  void addMeasurement(double, double, double, double, double, std::vector<double>);

  void readFolder();
  void checkMeasurements();
  void sortVariables();
  void checkVector(std::vector<int>);
  void checkVector(std::vector<double>);
  double calcDelta(std::vector<double>);
  std::vector<measuredBatch> getMeasurements() { return measurements; };

  std::vector<std::vector<std::shared_ptr<TGraphErrors>>> getMeanVariableGraph(std::string);
  std::vector<std::vector<std::shared_ptr<TGraph2D>>> getMean2VariableGraph(std::string);
  //std::vector<std::vector<std::shared_ptr<TH2F>>> getMean2VariableHist(std::string);
  std::vector<std::vector<double>> interpolateMean(std::string);

  void SaveImageOfGraphs(std::vector<std::vector<std::shared_ptr<TGraphErrors>>>, std::string);
  void SaveImageOfEachGraph(std::vector<std::vector<std::shared_ptr<TGraphErrors>>>, std::string);
  void SaveImageOf2DGraphs(std::vector<std::vector<std::shared_ptr<TGraph2D>>>, std::string);
  void SaveImageOfEach2DGraphs(std::vector<std::vector<std::shared_ptr<TGraph2D>>>, std::string);
  std::vector<double> getVariablesFromIndex(std::vector<int>);
  std::vector<double> getVariablesFromIndex(int);
  void removeDuplicates(std::vector<std::vector<double>> &);
  std::vector<std::vector<double>> getAllPoints(std::vector<double>, std::vector<double>);
  std::vector<double> interpolateMeanToFirstOrder(std::vector<double>);

  //LINEAR INTERPOLATION FUNCTIONS AND SUB FUNCTIONS
  std::vector<int> getIndex(std::vector<std::vector<double>>);
  std::vector<double> getDifferenceBetween2Vectors(std::vector<double> v1, std::vector<double>);
  int findFurthest(int vI, std::vector<std::vector<double>>);
  double calcDistance(std::vector<double> v, std::vector<double>);
  std::vector<double> returnSwap(std::vector<double>, std::vector<double>, std::vector<int>);
  std::vector<double> interpolateMeanLinearly(std::vector<double>);
  int findMeanPhotonList(std::vector<double>);
  int getUpperBoundaryIndex(std::vector<double>, double);

  //Bools
  bool isVectorPositive(std::vector<int>);
  bool isPositive(int);
  void printVariablLists();
  void readVariableLists();
  void testIndexFinder(std::vector<double>);
  void findVariableGivenMean(int, double,
                             double &, double &,
                             double &, double &,
                             double &,
                             std::vector<bool>);

  //SPline Correrct
  std::vector<std::vector<double>> getSplineValuesAlongAxis(std::vector<std::vector<double>>, std::vector<double>, std::vector<bool>);
  std::vector<std::vector<double>> removeVariable(std::vector<std::vector<double>>, std::vector<bool>);
  std::vector<std::vector<double>> getPointsAlongAxis(std::vector<bool>, std::vector<std::vector<double>>);
  std::vector<double> interpolateMeansLinearlyWithSpline(std::vector<double>, std::vector<bool>);

  //Gaussian Fit:
  double calcGaussian(double, double, double, double);
  std::vector<std::vector<double>> getGaussianValuesAlongAxis(std::vector<std::vector<double>>, std::vector<double>, std::vector<bool>);
  //std::vector<std::vector<double>> getPointsAlongAxis(std::vector<bool>, std::vector<std::vector<double>>);
  std::vector<TFitResultPtr> getGaussian(std::vector<double>, std::vector<int>, std::string);
  std::vector<double> interpolateMeansLinearlyWithGaussianFit(std::vector<double>, std::vector<bool>);

  //TribCubic X Y Theta
  //std::vector<double> interpolationTriCubicXYTheta(std::vector<double>);
  //std::vector<double> getTriCubicInterpolationValue(std::vector<int>, std::vector<double>, std::vector<bool>);

  //Interpolation Switch
  std::vector<double> interpolationSwitch(std::vector<double>, int, std::vector<bool>,
                                          double, double, double, double, double, double);

  double getMinX();
  double getMaxX();

  double getMinY();
  double getMaxY();

  double getMinTheta();
  double getMaxTheta();

  double getMinPhi();
  double getMaxPhi();

  double getMinBeta();
  double getMaxBeta();

  //Interpolation Algorthims
  void calcPixelPositionSteps(double, double, double, double &);
  void calcPixelThetaSteps(double, double, double, double &);
  void calcPixelPhiSteps(double, double &);

  void shiftBeamPosition(double detectorSizex, double detectorSizey,
                         double pixelLengthx, double pixelLengthy,
                         double xSteps, double ySteps,
                         std::vector<double> &variables);
  void shiftBeamTheta(double &);
  void shiftBeamPhi(double, double, double, double, double, std::vector<double> &, double &, double &);

  void shiftPixelMapPosition(double, double, double, std::vector<double> &);
  void shiftPixelMapTheta(double, double, double, double, std::vector<double> &);
  std::vector<double> shiftPixelMapPhiPosition(double, double, double, double, std::vector<double>);

  void transposeXY(std::vector<double> &, double &, double &);
  void flipVSteps(double &, double &, double, double);

  void transposePixelMap(double, double &, double &, std::vector<double> &);
  void flipPixelMapAlongAxis(double, double &, double &, std::vector<double> &, char);
};

#endif
