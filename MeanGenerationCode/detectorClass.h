#ifndef DETECTOR_CLASS_INCLUDE
#define DETECTOR_CLASS_INCLUDE

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

#include "aerogelClass.h"
#include "beamClass.h"
#include "interpolateLikelihoodClass.h"
#include "pmtClass.h"
#include "randomEventClass.h"

#include "Math/BrentRootFinder.h"
#include "Math/RootFinderAlgorithms.h"
#include "Math/WrappedTF1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TMinuit.h"
#include "TSpline.h"
class Detector {
private:
  std::vector<double> variables;
  std::shared_ptr<TSpline3> wSpline;
  std::shared_ptr<InterpolateLikelihood> inter;
  int interType = 0;
  void setVariables(std::vector<double>);
  
  bool w_Switch = false;

  std::vector<double> energyRangePMT;

  int plotCount;
  std::vector<double> pixelHits;
  double pHSingleEvent=0;
  std::vector<double> mean;
  double mSingleEvent=0;
  std::vector<std::shared_ptr<Aerogel>> aerogels;
  std::vector<std::vector<std::shared_ptr<RandomEvent>>> randomEvents;
  std::vector<std::vector<std::shared_ptr<DetectedEvent>>> detectedEvents;
  std::shared_ptr<PMT> pmt;
  double aeroDetDistance; // Distance from Final Aerogel to face of the
                          // multianodepmt

  void getBinEdges(std::vector< double> *, std::vector< double> *);

  // Setters
  void setCentralPositions();

  //TEST
  void testWrite(std::vector<double> ,std::vector<double>);
  // Event Generators)
  std::vector<std::shared_ptr<RandomEvent>> getEventFromBeam(std::shared_ptr<Beam>);
  std::vector<std::vector<std::shared_ptr<RandomEvent>>> getMultipleEventsFromBeam(std::shared_ptr<Beam>, int);

  // Calculate Private
  double calcLogOfProbabilityOfHit(double);
  double calcLikelihoodForEvent(double, double);
  double calcStepSize(double, double, double);

  //Minuit Setter and Getters
  TMinuit * m;
  void setMinuit(TMinuit *g){m=g;};


public:
  Detector(std::vector<std::shared_ptr<Aerogel>>, std::shared_ptr<PMT>, double);
  TMinuit * getMinuit(){return m;};

  //SplineGetter
  std::shared_ptr<TSpline3> getWSpline() { return wSpline; }
  double EvaluateSpline(TSpline3 *spline, int nSteps, double max, double min, double stepSize, double pos);

  //PixelHit Setter and Getter
  void setPixelHit(std::vector<double>);
  std::vector<double> getPixelHit();
  void setPixelHitSingleEvent(double);
  double getPixelHitSingleEvent();

  //Interpolator Setter and Getter;
  void setInterpolator(std::shared_ptr<InterpolateLikelihood> i,int t) {inter = i;interType = t;}
  std::shared_ptr<InterpolateLikelihood> getInterpolator() { return inter; }
  int getInterpolatorType() {return interType;}

  //W parameter Switch get
  void setWSwitch(bool w_ON){w_Switch=w_ON;}
  bool getWSwitch() {return w_Switch;}

  //Mean Setter and Getter;
  void setMean(std::vector<double>);
  std::vector<double> getMean();
  void setMeanSingleEvent(double);
  double getMeanSingleEvent();

  // Getters
  double getMaxWidth();
  double getMaxLength();
  std::vector<double> getVariables();

  // Simple (some calculation involved) Getters
  std::vector<std::shared_ptr<Aerogel>> getAerogelArray();
  std::shared_ptr<PMT> getPMT();
  double getAeroDetDistance();
  double getTotalDistance();
  std::vector<double> getGeneratedEnergyRange(int, std::vector<double>);

  // Add Object
  void addAerogel(std::shared_ptr<Aerogel>);

  // Pixel Hit Generators
  std::vector<std::vector<std::shared_ptr<RandomEvent>>> getRandomEvents();
  std::vector<std::vector<std::shared_ptr<DetectedEvent>>> getDetectedEvents(std::shared_ptr<Beam>, int);
  std::vector<double> getEventsForPixel(std::vector<std::vector<double>>, int);
  std::vector<std::vector<double>> getPixelHitsPerEvents(std::vector<std::vector<std::shared_ptr<DetectedEvent>>>, bool);

  // Pixel Hit Analysis
  std::vector<double> getSumOfHitsPerPixel(std::vector<std::vector<double>>);
  std::vector<std::shared_ptr<TH1D>> getArrayOfHistograms(std::vector<std::vector<double>>);
  std::vector<std::shared_ptr<TGraph>> getArrayOfExpectedPoissonHistogramsPerPixel(std::vector<std::vector<double>>);

  // Statistical Analysis Functions
  std::vector<double> getMeanPerPixel(std::vector<std::vector<double>>);
  std::vector<double> getStandardDeviationForPixel(std::vector<std::vector<double>>);
  std::vector<double> calcMeanRelativeUncertainty(std::vector<double>, int);
  std::vector<double> calcWeight(std::vector<double>, int);
  double calcWeightSum(std::vector<double>, int);
  double calcWeightSum(std::vector<double>);

  // Likelihood function;
  std::vector<double> getLogLikelihood(std::vector<std::vector<double>>, std::vector<double>, std::vector<double>);
  std::vector<double> getWeightedLogLikelihood(std::vector<std::vector<double>>, std::vector<double>, std::vector<double>);
  std::vector<std::vector<double>> getLogLikelihoodPerPixel(std::vector<std::vector<double>>, std::vector<double>, std::vector<double>);
  //std::vector<std::vector<double>> getWeightedLogLikelihoodPerPixel(std::vector<std::vector<double>>, std::vector<double>, std::vector<double>);
  //std::vector<std::vector<double>> getLogLikelihoodNuissancePerPixel(std::vector<std::vector<double>>, std::vector<double>, std::vector<double>);

  //Nuisance Parameter Source Program
  std::vector<double> getLogLikelihoodWithNuisanceParameters(std::vector<std::vector<double>>,
                                                             std::vector<double>,
                                                             std::vector<double>,
                                                             std::vector<double>,
                                                             std::vector<double>,
                                                             std::shared_ptr<InterpolateLikelihood>,
							     int, bool,
                                                             std::string);
  std::vector<double> getLogLikelihoodWithNuisanceParameters(std::vector<std::vector<double>>,
                                                             std::vector<double>,
                                                             std::vector<double>,
                                                             std::shared_ptr<Beam>,
                                                             std::shared_ptr<Beam>,
                                                             std::shared_ptr<InterpolateLikelihood>,
							     int, bool,
                                                             std::string);
  double calcLikelihoodForParametersNuisanceForParticle(TMinuit *, std::vector<double>, double);

  //Minuit Nuisance
  std::vector<double> getLogLikelihoodWithNuisanceMinuit(std::vector<std::vector<double>>,
                                                         std::vector<double>,
                                                         std::vector<double>);

  std::vector<double> getLogLikelihoodForParticleWithNuisanceMinuit(std::vector<std::vector<double>>,
                                                                    std::vector<double>,
                                                                    double);
  double calcLikelihoodForParticleNuisanceMinuit(std::vector<double>,
                                                 std::vector<double>,
                                                 double);

  double calcLikelihoodForPixelNuisance(std::vector<double>,
                                        std::vector<double>,
                                        double,
                                        double);

  double calcLikelihoodForEventNuisance(double, double, double, double);

  double calcLikelihoodWithParameterNuisanceMinuit(TMinuit *, double, std::vector<double>);

  //Read Parameters
  void quickReadWParameter();
  //Nuisance Parameter read and write
  void quickSaveWParameter(double, double, double);
  void quickHistogramWParameter(double);

  //void FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag);

  //delete this code
  std::vector<double> getLogLikelihoodStupid(std::vector<std::vector<double>>, std::vector<double>, std::vector<double>);

  //Nuisance
  //std::vector<double> getLogLikelihoodWithNuisance(std::vector<std::vector<double>>, std::vector<double>, std::vector<double>, bool);
  //double getNuisanceParameter(double, double, double, double);
  //double getNuisanceParameterMINUIT(double, double, double, double);

  // Log Likelihood Test Code [Scrap?]
  std::vector<double> getLogLikelihoodValuesForSingleParticle(std::vector<std::vector<double>>, std::vector<double>);
  double getLogLikelihoodValuesForSingleEventSingleMean(std::vector<double>, double);
  //void testSeperationLikelihoodOnVariable(int, std::shared_ptr<InterpolateLikelihood>, std::string,
  //std::vector<double>, std::vector<double>,
  //std::vector<std::vector<double>>,
  //std::vector<std::vector<double>>, std::vector<double> &,
  //std::vector<double> &, std::vector<double> &,
  //std::vector<double> &);

  void getListOfVariables(int, int, double, double, std::vector<bool>, std::vector<double>,
                          std::vector<double>, std::vector<double> &, std::vector<double> &,
                          std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);

  //TEST Code
  TMinuit *getWMinuit();
  void graphHistogramMinuit(int, TMinuit *);
  double getMinuitVariableAndWParameterMinimization(TMinuit *,
                                                    std::vector<double>,
                                                    std::vector<double>,
                                                    double,
                                                    std::vector<double>,
                                                    std::string);
  double getMinuitVariableAndWParameterMinimization(TMinuit *,
                                                    std::vector<double>,
                                                    std::vector<double>,
                                                    double,
                                                    std::shared_ptr<Beam>,
                                                    std::string);
  double getMinuitVariableAndWParameterMinimizationReturnParameter(std::shared_ptr<TMinuit>,
                                                                   std::vector<double>,
                                                                   std::vector<double>,
                                                                   double,
                                                                   std::vector<double> &,
                                                                   std::vector<double> &,
                                                                   std::vector<double> &,
                                                                   std::shared_ptr<Beam>,
                                                                   std::string);

  std::vector<double> getLogLikelihoodForVariableAndWparameter(std::vector<std::vector<double>>,
                                                               std::vector<double>,
                                                               std::vector<double>,
                                                               std::vector<double>,
                                                               std::vector<double>,
                                                               std::shared_ptr<InterpolateLikelihood>,
							       int,
							       bool,
                                                               std::string);

  std::vector<double> getLogLikelihoodForVariableAndWparameter(std::vector<std::vector<double>>,
                                                               std::vector<double>,
                                                               std::vector<double>,
                                                               std::shared_ptr<Beam>,
                                                               std::shared_ptr<Beam>,
                                                               std::shared_ptr<InterpolateLikelihood>,
							       int,
							       bool,
                                                               std::string);

  void getLogLikelihoodForVariableAndWparameterSingleParticle(std::vector<std::vector<double>>,
                                                              std::vector<double>,
                                                              std::vector<std::vector<double>> &,
                                                              std::vector<std::vector<double>> &,
                                                              std::vector<std::vector<double>> &,
                                                              std::vector<double> &,
                                                              std::shared_ptr<Beam>,
                                                              std::shared_ptr<InterpolateLikelihood>,
							      int,bool,
                                                              std::string);
};

#endif
