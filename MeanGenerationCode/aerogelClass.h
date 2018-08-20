#ifndef AEROGEL_CLASS_INCLUDE
#define AEROGEL_CLASS_INCLUDE

#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
#include <vector>
#include <memory>

//Root
#include "TMath.h"
#include "TRandom3.h"
#include "TSpline.h"

#include "beamClass.h"
#include "positionClass.h"
#include "randomEventClass.h"
class Aerogel {
private:
  double refractiveIndex=0.0;
  double length=0.0;
  double width=0.0;
  double thickness=0.0;
  std::shared_ptr<TRandom3> randomGenerate;
  std::shared_ptr<Position> position;

  std::vector<double> calcEnergyPdf(double, std::vector<double>);
  std::vector<double> calcEnergyCdf(double, std::vector<double>);
  std::vector<double> findCdf(std::vector<double>, std::vector<double>, std::vector<double>);
  std::vector<double> getRandomInteractionPosition(int, std::shared_ptr<Beam>, double);
  std::vector<double> getRandomTheta(int);
  std::vector<double> getRandomEnergy(int, double, double);
  double getCherenkovAngle(double);
  std::vector<double> getRandomInvCDF(int, double, std::vector<double>, std::vector<double>);

public:
  Aerogel(double, double, double, double);
  double getRefractiveIndex();
  double getLength();
  double getWidth();
  double getThickness();

  double calcN(std::shared_ptr<Beam>, std::vector<double>);
  void setPosition(std::shared_ptr<Position>);
  std::shared_ptr<Position> getPosition();

  std::vector<double> calcSellmeierRefractiveIndex(std::vector<double>);
  std::vector<double> calcRefractiveIndexEnergy(std::vector<double>); //should
  std::shared_ptr<RandomEvent> getEvent(double, std::shared_ptr<Beam>, std::vector<double>, std::shared_ptr<TSpline3>);
};

#endif
