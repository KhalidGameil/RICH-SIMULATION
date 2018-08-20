#ifndef PMT_CLASS_INCLUDE
#define PMT_CLASS_INCLUDE

#include "directionClass.h"
#include "positionClass.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

#include "TGraphErrors.h"
#include "TMath.h"
#include "TSpline.h"
struct crossTalk {
  double neighbor;
  double diagonal;
};

class PMT {
private:
  double width=0;
  double length=0;
  double pixelCount=0;
  crossTalk mask;

  double windowThickness=0;
  double deadSpace=0;
  double darkNoise=0;

  double timeResolution=0;

  std::shared_ptr<Position> position;
  std::vector<double> wavelength;
  std::vector<double> quantumEfficiency;

  void setEnergyAndQuantumEfficiency();
  void setCrossTalkValues();
  bool checkCrossTalk(int);

public:
  PMT(double, double, int);
  double getWidth();
  double getLength();

  double getWidthPixel();
  double getLengthPixel();

  int getPixelCount();

  void setWindowThickness(double);
  double getWindowThickness();

  void setDeadSpace(double);
  double getDeadSpace();

  void setDarkNoise(double);
  double getDarkNoise();

  void setTimeResolution(double);
  double getTimeResolution();

  void setPosition(std::shared_ptr<Position>);
  std::shared_ptr<Position> getPosition();

  int getLengthPixelCount();
  int getWidthPixelCount();

  std::vector<double> getEnergy();
  std::vector<double> getWavelength();
  std::vector<double> getQuantumEfficiency();

  std::shared_ptr<TSpline3> getQuantumEfficiencySpline();
  std::vector<double> applyCrossTalkMask(std::vector<double>);
  std::vector<double> applyCrossTalkTH2DMask(std::vector<double>);
};

#endif
