#ifndef BEAM_CLASS_INCLUDE
#define BEAM_CLASS_INCLUDE

#include "directionClass.h"
#include "positionClass.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

//Root Matrix
#include "Riostream.h"
#include "TCanvas.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "TH1D.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TRandom3.h"
#include "TVectorD.h"
class Beam {
private:
  std::shared_ptr<Position> position;
  std::shared_ptr<Direction> direction;
  double momentum;
  double particle;
  double beta;
  double errBetaMomentum;
  double savedBeta;
  double savedMomentum;

  //TEST CODE
  std::vector<double> betaList;
  std::vector<double> mList;
  std::vector<double> xList;
  std::vector<double> yList;
  std::vector<double> tList;
  std::vector<double> pList;

public:
  Beam(double, double, std::shared_ptr<Direction>, std::shared_ptr<Position>, double errBetaMomentum = 0);
  Beam(double, std::shared_ptr<Direction>, std::shared_ptr<Position>, double errBetaMomentum = 0);
  std::shared_ptr<Position> getPosition();
  std::shared_ptr<Direction> getDirection();
  double getMomentum();
  double getParticleMass();
  void setParticleMass(double);
  double getBeta();
  void getProjectedXYZ(double *, double *, double, double, double);
  double getLengthOfVector(double);
  void setRandom();
  void plotAll(std::string, std::string);
};

#endif
