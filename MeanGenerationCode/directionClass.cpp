// in myclass.cpp
#include "directionClass.h"

Direction::Direction(double theta, double phi, double errTheta, double errPhi) {
  this->theta = theta;
  this->phi = phi;
  this->errPhi = errPhi;
  this->errTheta = errTheta;
  savePhi = phi;
  saveTheta = theta;
}

double Direction::getPhi() {
  return phi;
}

double Direction::getTheta() {
  return theta;
}

void Direction::setRandom() {
  std::shared_ptr<TRandom3> randomGenerate(std::make_shared<TRandom3>());
  randomGenerate->SetSeed(1);
  phi = savePhi + randomGenerate->Gaus(0.0, errPhi);
  theta = saveTheta + randomGenerate->Gaus(0.0, errTheta);
}
