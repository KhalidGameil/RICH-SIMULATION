// in myclass.cpp
#include "positionClass.h"

Position::Position(double x, double y, double z,
                   double errX, double errY, double errZ) {
  this->x = x;
  this->y = y;
  this->z = z;
  this->errX = errX;
  this->errY = errY;
  this->errZ = errZ;
  savedX = x;
  savedY = y;
  savedZ = z;
}

double Position::getX() {
  return x;
}

double Position::getY() {
  return y;
}

double Position::getZ() {
  return z;
}


//Code SHOULD NOT BE IN USE
void Position::setRandom() {
  std::shared_ptr<TRandom3> randomGenerate(std::make_shared<TRandom3>());
  randomGenerate->SetSeed(0);
  x = savedX + randomGenerate->Gaus(0.0, errX);
  y = savedY + randomGenerate->Gaus(0.0, errY);
  z = savedZ + randomGenerate->Gaus(0.0, errZ);
}
