#ifndef POSITION_CLASS_INCLUDE
#define POSITION_CLASS_INCLUDE

#include "TRandom3.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
#include <vector>
#include <memory>

class Position {
private:
  double x=0;
  double y=0;
  double z=0;
  double errX=0;
  double errY=0;
  double errZ=0;
  double savedX=0;
  double savedY=0;
  double savedZ=0;

public:
  Position(double, double, double, double errX = 0, double errY = 0, double errZ = 0);
  double getX();
  double getY();
  double getZ();
  void setRandom();
};

#endif
