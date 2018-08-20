#ifndef DIRECTION_CLASS_INCLUDE
#define DIRECTION_CLASS_INCLUDE

#include "TRandom3.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
#include <vector>
#include <memory>

class Direction {
private:
  double theta=0;
  double phi=0;
  double errPhi=0;
  double errTheta=0;
  double savePhi=0;
  double saveTheta=0;

public:
  Direction(double, double, double errTheta = 0, double errPhi = 0);
  double getTheta();
  double getPhi();
  void setRandom();
};

#endif
