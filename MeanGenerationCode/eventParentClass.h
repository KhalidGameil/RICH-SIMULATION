#ifndef PARENT_EVENT_CLASS_INCLUDE
#define PARENT_EVENT_CLASS_INCLUDE


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <numeric>
#include <memory>
#include "beamClass.h"

class ParentEvent{
   protected:
    int N;
    double cherenkovAngle;
    std::vector<double> interactionPosition;
    std::vector<double> phi;
    std::vector<double> energy;
    std::shared_ptr<Beam> beam;
    std::vector<double> x;
    std::vector<double> y;
    void calcXY();
  public:
    //ParentEvent(){};
    ParentEvent(int ,double ,
                             std::vector<double> ,
                             std::vector<double> ,
                             std::vector<double> ,
                             std::shared_ptr<Beam>);
    int getN();
    double getCherenkovAngle();
    std::vector<double> getInteractionPosition();
    std::vector<double> getPhi();
    std::vector<double> getEnergy();
    std::vector<double> getWavelength();
    std::vector<double> getRadius();
    std::vector<double> getX();
    std::vector<double> getY();
    //std::vector<double> getZ();
};

#endif


