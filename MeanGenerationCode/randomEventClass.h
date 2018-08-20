#ifndef RANDOM_EVENT_CLASS_INCLUDE
#define RANDOM_EVENT_CLASS_INCLUDE


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <numeric>
#include<ctime>
#include <sys/time.h>
#include <memory>


#include "TSpline.h"
#include "TRandom3.h"


#include "eventParentClass.h"
#include "detectedEventClass.h"

class RandomEvent: public ParentEvent{
  private:
    //int N;
    //double cherenkovAngle;
    //std::vector<double> interactionPosition;
    //std::vector<double> theta;
    //std::vector<double> energy;
    std::shared_ptr<TSpline3> efficiencySpline;
    std::shared_ptr<TRandom3> randomGenerate;
    double evaluateSpline(double);
    
  public:
    using ParentEvent::ParentEvent;
    RandomEvent(int N,double cherenkovAngle,
                             std::vector<double> interactionPosition,
                             std::vector<double> theta,
                             std::vector<double> energy,
                             std::shared_ptr<Beam> beam,
                             std::shared_ptr<TSpline3> spl);
    //int getN();
    //double getCherenkovAngle();
    //std::vector<double> getInteractionPosition();
    //std::vector<double> getTheta();
    //std::vector<double> getEnergy();
    //std::vector<double> getWavelength();
    //std::vector<double> getRadius();
    std::vector<double> getQuantumEfficiency();
    std::shared_ptr<DetectedEvent> getDetectedHits();
};

#endif



