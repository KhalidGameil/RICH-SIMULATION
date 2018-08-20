#ifndef DETECTED_EVENT_CLASS_INCLUDE
#define DETECTED_EVENT_CLASS_INCLUDE


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <numeric>
#include <memory>

#include "eventParentClass.h"

class DetectedEvent: public ParentEvent{
  public:
    using ParentEvent::ParentEvent;
    DetectedEvent(int N,double cherenkovAngle,
        std::vector<double> interactionPosition,
        std::vector<double> theta,
        std::vector<double> energy,
        std::shared_ptr<Beam> beam);
};

#endif




