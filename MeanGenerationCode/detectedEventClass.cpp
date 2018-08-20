#include "detectedEventClass.h"

DetectedEvent::DetectedEvent(int N,double cherenkovAngle,
                             std::vector<double> interactionPosition,
                             std::vector<double> theta,
                             std::vector<double> energy,
                             std::shared_ptr<Beam> beam):
  ParentEvent(N,cherenkovAngle,interactionPosition,theta,energy,beam)
{ }
