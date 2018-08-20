#include "randomEventClass.h"

RandomEvent::RandomEvent(int N,double cherenkovAngle,
                             std::vector<double> interactionPosition,
                             std::vector<double> theta,
                             std::vector<double> energy,
                             std::shared_ptr<Beam> beam,
                             std::shared_ptr<TSpline3> spl):
  ParentEvent(N,cherenkovAngle,interactionPosition,theta,energy,beam)
{ this->efficiencySpline = spl;
  randomGenerate=std::make_shared<TRandom3>();
  randomGenerate->SetSeed(0);
}

//int RandomEvent::getN()
//{
  //return N;
//}
//double RandomEvent::getCherenkovAngle()
//{
  //return cherenkovAngle;
//}
//std::vector<double> RandomEvent::getInteractionPosition()
//{
  //return interactionPosition;
//}

//std::vector<double> RandomEvent::getTheta()
//{
  //return theta;
//}
//std::vector<double> RandomEvent::getEnergy()
//{
  //return energy;
//}
//std::vector<double> RandomEvent::getWavelength()
//{
  //double h = 4.135E-15; //eV*s
  //double c = 299792458; //m/s
  //std::vector<double> wvl;
  //for (int i = 0; i < energy.size();i++)
  //{
    //wvl.push_back( h*c/(energy.at(i)));
  //}
  //return wvl;
//}
//std::vector<double> RandomEvent::getRadius()
//{
  //std::vector<double> radius;
  //for ( int i = 0; i < interactionPosition.size(); i++)
  //{
    //radius.push_back(interactionPosition.at(i)*tan(cherenkovAngle));
  //}
  //return radius;
//}

double RandomEvent::evaluateSpline(double e)
{

  double h = 4.135E-15; //eV*s
  double c = 299792458; //m/s
  double wvl = 1E9*h*c/(e);
  //std::cout << wvl << "   ";
  return efficiencySpline->Eval(wvl);
}

std::vector<double> RandomEvent::getQuantumEfficiency()
{
  std::vector<double> qE;
  for (int i = 0; i < energy.size(); i++)
  {
    qE.push_back(evaluateSpline(energy.at(i))/100.0);
    //std::cout << qE.back() << std::endl;
  }
  return qE;
}

std::shared_ptr<DetectedEvent> RandomEvent::getDetectedHits()
{
  std::vector<double> qE = this->getQuantumEfficiency();
  std::vector<double> detectedIP;//interaction position
  std::vector<double> detectedTheta;
  std::vector<double> detectedEnergy;
  int energyCount = energy.size();
  detectedIP.reserve(energyCount);
  detectedTheta.reserve(energyCount);
  detectedEnergy.reserve(energyCount);

  for(int i =0; i < energyCount; i++)
  {
    double rG = randomGenerate->Uniform(0.0,1.0);
    
    //if (i%5 == 0){
    //std::cout << rG << "   " << qE.at(i) <<  "   " <<  (rG<=qE.at(i)) <<std::endl;}

    if (randomGenerate->Uniform(0.0,1.0)<= qE.at(i))
    {
      detectedIP.push_back(interactionPosition.at(i));
      detectedTheta.push_back(phi.at(i));
      detectedEnergy.push_back(energy.at(i));
    }
  }
  //std::cout << "N:" << detectedIP.size() << std::endl;//<< "   ";
  //std::cout << "time: " << t << std::endl;
  std::shared_ptr<DetectedEvent> dE(std::make_shared<DetectedEvent>(detectedIP.size(),
                                                                    cherenkovAngle,
                                                                    detectedIP,
                                                                    detectedTheta,
                                                                    detectedEnergy,
                                                                    beam));
  return dE;
}

