#include "eventParentClass.h"

ParentEvent::ParentEvent(int N,double cherenkovAngle,
                             std::vector<double> interactionPosition,
                             std::vector<double> phi,
                             std::vector<double> energy,
                             std::shared_ptr<Beam> beam)
{
  this->N = N;
  this->cherenkovAngle = cherenkovAngle;
  this->interactionPosition  = interactionPosition;
  this->phi = phi;
  this->energy = energy;
  //this->beamX = beamX;
  //this->beamY = beamY;
  //this->noise = noise;
  this->beam = beam;
  calcXY();
}

int ParentEvent::getN()
{
  return N;
}
double ParentEvent::getCherenkovAngle()
{
  return cherenkovAngle;
}
std::vector<double> ParentEvent::getInteractionPosition()
{
  return interactionPosition;
}

std::vector<double> ParentEvent::getPhi()
{
  return phi;
}
std::vector<double> ParentEvent::getEnergy()
{
  return energy;
}
std::vector<double> ParentEvent::getWavelength()
{
  double h = 4.135E-15; //eV*s
  double c = 299792458; //m/s
  std::vector<double> wvl;
  wvl.reserve(energy.size());
  for (int i = 0; i < energy.size();i++)
  {
    wvl.push_back( h*c/(energy.at(i)));
  }
  return wvl;
}
std::vector<double> ParentEvent::getRadius()
{
  std::vector<double> radius;
  radius.reserve(interactionPosition.size());
  for ( int i = 0; i < interactionPosition.size(); i++)
  {
    radius.push_back(interactionPosition.at(i)*tan(cherenkovAngle));
  }
  return radius;
}

std::vector<double> ParentEvent::getX()
{
  return x;
}

std::vector<double> ParentEvent::getY()
{
  return y;
}

//std::vector<double> ParentEvent::getZ()
//{

//}

void ParentEvent::calcXY()
{
  std::vector<double> r = this->getRadius();
  //std::cout << "TEST XY 1" << std::endl;
  x.reserve(r.size());
  y.reserve(r.size());
  std::vector<double> z = interactionPosition;
  
  for (int i =0; i < r.size(); i ++)
  {
    //std::cout << z.at(i) << "   " << r.at(i) << std::endl;
    double temp_x = r.at(i)*cos(phi.at(i));
    double temp_y = r.at(i)*sin(phi.at(i));
    //std::cout << "TEST XY 2 " << i  << std::endl;
    beam->getProjectedXYZ(&temp_x,&temp_y,z.at(i),cherenkovAngle,phi.at(i));
    //std::cout << "TEST XY 3 " << i  << std::endl;
    //std::cout << temp_x << "    " << temp_y << std::endl;
    x.push_back(temp_x);
    y.push_back(temp_y);
  }
}
