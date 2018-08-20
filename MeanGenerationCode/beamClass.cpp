// in myclass.cpp
#include "beamClass.h"

Beam::Beam(double momentum, double particleMass, std::shared_ptr<Direction> direction, std::shared_ptr<Position> position, double errBetaMomentum) {
  this->direction = direction;
  this->position = position;
  this->momentum = momentum;
  savedMomentum = momentum;
  this->particle = particleMass;
  this->beta = 0;
  savedBeta = 0;
  this->errBetaMomentum = errBetaMomentum;

  //TEST CODE
  betaList.reserve(10000);
  mList.reserve(10000);
  xList.reserve(10000);
  yList.reserve(10000);
  tList.reserve(10000);
  pList.reserve(10000);
}

Beam::Beam(double beta, std::shared_ptr<Direction> direction, std::shared_ptr<Position> position, double errBetaMomentum) {
  this->direction = direction;
  this->position = position;
  this->beta = beta;
  savedBeta = beta;
  this->momentum = 0;
  this->errBetaMomentum = errBetaMomentum;
  savedMomentum = momentum;

  //TEST CODE
  betaList.reserve(10000);
  mList.reserve(10000);
  xList.reserve(10000);
  yList.reserve(10000);
  tList.reserve(10000);
  pList.reserve(10000);
}

//DO NOT USE PREVIOUS TEST CODE
//
void Beam::setRandom() {
  position->setRandom();
  direction->setRandom();
  std::shared_ptr<TRandom3> randomGenerate(std::make_shared<TRandom3>());
  randomGenerate->SetSeed(0);
  beta = savedBeta + savedBeta * randomGenerate->Gaus(0.0, errBetaMomentum);
  momentum = savedMomentum + savedMomentum * randomGenerate->Gaus(0.0, errBetaMomentum);
  //std::cout << position->getX() << "   ";
  //std::cout << position->getY() << "   ";
  //std::cout << direction->getTheta() << "   ";
  //std::cout << direction->getPhi() << "   ";
  //std::cout << beta << "    ";
  //std::cout << momentum << std::endl;

  //TEST CODE
  betaList.push_back(this->getBeta());
  mList.push_back(this->getMomentum());
  xList.push_back(position->getX());
  yList.push_back(position->getY());
  tList.push_back(direction->getTheta());
  pList.push_back(direction->getPhi());
}

void Beam::plotAll(std::string particle = "", std::string path = "") {
  std::vector<std::string> variables;
  variables.push_back("beta");
  variables.push_back("momentum");
  variables.push_back("x");
  variables.push_back("y");
  variables.push_back("t");
  variables.push_back("p");

  std::vector<std::vector<double>> lists;
  lists.reserve(6);
  lists.push_back(betaList);
  lists.push_back(mList);
  lists.push_back(xList);
  lists.push_back(yList);
  lists.push_back(tList);
  lists.push_back(pList);

  for (int v = 0; v < lists.size(); v++) {
    std::string name;
    name = path + variables[v] + "RandomThrows" + particle;
    TCanvas *c = new TCanvas(name.c_str());
    TH1D *bH = new TH1D(name.c_str(), name.c_str(), 100, 0, 0);
    for (int i = 0; i < lists[v].size(); i++) {
      bH->Fill(lists[v][i]);
    }
    bH->Draw();
    name = name + ".png";
    c->SaveAs(name.c_str());
    //c->Delete();
  }
}

std::shared_ptr<Position> Beam::getPosition() {
  return position;
}

std::shared_ptr<Direction> Beam::getDirection() {
  return direction;
}

double Beam::getMomentum() {
  if (momentum != 0) {
    return momentum;
  } else {
    return particle * beta / (sqrt(1 - pow(beta, 2)));
  }
}

void Beam::setParticleMass(double particlemass) {
  this->particle = particlemass;
}

double Beam::getParticleMass() {
  return particle;
}

double Beam::getBeta() {
  if (this->beta == 0) {
    return momentum / sqrt(pow(momentum, 2) + pow(particle, 2));
  }
  return beta;
}

void Beam::getProjectedXYZ(double *x, double *y, double z, double cAngle, double phi) {
  double t = direction->getTheta();
  double p = direction->getPhi();
  //std::cout << "TEST XYZ 1 " << std::endl;
  double cos_t = cos(t);
  double cos_p = cos(p);
  double sin_t = sin(t);
  double sin_p = sin(p);

  TMatrixD R = TMatrixD(3, 3);
  R(0, 0) = 1 + (cos_t - 1) * pow(cos_p, 2);
  R(0, 1) = (cos_t - 1) * (sin_p * cos_p);
  R(0, 2) = -sin_t * cos_p;
  R(1, 0) = (cos_t - 1) * sin_p * cos_p;
  R(1, 1) = 1 + (cos_t - 1) * pow(sin_p, 2);
  R(1, 2) = -sin_t * sin_p;
  R(2, 0) = sin_t * cos_p;
  R(2, 1) = sin_t * sin_p;
  R(2, 2) = cos_t;

  //std::cout << z << std::endl;

  //std::cout << "TEST XYZ 2 " << std::endl;
  TVectorD v = TVectorD(3);
  v(0) = sin(cAngle) * cos(phi);
  v(1) = sin(cAngle) * sin(phi);
  v(2) = cos(cAngle);

  //std::cout << "TEST XYZ 3 " << std::endl;
  v = R.InvertFast() * v;
  //std::cout << "TEST XYZ 4 " << std::endl;
  double L = z * cos(t);
  *x = L * v(0) / v(2) + position->getX();
  *y = L * v(1) / v(2) + position->getY();
  //std::cout <<  position->getX() << "    " <<  position->getY() << std::endl;
  //std::cout << "TEST XYZ 5 " << std::endl;
}

double Beam::getLengthOfVector(double z) {
  return z / cos(direction->getTheta());
}
