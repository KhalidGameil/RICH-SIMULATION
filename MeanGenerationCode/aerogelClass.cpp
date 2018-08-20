#include "aerogelClass.h"

Aerogel::Aerogel(double width, double length, double thickness, double refractiveIndex) {
  randomGenerate=std::make_shared<TRandom3>();
  randomGenerate->SetSeed(0);
  this->width = width;
  this->length = length;
  this->thickness = thickness;
  this->refractiveIndex = refractiveIndex;
  this->setPosition(std::make_shared<Position>(0.0, 0.0, 0.0)); //Set position to Zero for testing
}

double Aerogel::getWidth() {
  return width;
}

double Aerogel::getLength() {
  return length;
}

double Aerogel::getThickness() {
  return thickness;
}

double Aerogel::getRefractiveIndex() {
  return refractiveIndex;
}

void Aerogel::setPosition(std::shared_ptr<Position> pos) {
  this->position = pos;
}

std::shared_ptr<Position> Aerogel::getPosition() {
  return position;
}

std::vector<double> Aerogel::calcRefractiveIndexEnergy(std::vector<double> energy) {
  std::vector<double> n(energy.size(), this->refractiveIndex);
  return n;
}

std::vector<double> Aerogel::calcSellmeierRefractiveIndex(std::vector<double> energy) {
  double h = 4.136E-15; //eVs
  double c = 2.998E8;   //m/s
  std::vector<double> n;
  double a_0 = 0.0;
  double wvl_0 = 0.0;
  for (int i = 0; i < energy.size(); i++) {
    double wvl = h * c / (energy.at(i));
    n.push_back(sqrt(a_0 * pow(wvl, 2) / (pow(wvl, 2) - pow(wvl_0, 2)) + 1));
  }
  return n;
}

double Aerogel::calcN(std::shared_ptr<Beam> beam, std::vector<double> energy) {
  double beta = beam->getBeta();
  double dx = beam->getLengthOfVector(this->thickness) * 100; //has to be in cm for the equation below
  double z = 1.0;
  double dNdX = 0.0;
  std::vector<double> n = calcRefractiveIndexEnergy(energy);
  double dEnergy = energy.at(1) - energy.at(0);
  //std::cout << dx << "   " << beta << "   " << dEnergy << "   " << std::endl;
  for (int i = 0; i < n.size(); i++) {
    dNdX = dNdX + dx * (370 * pow(z, 2) * (1.0 - 1.0 / pow((n.at(i) * beta), 2))) * (dEnergy);
  }
  return dNdX;
}

std::vector<double> Aerogel::calcEnergyPdf(double beta, std::vector<double> energy) {
  std::vector<double> n = calcRefractiveIndexEnergy(energy);
  std::vector<double> pdf;
  pdf.reserve(n.size());
  for (int i = 0; i < n.size(); i++) {
    pdf.push_back((1 - 1.0 / (pow(n.at(i) * beta, 2))));
  }
  return pdf;
}

std::vector<double> Aerogel::calcEnergyCdf(double beta, std::vector<double> energy) {
  std::vector<double> n = calcRefractiveIndexEnergy(energy);
  std::vector<double> pdf = calcEnergyPdf(beta, energy);
  std::vector<double> cdf;
  cdf.reserve(pdf.size());
  double norm = std::accumulate(pdf.begin(), pdf.end(), 0.0);

  //std::cout << "NORM" << norm << std::endl;
  for (int i = 0; i < pdf.size(); i++) {
    std::vector<double> pdfToI(&pdf[0], &pdf[i]);
    cdf.push_back(std::accumulate(pdfToI.begin(), pdfToI.end(), 0.0) / norm);
    //std::cout << norm << "   " <<pdfToI.size() << "   " << cdf.back() << std::endl;
  }
  cdf.back() = 1.0;
  return cdf;
}

std::vector<double> Aerogel::findCdf(std::vector<double> r, std::vector<double> energy, std::vector<double> cdf) {
  std::vector<double> x_r;
  x_r.reserve(r.size());
  //std::cout << r.size() << "   " << energy.size() << "   " << cdf.size() << std::endl;
  //std::cout << cdf.back() << std::endl;
  for (int i = 0; i < r.size(); i++) {
    //std::vector<double> temp(energy.size(),0.0);
    //std::transform(cdf.begin(), cdf.end(), myvec.begin(),
    //bind2nd(std::minus<double>(), r.at(i)))
    //std::cout << i << "    "<< r.at(i) <<"  ";
    int eIndex = std::upper_bound(cdf.begin(), cdf.end(), r.at(i)) - cdf.begin();
    //std::cout << eIndex <<"    " << energy.at(eIndex) <<std::endl;
    x_r.push_back(energy.at(eIndex));
    //std::cout << "TEST2: " << x_r.back() << std::endl;
  }
  return x_r;
}

std::vector<double> Aerogel::getRandomInvCDF(int N, double beta, std::vector<double> energy, std::vector<double> e_cdf) {
  //std::shared_ptr<TRandom3> randomGenerate(std::make_shared<TRandom3>());
  //randomGenerate->SetSeed(1);

  std::vector<double> random_cdfValue;
  random_cdfValue.reserve(N);
  //std::cout << "#########N#######    " << N << std::endl;
  for (int i = 0; i < N; i++) {
    double rG = randomGenerate->Uniform(0.0, 1.0);
    //std::cout << i << "   " << rG << std::endl;
    random_cdfValue.push_back(rG);
  }
  return findCdf(random_cdfValue, energy, e_cdf);
}

std::vector<double> Aerogel::getRandomInteractionPosition(int N, std::shared_ptr<Beam> beam, double tD) {

  std::vector<double> interactionPosition;
  interactionPosition.reserve(N);
  double start = beam->getLengthOfVector(position->getZ());
  double end = start + beam->getLengthOfVector(thickness);
  //std::shared_ptr<TRandom3> randomGenerate(std::make_shared<TRandom3>());
  //randomGenerate->SetSeed(0);
  tD = beam->getLengthOfVector(tD);
  for (int i = 0; i < N; i++) {
    //std::cout << tD << "   " << start << "  "  << end << "  " << tD - (start+end)/2.0 << std::endl;
    interactionPosition.push_back(tD - randomGenerate->Uniform(start, end));
  }
  return interactionPosition;
}

std::vector<double> Aerogel::getRandomTheta(int N) {
  std::vector<double> theta;
  theta.reserve(N);
  //std::shared_ptr<TRandom3> randomGenerate(std::make_shared<TRandom3>());
  //randomGenerate->SetSeed(0);
  for (int i = 0; i < N; i++) {
    theta.push_back(randomGenerate->Uniform(0.0, 2 * TMath::Pi()));
  }
  return theta;
}

double Aerogel::getCherenkovAngle(double beta) {
  return acos(1.0 / (refractiveIndex * beta));
}

std::vector<double> Aerogel::getRandomEnergy(int N, double energyMin, double energyMax) {
  std::vector<double> energy;
  energy.reserve(N);
  //std::shared_ptr<TRandom3> randomGenerate(std::make_shared<TRandom3>());
  //randomGenerate->SetSeed(1);
  for (int i = 0; i < N; i++) {
    energy.push_back(randomGenerate->Uniform(energyMin, energyMax));
  }
 
  return energy;
}

std::shared_ptr<RandomEvent> Aerogel::getEvent(double totalDistance, std::shared_ptr<Beam> beam,
                                               std::vector<double> energy, std::shared_ptr<TSpline3> spl) {
  //beamX = 0.0;
  //beamY = 0.0;
  //std::cout << "TEST" << std::endl;
  //randomGenerate(std::make_shared<TRandom3>());
  //randomGenerate=std::make_shared<TRandom3>();
  //randomGenerate->SetSeed(0);
  //beam->setRandom();

  //*** REMOVED FOR SPEED TEST***
  
  std::clock_t begin_single = std::clock();
  double beta = beam->getBeta();
  int randomN = randomGenerate->Poisson(calcN(beam, energy)); //Number of photons
  double chernAngle = getCherenkovAngle(beam->getBeta()); //removed for speed

  //REMOVED FOR CALCULATION TIME
  //std::vector<double> e_cdf = this->calcEnergyCdf(beta,energy);
  //double emin = *std::min_element(energy.begin(), energy.end());
  //double emax = *std::max_element(energy.begin(), energy.end());

  //*** REMOVED FOR SPEED TEST***
  std::vector<double> photonEnergy = getRandomEnergy(randomN, energy.front(), energy.back()); //getRandomInvCDF(randomN,beta,energy,e_cdf);
  std::vector<double> interactionPosition = getRandomInteractionPosition(randomN, beam, totalDistance);
  std::vector<double> theta = getRandomTheta(randomN);

  //SPEED TEST
  std::shared_ptr<RandomEvent> event(std::make_shared<RandomEvent>(randomN, 
								   chernAngle,//getCherenkovAngle(beam->getBeta()),
                                                                   interactionPosition,//getRandomInteractionPosition(randomN, beam, totalDistance), 
								   theta,//getRandomTheta(randomN),
                                                                   photonEnergy,//getRandomEnergy(randomN, energy.front(), energy.back()), 
								   beam, spl));
 std::clock_t end_single = std::clock();
  double elapsed_secs_single =  double(end_single- begin_single) / CLOCKS_PER_SEC;

  return event;
}
