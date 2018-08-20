#include "TH2D_SMOOTHING_CUSTOM.h"
#include "pmtClass.h"

PMT::PMT(double width, double length, int pixelCount) {
  this->width = width;
  this->length = length;
  this->pixelCount = pixelCount;
  this->windowThickness = 0.0;
  this->deadSpace = 0.0;
  this->timeResolution = 1E-9;            //s
  this->darkNoise = 8.0 * timeResolution; //counts
  this->setEnergyAndQuantumEfficiency();
  this->setCrossTalkValues();
}

double PMT::getWidth() {
  return width;
}

double PMT::getLength() {
  return length;
}

int PMT::getPixelCount() {
  return pixelCount;
}

void PMT::setWindowThickness(double wT) {
  this->windowThickness = wT;
}

double PMT::getWindowThickness() {
  return windowThickness;
}

void PMT::setDeadSpace(double dS) {
  this->deadSpace = dS;
}

double PMT::getDeadSpace() {
  return deadSpace;
}

void PMT::setDarkNoise (double dN) {
  this->darkNoise = dN;
}

double PMT::getDarkNoise() {
  return darkNoise;
}

void PMT::setTimeResolution(double tR) {
  this->timeResolution = tR;
}

double PMT::getTimeResolution() {
  return timeResolution;
}

void PMT::setPosition(std::shared_ptr<Position> pos) {
  this->position = pos;
}

std::shared_ptr<Position> PMT::getPosition() {
  return position;
}

int PMT::getLengthPixelCount() {
  return sqrt(pixelCount);
}

int PMT::getWidthPixelCount() {
  return sqrt(pixelCount);
}

std::vector<double> PMT::getEnergy() {
  double h = 4.135E-15; //eV*s
  double c = 299792458; //m/s
  std::vector<double> energy;
  energy.reserve(wavelength.size());
  for (int i = 0; i < wavelength.size(); i++) {
    energy.push_back(h * c / (1E-9 * wavelength.at(i)));
    //std::cout <<energy.back() << std::endl;
  }
  return energy;
}

void PMT::setEnergyAndQuantumEfficiency() {
  //std::string symFile ="quantumEfficiency.csv";
  //std::string absolutePath = realpath(symFile.c_str(),NULL);
  //std::cout << absolutePath.c_str() << std::endl;
  std::ifstream infile("quantumEfficiency.csv");
  //std::cout << "TEST" << std::endl;
  if (infile.is_open()) {
    int i = 0;
    while (infile.good()) {
      std::string line;
      getline(infile, line);
      std::string delimeter = ",";
      int pos = line.find(delimeter);
      std::string temp = line.substr(0, pos);
      line.erase(0, pos + delimeter.size());
      wavelength.push_back(stod(temp));
      quantumEfficiency.push_back(stod(line));
    }
    infile.close();
  } else {
    std::cout << "ERROR NO QUANTUM EFFICIENCY FILE" << std::endl;
  }
  return;
}

std::vector<double> PMT::getWavelength() {
  return wavelength;
}

std::vector<double> PMT::getQuantumEfficiency() {
  return quantumEfficiency;
}

std::shared_ptr<TSpline3> PMT::getQuantumEfficiencySpline() {

  std::vector<double> x = this->wavelength;
  std::vector<double> y = quantumEfficiency;
  std::vector<double> x_d(x.begin(),x.end());
  std::vector<double> y_d(y.begin(),y.end());
  std::shared_ptr<TGraphErrors> gr(std::make_shared<TGraphErrors>(x_d.size(), &(x_d[0]), &(y_d[0])));
  gr->Sort();
  std::shared_ptr<TSpline3> grs(std::make_shared<TSpline3>("grs", gr.get()));
  return grs;
}

void PMT::setCrossTalkValues() {
  std::ifstream infile("crosstalkMask.txt");
  //std::cout << "TEST" << std::endl;
  if (infile.is_open()) {
    while (infile.good()) {
      std::string line;
      getline(infile, line);
      if (line.find("#") == std::string::npos) {
        std::istringstream iss(line);
        std::string name;
        double value;
        iss >> name >> value;
        if (name.find("neighbor") != std::string::npos) {
          mask.neighbor = value / 100.0;
        } else if (name.find("diagonal") != std::string::npos) {
          mask.diagonal = value / 100.0;
        }
      }
    }
    infile.close();
  } else {
    mask.neighbor = double(1) / 100;
    mask.diagonal = double(.1) / 100;
    std::cout << "ERROR NO CROSSTALK FILE" << std::endl;
  }
}

bool PMT::checkCrossTalk(int index) {
  return ((index < getPixelCount()) && (index >= 0));
}

std::vector<double> PMT::applyCrossTalkMask(std::vector<double> hitsPerPixel) {
  std::vector<double> temp;
  temp.reserve(hitsPerPixel.size());
  double N = getWidthPixelCount();
  for (int i = 0; i < hitsPerPixel.size(); i++) {
    //std::vector<bool> applyMask;
    //std::vector<double> index;
    double totalValue = 0;
    for (int j = -1; j < 2; j++) {
      for (int k = -1; k < 2; k++) {
        if (int(i + k * N + j) % 2 == 0 && checkCrossTalk(i + k * N + j)) {
          totalValue = totalValue + mask.diagonal * hitsPerPixel.at(i + k * N + j);
        } else if (int(i + k * N + j) % 2 != 0 && checkCrossTalk(i + k * N + j)) {
          totalValue = totalValue + mask.neighbor * hitsPerPixel.at(i + k * N + j);
        }
      }
    }
    temp.push_back(floor(hitsPerPixel.at(i) + totalValue) + 1E-20);
    //totalValue = 0;
  }
  return temp;
}

std::vector<double> PMT::applyCrossTalkTH2DMask(std::vector<double> pixelMeans){
  
	TH2D * h2 = new TH2D("h2","",10,-5,5,10,-5,5);
	for( int pixel = 0; pixel <= pixelMeans.size(); pixel++)
	{
		h2->SetBinContent(pixel,pixelMeans[pixel]);
	}
	TH2D_Custom::Smooth(h2,"k3a");
	std::vector<double> crossTalkMeans;
	crossTalkMeans.reserve(pixelMeans.size());
	for (int pixel = 0; pixel <= pixelMeans.size();pixel++)
	{
		crossTalkMeans.push_back(h2->GetBinContent(pixel,pixelMeans[pixel]));
	}	
	return crossTalkMeans;
}

double PMT::getWidthPixel() {
  return width / getWidthPixelCount();
}

double PMT::getLengthPixel() {
  return length / getLengthPixelCount();
}
