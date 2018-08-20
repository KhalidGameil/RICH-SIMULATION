#include "detectorClass.h"

Detector::Detector(std::vector<std::shared_ptr<Aerogel>> aerogels, std::shared_ptr<PMT> pmt,
                   double aeroDetDistance) {
  this->aerogels = aerogels;
  this->pmt = pmt;
  this->aeroDetDistance = aeroDetDistance; // Distance from Final Aerogel to face of the multianodepmt
  this->setCentralPositions();

  double eCount = 1000;
  energyRangePMT = getGeneratedEnergyRange(eCount, pmt->getEnergy()); //calculate the new range of the pmt
  
  plotCount = 0;
  quickReadWParameter();
}

//Detector::Detector(Detector *d) {
//this = d;
//}

/*
 * Private Getters
 */
double Detector::getMaxWidth() {
  std::vector<double> widths;
  widths.reserve(aerogels.size());
  for (int i = 0; i < aerogels.size(); i++) {
    widths.push_back(aerogels.at(i)->getWidth());
  }
  widths.push_back(pmt->getWidth());
  return *std::max_element(widths.begin(), widths.end());
}
double Detector::getMaxLength() {
  std::vector<double> lengths;
  lengths.reserve(aerogels.size());
  for (int i = 0; i < aerogels.size(); i++) {
    lengths.push_back(aerogels.at(i)->getLength());
  }
  lengths.push_back(pmt->getWidth());
  return *std::max_element(lengths.begin(), lengths.end());
}

/*
 *Simple (some calculation involved) Getters
 */
std::vector<std::shared_ptr<Aerogel>> Detector::getAerogelArray() { return this->aerogels; } //return List of Aerogels In detector
std::shared_ptr<PMT> Detector::getPMT() { return this->pmt; }                                //return pmt
double Detector::getAeroDetDistance() { return this->aeroDetDistance; }                      // return distance from aerogel and detector
double Detector::getTotalDistance() {
  /*
 *Calculate the total distance between the first aerogel to the face of detector
 *Input: 
 *Output: totalDistance 
 */
  double totalDistance = 0.0;
  for (int i = 0; i < aerogels.size(); i++) { //loops through list of aerogels
    totalDistance = totalDistance + aerogels.at(i)->getThickness();
  }
  totalDistance = totalDistance + aeroDetDistance;
  return totalDistance;
}

/*
 * Add Object
 */
void Detector::addAerogel(std::shared_ptr<Aerogel> aerogel) { this->aerogels.push_back(aerogel); } //add Aergeol to list

/*
 * Setters
 */
void Detector::setCentralPositions() {
  /*
 * Using all the position information from all the objectggs (aerogels, pmts, etc?) 
 * calculate and set the central most position for all the objects
 * Input: 
 * Output: 
 */
  double maxW = getMaxWidth();  //highest max width
  double maxL = getMaxLength(); //highest max length
  double zPosition = 0.0;

  for (int i = 0; i < this->aerogels.size(); i++) {                          //Iterate the list of aerogels
    double xPosition = maxW / 2.0 - this->aerogels.at(i)->getWidth() / 2.0;  //calculate the new x position
    double yPosition = maxL / 2.0 - this->aerogels.at(i)->getLength() / 2.0; //calculate the new y position
    std::shared_ptr<Position> position(std::make_shared<Position>(xPosition, yPosition, zPosition));
    this->aerogels.at(i)->setPosition(position);
    zPosition = zPosition + this->aerogels.at(i)->getThickness();
  }
  double xPmt = maxW / 2.0 - pmt->getWidth() / 2.0;
  double yPmt = maxL / 2.0 - pmt->getLength() / 2.0;
  zPosition = zPosition + aeroDetDistance;
  std::shared_ptr<Position> aeroPos(std::make_shared<Position>(xPmt, yPmt, zPosition));
  pmt->setPosition(aeroPos);
}

/*)
 * Event Generators
 */
std::vector<std::shared_ptr<RandomEvent>> Detector::getEventFromBeam(std::shared_ptr<Beam> beam) {
  /*
 * Wrapper function for the aerogel getEvent function and pass into the beam class variables
 * Input: beam  = beam object used to get events
 * Output: vector of random event objects per pixel (size might be 256 or (8n)^2)
 */
  std::vector<std::shared_ptr<RandomEvent>> eventsForAerogels;
  eventsForAerogels.reserve(aerogels.size());
  //std::vector<double> tempEnergy = pmt->getEnergy(); //get energy range of the pmt

  std::vector<RandomEvent> event;
  double tD = getTotalDistance();             //total distance
  for (int i = 0; i < aerogels.size(); i++) { //iterate over all aerogel
    eventsForAerogels.push_back(aerogels.at(i)->getEvent(tD, beam, energyRangePMT, pmt->getQuantumEfficiencySpline()));
  }
  return eventsForAerogels;
}

std::vector<std::vector<std::shared_ptr<RandomEvent>>> Detector::getMultipleEventsFromBeam(std::shared_ptr<Beam> beam, int eventCounts = 1) {
  /*
 * Wrapper function for getEventFromBeam to iterate over the multiple runs of the beams
 * Input: beam  = beam object used to get events
 * Ouptut: vector of random event objects per pixel (size might be 256 or (8n)^2) per event (10,000 or eventCounts)
 */
  std::vector<std::vector<std::shared_ptr<RandomEvent>>> rE;
  rE.reserve(eventCounts);
  for (int event = 0; event < eventCounts; event++) {
    rE.push_back(getEventFromBeam(beam));
  }
  this->randomEvents = rE; //Stores randomEvents
  return rE;
}

std::vector<double> Detector::getGeneratedEnergyRange(int eCount, std::vector<double> tempEnergy) {
  /*
 * From the range of energies and the number of steps desired return a new vector containing those steps
 * Input: 
 *  eCount = the number of steps desired for the new vector
 *  tempEnergy = the input vector of energies
 * Output: a vector of energies the size of eCount
 */
  double emin = *std::min_element(tempEnergy.begin(), tempEnergy.end());
  double emax = *std::max_element(tempEnergy.begin(), tempEnergy.end());
  double de = (emax - emin) / (eCount - 1); //step size

  std::vector<double> energy;
  energy.reserve(eCount);
  double nextE = emin;
  for (int i = 0; i < eCount; i++) {
    energy.push_back(nextE);
    nextE = nextE + de;
  }
  return energy;
}

std::vector<std::vector<std::shared_ptr<RandomEvent>>> Detector::getRandomEvents() {
  /*
 * Returns the random events
 */
  return randomEvents;
}

std::vector<std::vector<std::shared_ptr<DetectedEvent>>> Detector::getDetectedEvents(std::shared_ptr<Beam> beam, int eventCounts = 1) {
 std::vector<std::vector<std::shared_ptr<RandomEvent>>> rE = getMultipleEventsFromBeam(beam, eventCounts);
  std::vector<std::vector<std::shared_ptr<DetectedEvent>>> dE;
  dE.reserve(rE.size());
  for (int i = 0; i < rE.size(); i++) {
    std::vector<std::shared_ptr<DetectedEvent>> temp;
    std::vector<std::shared_ptr<RandomEvent>> rEAerogels = rE.at(i);
    temp.reserve(rEAerogels.size());
    for (int j = 0; j < rEAerogels.size(); j++) {
	  temp.push_back(rEAerogels.at(j)->getDetectedHits());
    }
    dE.push_back(temp);
  }
  return dE;
}

void Detector::testWrite(std::vector<double> x, std::vector<double> y)
{

	std::string fileName="TEST_X_AND_Y_POSITIONS.txt";
	FILE *f = fopen(fileName.c_str(), "a");
	//assert(f); //does f exist
	for (int variable = 0; variable < x.size(); variable++) {
		fprintf(f, "%g ", double( x[variable]));
		fprintf(f, "%g ", double( y[variable]));
		fprintf(f,"\n");
	}
	fclose(f);

}

void Detector::getBinEdges(std::vector<double> *binEdgesx, std::vector<double> *binEdgesy) {
  // Calculated the binEdges given the PMT value
  // Input:
  // Pointer to vector of doubles for binEdgesx
  // Pointer to vector of doubles for binEdgesy
  double dW = pmt->getWidth();
  double dL = pmt->getLength();

  int pixelN = pmt->getPixelCount();
  int binEdgesN = pmt->getLengthPixelCount() + 1; // where getLengthPixelCount = getWidthPixelCount

  double xD = pmt->getPosition()->getX();
  double yD = pmt->getPosition()->getY();

  double deltaXD = dW / (binEdgesN - 1); //((xD+dW)-(xD))/N;
  double deltaYD = dL / (binEdgesN - 1); //((yD+dL)-(yD))/N;

  double totalX = xD;
  double totalY = yD;

  for (int i = 0; i < binEdgesN; i++) {
    binEdgesx->push_back(totalX);
    binEdgesy->push_back(totalY);
    totalX = totalX + deltaXD;
    totalY = totalY + deltaYD;
  }
}

std::vector<std::vector<double>> Detector::getPixelHitsPerEvents(
    std::vector<std::vector<std::shared_ptr<DetectedEvent>>> events,
    bool testSignal = false) {
  // Calculates the 2D vector of Histogram hits per pixel for each event
  // INPUT:
  // 2D Vector of Detected Hits
  // PMT Object
  double pixelN = pmt->getPixelCount();
  std::vector<double> binEdgesx;
  std::vector<double> binEdgesy;
  getBinEdges(&binEdgesx, &binEdgesy);

  // Generates random uniform noise
  //std::shared_ptr<TRandom3> randomgenerate(std::make_shared<TRandom3>());
  //randomgenerate->SetSeed(1);

  // Loops through all the the events (for each aerogel)
  // Then "merges" all the x values/y values for each aerogel (per an event)
  // Then (for that event) calculates the histogram values of x and y
  std::vector<std::vector<double>> hit_on_pixels_perEvent;
  double maxX = 0;
  double meanX = 0;
  double upperMeanX = 0;
  int countMean = 0;
  int countUpperMean = 0;

    std::shared_ptr<TH2D> h2(std::make_shared<TH2D>("PMT BINNING HISTOGRAM", "PMT BINNING HISTOGRAM", binEdgesx.size() - 1, &binEdgesx[0],
        binEdgesy.size() - 1, &binEdgesy[0]));

  for (int beamRuns = 0; beamRuns < events.size(); beamRuns++) {
    std::vector<std::shared_ptr<DetectedEvent>> eventsForAerogels = events.at(beamRuns);
    double totalN = 0;
    //std::vector<double> x_perEvent;
    //std::vector<double> y_perEvent;
    //Merge events from each aerogels for future histogram sorting
    for (int aerogels = 0; aerogels < eventsForAerogels.size(); aerogels++) {
      std::shared_ptr<DetectedEvent> dE = eventsForAerogels.at(aerogels);
      std::vector<double> tempX = dE->getX();
      std::vector<double> tempY = dE->getY();
      //x_perEvent.insert(x_perEvent.end(), tempX.begin(), tempX.end());
      //y_perEvent.insert(y_perEvent.end(), tempY.begin(), tempY.end());
      if (!testSignal) //Normal signal
      {
	for (int hit = 0; hit < tempX.size(); hit++)
	{
	  h2->Fill(tempX[hit],tempY[hit]);
	}
      }
    }
    if (testSignal != true) {//NORMAL SIGNAL
      h2->Draw();
    } else {//TEST SIGNAL FOR DEBUGING
      int testSignal = 0;
      for (int j = 1; j < binEdgesy.size(); j++) {
        for (int k = 1; k < binEdgesx.size(); k++) {
          h2->SetBinContent(k, j, testSignal);
          testSignal++;
        }
      }
    }
    //END OF TEST SIGNAL

    // Creates a new vector out of each hit per event
    std::vector<double> pixelHits;
    pixelHits.reserve(pixelN);
    /*for (int j = 1; j < binEdgesy.size(); j++) {
      for (int k = 1; k < binEdgesx.size();
           k++) {
       if (true)//((j-1)*(binEdgesy.size()-1)+k-1 == 588)
       	  {//std::cout << h2->GetBinContent(k,j) << std::endl;
	   pixelHits.push_back(h2->GetBinContent(k, j));}
       else {pixelHits.push_back(0.0);}
	}
    }*/
    for (int j = 0; j < h2->GetSize(); j++)
    {
	pixelHits.push_back(h2->GetBinContent(j));
    }
    //pixelHits = pmt->applyCrossTalkMask(pixelHits); ##THIS IS A REMINDER THERE IS NO REASON TO PUT CROSSTALK IN THIS STEP
    hit_on_pixels_perEvent.push_back(pixelHits);
    h2->Reset("ICESM");
  }
  return hit_on_pixels_perEvent;
}

std::vector<double> Detector::getSumOfHitsPerPixel(std::vector<std::vector<double>> hit_on_pixels_perEvent) {
  std::vector<double> sums(hit_on_pixels_perEvent.at(0).size());
  for (int i = 0; i < hit_on_pixels_perEvent.size(); i++) {
    std::vector<double> pHRow = hit_on_pixels_perEvent.at(i);
    for (int j = 0; j < pHRow.size(); j++) {
      sums.at(j) = sums.at(j) + pHRow.at(j);
    }
  }
  return sums;
}

std::vector<std::shared_ptr<TH1D>> Detector::getArrayOfHistograms(std::vector<std::vector<double>> pH) {
  std::vector<std::shared_ptr<TH1D>> hists;
  std::vector<double> currentMean;
  hists.reserve(pH.at(0).size());
  ////std::cout << pH.at(0).size() << std::endl;
  for (size_t colNo = 0; colNo < pH.at(0).size(); colNo++) {
    std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>(
        (std::to_string(colNo) + "Number" + std::to_string(plotCount)).c_str()));
    std::shared_ptr<TH1D> h1(std::make_shared<TH1D>(
        (std::to_string(colNo) + "Number" + std::to_string(plotCount)).c_str(),
        ("Histogram of pixel" + std::to_string(colNo)).c_str(), 10, 0, 10));
    for (const auto &row : pH) {
      h1->Fill(row.at(colNo));
    }
    h1->SetStats();
    h1->SetFillColor(38);
    h1->SetXTitle("Number of Detected Photons");
    h1->SetYTitle("Counts");
    // TF1 * f1 =  new TF1("pois","[0]*TMath::Poisson(x,[1])",0,10);
    // h1->Fit("pois");
    double area = h1->Integral();
    if (area != 0) {
      h1->Scale(1 / area);
    }
    h1->Draw();
    currentMean.push_back(h1->GetMean());
    hists.push_back(h1);
  }
  mean = currentMean;
  // //std::cout <<  hists.size() << std::endl;

  plotCount++;
  return hists;
}

std::vector<std::shared_ptr<TGraph>> Detector::getArrayOfExpectedPoissonHistogramsPerPixel(std::vector<std::vector<double>> pH) {
  std::vector<std::shared_ptr<TGraph>> graphs;
  graphs.reserve(pH.at(0).size());

  // Create the increments required to fill the dependent variable for the
  // poisson distribution
  double start = 0.0;
  double end = 10.0;
  double steps = 100;
  double delta = (end - start) / (steps - 1);

  for (int i = 0; i < pH.at(0).size(); i++) {
    std::vector<double> x;
    std::vector<double> y;
    start = 0.0;
    for (int j = 0; j < steps; j++) {
      x.push_back(start);
      start = start + delta;
      y.push_back(TMath::Poisson(x.back(), mean.at(i)));
    }
    std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>(
        (std::to_string(i) + "Number" + std::to_string(plotCount)).c_str()));
    std::shared_ptr<TGraph> g1(std::make_shared<TGraph>(x.size(), &x[0], &y[0]));
    g1->Draw();
    graphs.push_back(g1);
  }

  // plotCount++;
  return graphs;
}

/*
 * Log Likelihood Function
 */

//std::vector<double> Detector::getLogLikelihoodWithNuisance(std::vector<std::vector<double>> pH,
//std::vector<double> mean1,
//std::vector<double> mean2,
//bool minuit = false) {

//std::clock_t begin = std::clock();
//std::vector<double> ll_pH;
//for (int rowNo = 0; rowNo < pH.size(); rowNo++) {
//std::vector<double> pH_event = pH.at(rowNo);
//double ll_event = 0.0;
//for (int colNo = 0; colNo < pH_event.size(); colNo++) {
//double w1;
//double sig1 = .003502;
//if (!minuit) {
//w1 = getNuisanceParameter(pH_event.at(colNo), mean1.at(colNo), sig1, pH.size());
////std::cout << "w1: " << w1 << "    ";
//} else {
//w1 = getNuisanceParameterMINUIT(pH_event.at(colNo), mean1.at(colNo), sig1, pH.size());
////std::cout << "w1n: " << w1 << std::endl;
//}
//double w2;
//double sig2 = .002851;
//if (!minuit) {
//w2 = getNuisanceParameter(pH_event.at(colNo), mean2.at(colNo), sig2, pH.size());
////std::cout << "w2: " << w2 << "    ";
//} else {
//w2 = getNuisanceParameterMINUIT(pH_event.at(colNo), mean2.at(colNo), sig2, pH.size());
////std::cout << "w2n: " << w2 << std::endl;
//}
//ll_event = ll_event - 2 * calcLikelihoodForEvent(pH_event.at(colNo), mean1.at(colNo) + w1) + pow(w1 / sig1, 2) +
//2 * calcLikelihoodForEvent(pH_event.at(colNo), mean2.at(colNo) + w2) - pow(w2 / sig2, 2);
//}
//ll_pH.push_back(ll_event);
//}
//std::clock_t end = std::clock();
//double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//std::cout << "TIME THE CODE TOOK NUISANCE:" << elapsed_secs << std::endl;
//return ll_pH;
//}

void Detector::setPixelHit(std::vector<double> pH) {
  pixelHits = pH;
}

std::vector<double> Detector::getPixelHit() {
  return pixelHits;
}

void Detector::setMean(std::vector<double> mean) {
  this->mean = mean;
}

std::vector<double> Detector::getMean() {
  return mean;
}

void Detector::setPixelHitSingleEvent(double pH) {
  this->pHSingleEvent = pH;
}

double Detector::getPixelHitSingleEvent() {
  return pHSingleEvent;
}

void Detector::setMeanSingleEvent(double mean) {
  this->mSingleEvent = mean;
}

double Detector::getMeanSingleEvent() {
  return mSingleEvent;
}

void FCN(Int_t &npar, Double_t *gin, Double_t &f,
         Double_t *u, Int_t flag) {
  Detector *d = (Detector *)gMinuit->GetObjectFit();
  //f = d->calcLikelihoodForPixelNuisance(d->getPixelHit(), d->getMean(), u[0], u[1]);
  f = -2 * d->calcLikelihoodForEventNuisance(d->getPixelHitSingleEvent(), d->getMeanSingleEvent(), u[0], u[1]);
  //std::cout << u[0] << "    ";
  //std::cout << f << std::endl;
}

void Detector::setVariables(std::vector<double> v) {
  variables = v;
}
std::vector<double> Detector::getVariables() {
  return variables;
}

void FCN2(Int_t &npar, Double_t *gin, Double_t &f,
          Double_t *u, Int_t flag) {
  Detector *d = (Detector *)gMinuit->GetObjectFit();
  std::vector<double> variables = d->getVariables();
  f = 0;
  for (int i = 0; i < variables.size(); i++) {
    f = f + ((-variables[i] + u[i]) * (-variables[i] + u[i])) / (u[i + variables.size()] * u[i + variables.size()]);
  }
  f = f + u[10];
}

void FCN3(Int_t &npar, Double_t *gin, Double_t &f,
          Double_t *u, Int_t flag) {

  Detector *d = (Detector *)gMinuit->GetObjectFit();
  std::vector<double> variables = d->getVariables();
  std::shared_ptr<InterpolateLikelihood> i = d->getInterpolator();
  f = 0;
  std::vector<double> pH = d->getPixelHit();

  std::vector<double> variablesTrue;
  for (int v = 0; v < variables.size(); v++) {
    variablesTrue.push_back(u[v]);
    //std::cout << variablesTrue[v] << "   ";
  }
  std::vector<bool> variableSwitch;
  variableSwitch.push_back(false);
  variableSwitch.push_back(false);
  variableSwitch.push_back(false);
  variableSwitch.push_back(false);  
  variableSwitch.push_back(true);
  std::vector<double> vT_double(variables.begin(),variables.end());
  std::vector<double> mean = i->interpolationSwitch(variablesTrue, 2, variableSwitch,
                                                    d->getPMT()->getWidthPixelCount(), d->getTotalDistance(),
                                                    d->getMaxWidth(), d->getMaxLength(),
                                                    d->getPMT()->getWidthPixel(), d->getPMT()->getLengthPixel());
  std::shared_ptr<TSpline3> s = d->getWSpline();

  for (int pixelCount = 0; pixelCount < pH.size(); pixelCount++) {
    double m = mean[pixelCount];
    double w = s->Eval(m);
    f = f - 2 * d->calcLikelihoodForEventNuisance(pH[pixelCount], m, w, u[10]);
    //std::cout << pixelCount << "    " << pH[pixelCount] << "    ";
    //std::cout << m << "  " << w << "    " << u[10] << "   " << f << std::endl;
  }
  //std::cout << f << "    ";
  for (int v = 0; v < variables.size(); v++) {
    f = f + ((-variables[v] + u[v]) * (-variables[v] + u[v])) / (u[v + variables.size()] * u[v + variables.size()]);
  }
  //std::cout << f << std::endl;
}

std::vector<double> Detector::getLogLikelihoodWithNuisanceParameters(std::vector<std::vector<double>> pH,
                                                                     std::vector<double> mean1,
                                                                     std::vector<double> mean2,
                                                                     std::shared_ptr<Beam> beam1 = NULL,
                                                                     std::shared_ptr<Beam> beam2 = NULL,
                                                                     std::shared_ptr<InterpolateLikelihood> i = NULL,
								     int iType = 0, bool wSwitch = false,
                                                                     std::string name = "") {
  //std::cout << (variables1.size() == 0 || variables2.size() == 0) << std::endl;
  //std::cout << bool(wSpline == NULL) << std::endl;
  if ((beam1 == NULL || beam2 == NULL) && wSpline == NULL) {
    std::cout << "NO SPLINE NO 5 PAR" << std::endl;
    return getLogLikelihoodWithNuisanceMinuit(pH, mean1, mean2);
  } else if ((beam1 == NULL || beam2 == NULL) && wSpline != NULL) {
    std::cout << "SPLINE" << std::endl;
    double sig1 = .00315; //These are taken to be constants
    double sig2 = .002851;

    std::vector<double> ll;
    for (int eventCount = 0; eventCount < pH.size(); eventCount++) {
      double ll1, ll2 = 0;
      for (int pixelCount = 0; pixelCount < pH[eventCount].size(); pixelCount++) {
        if (pH[eventCount][pixelCount] >= 1) {
          ll1 = ll1 + -2 * calcLikelihoodForEventNuisance(pH[eventCount][pixelCount], mean1[pixelCount], wSpline->Eval(mean1[pixelCount]), sig1);
          ll2 = ll2 + -2 * calcLikelihoodForEventNuisance(pH[eventCount][pixelCount], mean2[pixelCount], wSpline->Eval(mean2[pixelCount]), sig2);
        } else {
          double w = -pow(sig1, 2);
          ll1 = ll1 + -2 * calcLikelihoodForEventNuisance(pH[eventCount][pixelCount], mean1[pixelCount], w, sig1);
          w = -pow(sig2, 2);
          ll2 = ll2 + -2 * calcLikelihoodForEventNuisance(pH[eventCount][pixelCount], mean2[pixelCount], w, sig2);
        }
      }
      ll.push_back(ll2 - ll1);
    }
    return ll;
  } else {
    std::cout << "5 variable" << std::endl;
    return getLogLikelihoodForVariableAndWparameter(pH, mean1, mean2, beam1, beam2, i,iType,wSwitch, name);
  }
}

std::vector<double> Detector::getLogLikelihoodWithNuisanceParameters(std::vector<std::vector<double>> pH,
                                                                     std::vector<double> mean1,
                                                                     std::vector<double> mean2,
                                                                     std::vector<double> variables1 = std::vector<double>(),
                                                                     std::vector<double> variables2 = std::vector<double>(),
                                                                     std::shared_ptr<InterpolateLikelihood> i = NULL,
								     int iType = 0, bool wSwitch = false,
                                                                     std::string name = "") {

  std::cout << (variables1.size() == 0 || variables2.size() == 0) << std::endl;
  std::cout << bool(wSpline == NULL) << std::endl;
  if ((variables1.empty() || variables2.empty()) && wSpline == NULL) {
    std::cout << "NO SPLINE NO 5 PAR" << std::endl;
    return getLogLikelihoodWithNuisanceMinuit(pH, mean1, mean2);
  } else if ((variables1.empty() || variables2.empty()) && wSpline != NULL) {
    std::cout << "SPLINE" << std::endl;
    double sig1 = .00315; //These are taken to be constants
    double sig2 = .002851;

    std::vector<double> ll;
    for (int eventCount = 0; eventCount < pH.size(); eventCount++) {
      double ll1, ll2 = 0;
      for (int pixelCount = 0; pixelCount < pH[eventCount].size(); pixelCount++) {
        if (pH[eventCount][pixelCount] >= 1) {
          ll1 = ll1 + -2 * calcLikelihoodForEventNuisance(pH[eventCount][pixelCount], mean1[pixelCount], wSpline->Eval(mean1[pixelCount]), sig1);
          ll2 = ll2 + -2 * calcLikelihoodForEventNuisance(pH[eventCount][pixelCount], mean2[pixelCount], wSpline->Eval(mean2[pixelCount]), sig2);
        } else {
          double w = -pow(sig1, 2);
          ll1 = ll1 + -2 * calcLikelihoodForEventNuisance(pH[eventCount][pixelCount], mean1[pixelCount], w, sig1);
          w = -pow(sig2, 2);
          ll2 = ll2 + -2 * calcLikelihoodForEventNuisance(pH[eventCount][pixelCount], mean2[pixelCount], w, sig2);
        }
      }
      ll.push_back(ll2 - ll1);
    }
    return ll;
  } else {
    std::cout << "5 variable" << std::endl;
    return getLogLikelihoodForVariableAndWparameter(pH, mean1, mean2, variables1, variables2, i, iType, wSwitch, name);
  }
}

double Detector::calcLikelihoodForParametersNuisanceForParticle(TMinuit *gMinuit, std::vector<double> variables, double sig) {
  gMinuit->SetPrintLevel(-1);
  gMinuit->SetFCN(FCN3);
  gMinuit->SetObjectFit((TObject *)this);

  double arglist[10];
  int ierflg = 0;

  arglist[0] = 2; //Migrad and Hesse (correctly account for correlation)
  gMinuit->mnexcm("SET STR", arglist, 1, ierflg);

  arglist[0] = 0;
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  static Double_t vstart[5] = {variables[0], variables[1], variables[2], variables[3], variables[4]};
  static Double_t step[5] = {0.0001, 0.0001, 0.00001, 1E-9, 0.09};

  gMinuit->mnparm(0, "x", vstart[0], step[0], 0, 0, ierflg);
  gMinuit->mnparm(1, "y", vstart[1], step[1], 0, 0, ierflg);
  gMinuit->mnparm(2, "t", vstart[2], step[2], 0, 0, ierflg);
  gMinuit->mnparm(3, "p", vstart[3], step[3], 0, 0, ierflg);
  gMinuit->mnparm(4, "bm", vstart[4], step[4], 0, 0, ierflg);

  gMinuit->mnparm(5, "errX", 5E-7, 0, 0, 0, ierflg);
  gMinuit->mnparm(6, "erry", 5E-7, 0, 0, 0, ierflg);
  gMinuit->mnparm(7, "errt", 0.5 / 1000.0, 0, 0, 0, ierflg);
  gMinuit->mnparm(8, "errp", 7E-7, 0, 0, 0, ierflg);
  gMinuit->mnparm(9, "errBM", .13 * variables[4], 0, 0, 0, ierflg);

  gMinuit->FixParameter(5);
  gMinuit->FixParameter(6);
  gMinuit->FixParameter(7);
  gMinuit->FixParameter(8);
  gMinuit->FixParameter(9);

  setVariables(variables);
  gMinuit->mnparm(10, "sig", sig, 0, 0, 0, ierflg);
  gMinuit->FixParameter(10);
  arglist[0] = 500;
  gMinuit->mnexcm("MIGRAD", arglist, 1, ierflg);
  gMinuit->mnexcm("HESSE", arglist, 1, ierflg);
  //gMinuit->Command("SCAn");
  //arglist[0] = 0;
  // get the graph and write it to the root file:
  //for (int i = 0; i < variables.size(); i++) {
  //arglist[0] = i;
  //arglist[1] = 90;
  //gMinuit->mnexcm("SCAN", arglist, 1, ierflg);
  //TGraph *gr = (TGraph *)gMinuit->GetPlot();
  //gr->Write();
  //}
  double ll, edm, errdef=0;
  int nvpar, nparx, icstat=0;
  gMinuit->mnstat(ll, edm, errdef, nvpar, nparx, icstat);

  double x, y, t, p, b = 0;
  double xErr, yErr, tErr, pErr, bErr = 0;

  gMinuit->GetParameter(0, x, xErr);
  gMinuit->GetParameter(1, y, yErr);
  gMinuit->GetParameter(2, t, tErr);
  gMinuit->GetParameter(3, p, pErr);
  gMinuit->GetParameter(4, b, bErr);

  std::cout << x - variables[0] << "   ";
  std::cout << y - variables[1] << "   ";
  std::cout << t - variables[2] << "   ";
  std::cout << p - variables[3] << "   ";
  std::cout << b - variables[4] << "   ";
  std::cout << ll << std::endl;
  return ll;
}

std::vector<double> Detector::getLogLikelihoodWithNuisanceMinuit(std::vector<std::vector<double>> pH,
                                                                 std::vector<double> mean1,
                                                                 std::vector<double> mean2) {

  std::clock_t begin = std::clock();
  double sig1 = .00315; //These are taken to be constants
  double sig2 = .002851;

  std::vector<double> ll_pH;
  //root file that shall store the graph :
  //TFile *f = new TFile("MinuitGraph.root", "RECREATE");
  for (int eventCount = 0; eventCount < pH.size(); eventCount++) {
    double ll_particle_1 = calcLikelihoodForParticleNuisanceMinuit(pH[eventCount], mean1, sig1);
    double ll_particle_2 = calcLikelihoodForParticleNuisanceMinuit(pH[eventCount], mean2, sig2);
    ll_pH.push_back(ll_particle_1 - ll_particle_2);
  }
  std::clock_t end = std::clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout << "MINUIT: " << elapsed_secs << std::endl;
  return ll_pH;
}

double Detector::calcLikelihoodWithParameterNuisanceMinuit(TMinuit *gMinuit, double ll, std::vector<double> variables) {
  double arglist[10];
  int ierflg = 0;

  arglist[0] = 0; //Migrad and Hesse (correctly account for correlation)
  gMinuit->mnexcm("SET STR", arglist, 1, ierflg);

  arglist[0] = 0;
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  static Double_t vstart[5] = {0.097, 0.097, 0.0, 0.0, 1};
  static Double_t step[5] = {0.0001, 0.0001, 0.00001, 1E-9, 1};
  gMinuit->mnparm(0, "x", vstart[0], step[0], 0, 0, ierflg);
  gMinuit->mnparm(1, "y", vstart[1], step[1], 0, 0, ierflg);
  gMinuit->mnparm(2, "t", vstart[2], step[2], 0, 0, ierflg);
  gMinuit->mnparm(3, "p", vstart[3], step[3], 0, 0, ierflg);
  gMinuit->mnparm(4, "bm", vstart[4], step[4], 0, 0, ierflg);

  gMinuit->mnparm(5, "errX", 5E-7, 0, 0, 0, ierflg);
  gMinuit->mnparm(6, "erry", 5E-7, 0, 0, 0, ierflg);
  gMinuit->mnparm(7, "errt", 0.5 / 1000.0, 0, 0, 0, ierflg);
  gMinuit->mnparm(8, "errp", 5E-7, 0, 0, 0, ierflg);
  gMinuit->mnparm(9, "errBM", .13 * variables[4], 0, 0, 0, ierflg);

  gMinuit->FixParameter(5);
  gMinuit->FixParameter(6);
  gMinuit->FixParameter(7);
  gMinuit->FixParameter(8);
  gMinuit->FixParameter(9);

  gMinuit->mnparm(10, "ll", ll, 0, 0, 0, ierflg);
  gMinuit->FixParameter(10);

  arglist[0] = 500;
  gMinuit->mnexcm("MIGRAD", arglist, 1, ierflg);
  gMinuit->mnexcm("HESSE", arglist, 1, ierflg);

  double ll_new, edm, errdef;
  int nvpar, nparx, icstat;
  gMinuit->mnstat(ll_new, edm, errdef, nvpar, nparx, icstat);
  return ll_new;
}

void graphHistogramMinuit(int i, TMinuit *gMinuit) {
  double arglist[10];
  int ierflg = 0;

  arglist[0] = i;
  arglist[1] = 90;
  gMinuit->mnexcm("SCAN", arglist, 1, ierflg);
  TGraph *gr = (TGraph *)gMinuit->GetPlot();
  gr->Write();
}

std::vector<double> Detector::getLogLikelihoodForParticleWithNuisanceMinuit(std::vector<std::vector<double>> pH,
                                                                            std::vector<double> mean,
                                                                            double sig) {
  std::clock_t begin = std::clock();
  std::vector<double> ll_pH;

  for (int eventCount = 0; eventCount < pH.size(); eventCount++) {
    double ll_particle_1 = calcLikelihoodForParticleNuisanceMinuit(pH[eventCount], mean, sig);
    ll_pH.push_back(ll_particle_1);
  }
  return ll_pH;
}

double Detector::calcLikelihoodForParticleNuisanceMinuit(std::vector<double> pH,
                                                         std::vector<double> mean,
                                                         double sig) {
  double ln = 0;
  //if (wSpline != NULL) {
  //for (int i = 0; i < mean.size(); i++) {
  //ln = ln + -2 * calcLikelihoodForEventNuisance(pH[i], mean[i], wSpline->Eval(mean[i]), sig);
  //}
  //return ln;
  //}
  TMinuit *gMinuit = new TMinuit(1);
  //gMinuit->SetPrintLevel(-1);
  gMinuit->SetFCN(FCN);
  gMinuit->SetObjectFit((TObject *)this);
  double arglist[10];
  int ierflg = 0;

  arglist[0] = 0; //Migrad and Hesse (correctly account for correlation)
  gMinuit->mnexcm("SET STR", arglist, 1, ierflg);

  arglist[0] = 0;
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  //TH1D *h = new TH1D("wParameter", "wParameter", 1000, -1, 1);
  for (int i = 0; i < mean.size(); i++) {
    setMeanSingleEvent(mean[i]);
    setPixelHitSingleEvent(pH[i]);

    static Double_t vstart[2] = {0.0005005, sig};
    static Double_t step[2] = {0.000001, 0};
    gMinuit->mnparm(0, "w", vstart[0], step[0], 0, 0, ierflg);
    gMinuit->mnparm(1, "sig", vstart[1], step[1], 0, 0, ierflg);

    gMinuit->FixParameter(1);

    arglist[0] = 500;
    gMinuit->mnexcm("MIGRAD", arglist, 1, ierflg);
    //gMinuit->mnexcm("HESSE", arglist, 1, ierflg);
    //gMinuit->Command("SCAn");
    //arglist[0] = 0;
    // get the graph and write it to the root file:
    //arglist[0] = 0;
    //arglist[1] = 90;
    //gMinuit->mnexcm("SCAN", arglist, 1, ierflg);
    //TGraph *gr = (TGraph *)gMinuit->GetPlot();
    //gr->Write();

    double ln0, edm, errdef;
    int nvpar, nparx, icstat;

    double w, wErr = 0;
    gMinuit->GetParameter(0, w, wErr);
    //std::cout << i << "    " << w << "   " << std::endl;
    //if (pH[i] >= 1) {
    //std::cout << "TEST" << std::endl;
    quickSaveWParameter(mean[i], sig, w);
    //}

    //h->Fill(w);

    //std::cout << i << "   ";
    //std::cout << mean[i] << "   ";
    //std::cout << sig << "   ";
    //std::cout << w << "    ";
    //std::cout << ln0 << std::endl;

    gMinuit->mnstat(ln0, edm, errdef, nvpar, nparx, icstat);
    ln = ln + ln0;
  }
  //h->Print("wParameterHistogram");
  return ln;
}

void Detector::quickSaveWParameter(double m, double s, double w) {
  std::string pathName = "./Textfiles/";
  mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  pathName += "wParameter.txt";
  FILE *f = fopen(pathName.c_str(), "at");
  assert(f);
  fprintf(f, "%g ", m);
  fprintf(f, "%g ", s);
  fprintf(f, "%g ", w);
  fwrite("\n", sizeof(char), 1, f);
  fclose(f);
  return;
}

void Detector::quickReadWParameter() {
  std::ifstream infile("wParameter_sorted_removeDuplicates.txt");
  if (!infile.is_open()) {
    std::cout << "No File " << std::endl;
    return;
  }
  std::vector<double> x;
  std::vector<double> y;
  //x.reserve(100);
  //y.reserve(100);
  int i = 0;
  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty()) // be careful: an empty line might be read
    {
      continue;
    }
    std::istringstream iss(line);
    double mean;
    double sig;
    double w;
    iss >> mean >> w;
    x.push_back(mean);
    y.push_back(w);
    i++;
  }
  
  std::shared_ptr<TGraphErrors> gr(std::make_shared<TGraphErrors>(x.size(), &(x[0]), &(y[0])));
  gr->Sort();
  std::shared_ptr<TSpline3> grs(std::make_shared<TSpline3>("grs", gr.get()));
  wSpline = grs;

  //double steps = 99;
  //double max = wSpline->GetXmax();
  //double min = wSpline->GetXmin();
  //double stepSize = (max - min) / (steps - 1);
  ////std::cout << steps << "  " << max << "  " << min << "  " << stepSize << std::endl;

  double stepSize =  double(1.126) / double(9.0);

  for (int i = 0; i < 10; i++) {
    std::cout << stepSize * i << "   ";
    std::cout << wSpline->Eval(stepSize * i) << "    ";
    //double w = EvaluateSpline(wSpline.get(), steps, max, min, stepSize, stepSize * i);
    std::cout << exp(-1.05 * log(-2.965 + 2.975 * exp(0.8316 * stepSize * i)) - 10.64) << std::endl;
  }

  return;
}

double Detector::calcLikelihoodForPixelNuisance(std::vector<double> pH,
                                                std::vector<double> mean,
                                                double w,
                                                double sig) {
  double ll_pixel = 0;
  for (int pixelCount = 0; pixelCount < pH.size(); pixelCount++) {
    ll_pixel = ll_pixel - 2 * calcLikelihoodForEventNuisance(pH[pixelCount], mean.at(pixelCount), w, sig);
  }
  return ll_pixel;
}

double Detector::calcLikelihoodForEventNuisance(double pH_event, double mean_pixel, double w, double sig) {
  double threshold = 1;
  mean_pixel = mean_pixel + pmt->getDarkNoise();
  //std::cout << pH_event << "    ";
  //std::cout << mean_pixel << "    ";
  //std::cout << w << "    ";
  //std::cout << sig << std::endl;
  // std::cout << mean_pixel << std::endl;
  if (pH_event >= threshold) {
    //if (mean_pixel == w) {
    //w = w + mean_pixel / 10; //Penalty term if equal
    //}
    //double x = mean_pixel + w;
    return log(1 - exp(-mean_pixel - w)) - (w * w) / (2 * sig * sig);
  } else {
    //w = -pow(sig, 2);
    return -1.0 * (mean_pixel + w) - (w * w) / (2 * sig * sig);
  }
}

void minuitFunctionWVariables(Int_t &npar, Double_t *gin, Double_t &f,
                              Double_t *u, Int_t flag) {
  f = 0;

  //double vMin[5] = {0, 0, -0.4, 0, 0};
  //double vMax[5] = {2 * 0.097, 2 * 0.097, 0.4, 2 * TMath::Pi(), 1};
  
  Detector *d = (Detector *)gMinuit->GetObjectFit();
  std::vector<double> variables = d->getVariables();
  std::vector<double> variablesTrue;
  variablesTrue.reserve(5);
  for (int v = 0; v < variables.size(); v++) {
    if (u[v] > .9999) {
      u[v] = .9999;
      f = 50;
    }
    variablesTrue.push_back(u[v]);
    //if (variablesTrue[v] < vMin[v] || variablesTrue[v] > vMax[v]) {
    //f = 1E4;
    //}
  }

  std::shared_ptr<InterpolateLikelihood> i = d->getInterpolator();
  int iType = d->getInterpolatorType();
  bool w_on = d->getWSwitch();
  std::vector<double> pH = d->getPixelHit();
  std::vector< double> mean;
  double sig = u[10];
  if (i == NULL) {
     std::vector<double> m = d->getMean();
     mean= std::vector<double>(m.begin(),m.end());
  } else {
    std::vector<bool> variableSwitch;
    variableSwitch.reserve(5);
    variableSwitch.push_back(false); 
    variableSwitch.push_back(false);
    variableSwitch.push_back(false);
    variableSwitch.push_back(false);
    variableSwitch.push_back(true); 
    mean = i->interpolationSwitch(variablesTrue, iType, variableSwitch,
                                  d->getPMT()->getWidthPixelCount(), d->getTotalDistance(),
                                  d->getMaxWidth(), d->getMaxLength(),
                                  d->getPMT()->getWidthPixel(), d->getPMT()->getLengthPixel());
  }
  std::shared_ptr<TSpline3> wSpline = d->getWSpline();
  if (wSpline == NULL) {
    f = d->calcLikelihoodForParticleNuisanceMinuit(pH, 
					std::vector<double>(mean.begin(),mean.end()),
					 sig);
  } else {
    double ll = 0;
    double meanSum = 0;
    double delta = 0;
    double wSum = 0;
    for (int pixelCount = 0; pixelCount < pH.size(); pixelCount++) {
      delta += fabs(mean[pixelCount] - pH[pixelCount]);
      if (pH[pixelCount] >= 1) {
        //double steps = 100;
        //double max = wSpline->GetXmax();
        //double min = wSpline->GetXmin();
        //double stepSize = (max - min) / (steps - 1);
        //std::cout << steps << "  " << max << "  " << min << "  " << stepSize << std::endl;
        //double w = d->EvaluateSpline(wSpline.get(), steps, max, min, stepSize, mean[pixelCount]);
        double w = exp(-1.05 * log(-2.965 + 2.975 * exp(0.8316 * mean[pixelCount])) - 10.64);
        if (!w_on) {w = 0;}
        ll += -2 * d->calcLikelihoodForEventNuisance(pH[pixelCount], mean[pixelCount], w, sig);
        wSum += w;
        meanSum += mean[pixelCount];
      } else {
        double w = -pow(sig, 2);
        if(!w_on){w = 0;}
        ll += -2 * d->calcLikelihoodForEventNuisance(pH[pixelCount], mean[pixelCount], w, sig);
        wSum += w;
        meanSum += mean[pixelCount];
      }
    }
    //std::cout << "dSum: " << delta << "   ";
    //std::cout << "meanSum: " << meanSum << "   ";
    //std::cout << "wSum:" << wSum << "    ";
    f += ll;
  }
  
  gMinuit = d->getMinuit();

  for (int v = 0; v < variables.size(); v++) {
    f = f + ((-variables[v] + u[v]) * (-variables[v] + u[v])) / (u[v + variables.size()] * u[v + variables.size()]);
    //std::cout << u[v] << "+/-" << u[v + variables.size()] << "    ";
  }
  //std::cout << f << std::endl;
}

TMinuit *Detector::getWMinuit() {
  TMinuit *gMinuit = new TMinuit(1);
  //  gMinuit->SetPrintLevel(-1);
  gMinuit->SetFCN(minuitFunctionWVariables);
  gMinuit->SetObjectFit((TObject *)this);
  return gMinuit;
}

std::vector<double> Detector::getLogLikelihoodForVariableAndWparameter(std::vector<std::vector<double>> pH,
                                                                       std::vector<double> mean1,
                                                                       std::vector<double> mean2,
                                                                       std::vector<double> variables1 = std::vector<double>(),
                                                                       std::vector<double> variables2 = std::vector<double>(),
                                                                       std::shared_ptr<InterpolateLikelihood> i = NULL,
								       int iType = 0,
								       bool wSwitch =false,
                                                                       std::string name = "") {
  std::clock_t begin = std::clock();
  double sig1 = .00315; //These are taken to be constants
  double sig2 = .002851;

  std::vector<double> ll_pH;

  TMinuit *gMinuit = new TMinuit(5);
  gMinuit->SetPrintLevel(-1);
  gMinuit->SetFCN(minuitFunctionWVariables);
  gMinuit->SetObjectFit((TObject *)this);

  this->setInterpolator(i,iType);
  this->setWSwitch(wSwitch);

  for (int eventCount = 0; eventCount < pH.size(); eventCount++) {
    std::cout << eventCount << std::endl;
    double ll_particle_1 = getMinuitVariableAndWParameterMinimization(gMinuit, pH[eventCount], mean1, sig1, variables1, name + "1");
    double ll_particle_2 = getMinuitVariableAndWParameterMinimization(gMinuit, pH[eventCount], mean2, sig2, variables2, name + "2");
    ll_pH.push_back(ll_particle_1 - ll_particle_2);
  }

  std::clock_t end = std::clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout << "  TIME: " << elapsed_secs << "\n";
  return ll_pH;
}

std::vector<double> Detector::getLogLikelihoodForVariableAndWparameter(std::vector<std::vector<double>> pH,
                                                                       std::vector<double> mean1,
                                                                       std::vector<double> mean2,
                                                                       std::shared_ptr<Beam> beam1 = NULL,
                                                                       std::shared_ptr<Beam> beam2 = NULL,
                                                                       std::shared_ptr<InterpolateLikelihood> i = NULL,
 								       int iType = 0,
								       bool wSwitch =false,
								       std::string name = "") {
  std::clock_t begin = std::clock();
  double sig1 = .00315; //These are taken to be constants
  double sig2 = .002851;

  std::vector<double> ll_pH;

  TMinuit *gMinuit = new TMinuit(5);
  gMinuit->SetPrintLevel(-1);
  gMinuit->SetFCN(minuitFunctionWVariables);
  gMinuit->SetObjectFit((TObject *)this);

  this->setInterpolator(i,iType);

  for (int eventCount = 0; eventCount < pH.size(); eventCount++) {
    std::cout << eventCount << std::endl;
    double ll_particle_1 = getMinuitVariableAndWParameterMinimization(gMinuit, pH[eventCount], mean1, sig1, beam1, name + "1");
    double ll_particle_2 = getMinuitVariableAndWParameterMinimization(gMinuit, pH[eventCount], mean2, sig2, beam2, name + "2");
    ll_pH.push_back(ll_particle_1 - ll_particle_2);
  }

  std::clock_t end = std::clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout << " TIME: " << elapsed_secs << "\n";
  return ll_pH;
}

void Detector::getLogLikelihoodForVariableAndWparameterSingleParticle(std::vector<std::vector<double>> pH,
                                                                      std::vector<double> mean,
                                                                      std::vector<std::vector<double>> &v,
                                                                      std::vector<std::vector<double>> &v_err,
                                                                      std::vector<std::vector<double>> &m_war,
                                                                      std::vector<double> &v_ll,
                                                                      std::shared_ptr<Beam> beam = NULL,
                                                                      std::shared_ptr<InterpolateLikelihood> i = NULL,
								      int interpolatorType = 0,
								      bool w_Switch = false,
                                                                      std::string name = "") {
  double sig = .00315; //These are taken to be constants

  std::vector<double> ll_pH;

  std::shared_ptr<TMinuit> gMinuit(std::make_shared<TMinuit>(5));
  gMinuit->SetPrintLevel(-1);
  gMinuit->SetFCN(minuitFunctionWVariables);
  gMinuit->SetObjectFit((TObject *)this);
  this->setInterpolator(i,interpolatorType);
  this->setWSwitch(w_Switch);
  std::vector<double> v_perEvent;
  std::vector<double> v_perEvent_err;
  v_perEvent.reserve(5);
  v_perEvent_err.reserve(5);
  std::vector<double> m_war_perEvent;
  m_war_perEvent.reserve(5);
  for (int eventCount = 0; eventCount < pH.size(); eventCount++) {
    std::cout << eventCount << "  ";
    std::clock_t begin = std::clock();
    v_ll.push_back(getMinuitVariableAndWParameterMinimizationReturnParameter(gMinuit, pH[eventCount], mean, sig,
                                                                             v_perEvent, v_perEvent_err, m_war_perEvent,
                                                                             beam, name+std::to_string(eventCount)));
    v.push_back(v_perEvent);
    v_err.push_back(v_perEvent_err);
    m_war.push_back(m_war_perEvent);
    std::clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << " TIME: " << elapsed_secs << "\n";
    v_perEvent.clear();
    v_perEvent_err.clear();
    m_war_perEvent.clear();
  }
  return;
}

double Detector::EvaluateSpline(TSpline3 *spline, int nSteps, double max, double min, double stepSize, double pos) {

  double tmp;
  double x, y, b, c, d, num;
  int l;
  if (spline == NULL)
    return 1.0;
  //    std::cout << "pos = " << pos << std::endl;
  //    std::cout << "max = " << max << std::endl;
  //    std::cout << "min = " << min << std::endl;

  if (pos < min) {
    spline->GetCoeff(0, x, y, b, c, d);
    num = pos - x;
    tmp = (y + num * b + num * num * c + num * num * num * d);
  } else if (pos > max) {
    spline->GetCoeff(nSteps - 1, x, y, b, c, d);
    num = pos - x;
    tmp = (y + num * b + num * num * c + num * num * num * d);
  } else {
    l = int((pos - min) / stepSize) + 1;
    spline->GetCoeff(l, x, y, b, c, d);
    num = pos - x;
    if (num < 0) {
      l = l - 1;
      spline->GetCoeff(l, x, y, b, c, d);
      num = pos - x;
    }
    tmp = (y + num * b + num * num * c + num * num * num * d);
    //        if(fabs(tmp - spline->Eval(pos)) > 0.00001){
    //            std::cout << l << std::endl;
    //            std::cout << num << std::endl;
    //            std::cout << pos << std::endl;
    //            std::cout << stepSize << std::endl;
    //        }
  }
  if (tmp == tmp)
    return tmp;
  else
    std::cout << "Spline calculation error" << std::endl;
  return 0.0;
}

double Detector::getMinuitVariableAndWParameterMinimization(TMinuit *gMinuit,
                                                            std::vector<double> pH,
                                                            std::vector<double> mean,
                                                            double sig,
                                                            std::vector<double> v,
                                                            std::string name = "") {
  gMinuit->SetFCN(minuitFunctionWVariables);
  gMinuit->SetObjectFit((TObject *)this);
  double arglist[10];
  int ierflg = 0;

  arglist[0] = 0; //Migrad and Hesse (correctly account for correlation)
  gMinuit->mnexcm("SET STR", arglist, 1, ierflg);

  //arglist[0] = 1; //Has to be set at -2log likelihood
  //gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  //arglist[0] = 0.01; //Set Precision of Minuit to 1E-4
  //gMinuit->mnexcm("SET EPS", arglist, 1, ierflg);

  //std::vector<double> v;
  //v.reserve(5);
  //v.push_back(beam->getPosition()->getX());
  //v.push_back(beam->getPosition()->getY());
  //v.push_back(beam->getDirection()->getTheta());
  //v.push_back(beam->getDirection()->getPhi());
  //v.push_back(beam->getBeta());

  setMean(mean);
  setPixelHit(pH);
  setVariables(v);
  double vSig[5] = {5E-7, 5E-7, 0.5 / 1000, 7E-7, 0.13 * v[4]};
  //double particle = beam->getParticleMass();
  //vSig[4] = vSig[4] / sqrt(pow(vSig[4], 2) + pow(particle, 2)); //beta unceratainty
  //static Double_t step[5] = {, , 1, 1, 1};

  gMinuit->mnparm(0, "x", v[0], 0.01, 0, 0, ierflg);  //0, 2 * 0.097
  gMinuit->mnparm(1, "y", v[1], 0.01, 0, 0, ierflg);  //0, 2 * 0.097
  gMinuit->mnparm(2, "t", v[2], 0.01, 0, 0, ierflg);  //-0.4, 0.4,
  gMinuit->mnparm(3, "p", v[3], 0.01, 0, 0, ierflg);  //0, 2 * TMath::Pi()
  gMinuit->mnparm(4, "bm", v[4], 0.01, 0, 0, ierflg); // 1 / 1.012, 1.001, ierflg); //0, 1

  gMinuit->mnparm(5, "errX", vSig[0], 0, 0, 0, ierflg);
  gMinuit->mnparm(6, "errY", vSig[1], 0, 0, 0, ierflg);
  gMinuit->mnparm(7, "errT", vSig[2], 0, 0, 0, ierflg);
  gMinuit->mnparm(8, "errP", vSig[3], 0, 0, 0, ierflg);
  gMinuit->mnparm(9, "errBM", vSig[4], 0, 0, 0, ierflg);

  //gMinuit->FixParameter(1);
  //gMinuit->FixParameter(2);
  //gMinuit->FixParameter(3);
  //  gMinuit->FixParameter(4);

  gMinuit->FixParameter(5);
  gMinuit->FixParameter(6);
  gMinuit->FixParameter(7);
  gMinuit->FixParameter(8);
  gMinuit->FixParameter(9);

  gMinuit->mnparm(10, "sig", sig, 0, 0, 0, ierflg);
  gMinuit->FixParameter(10);

  double ln, edm, errdef;
  int nvpar, nparx, icstatMigrad, icstatHesse;

  arglist[0] = 500;
  gMinuit->mnexcm("MIGRAD", arglist, 1, ierflg);
  gMinuit->mnstat(ln, edm, errdef, nvpar, nparx, icstatMigrad);
  gMinuit->mnexcm("HESSE", arglist, 1, ierflg);
  gMinuit->mnstat(ln, edm, errdef, nvpar, nparx, icstatHesse);

  //gMinuit->mnparm(4, "bm", v[4], 0.01, 0, 0, ierflg);
  //gMinuit->mnexcm("MIGRAD", arglist, 1, ierflg);
  //gMinuit->mnstat(ln, edm, errdef, nvpar, nparx, icstatMigrad);
  //gMinuit->mnexcm("HESSE", arglist, 1, ierflg);
  //gMinuit->mnstat(ln, edm, errdef, nvpar, nparx, icstatHesse);

  double x, y, t, p, b = 0;
  double xErr, yErr, tErr, pErr, bErr = 0;

  gMinuit->GetParameter(0, x, xErr);
  gMinuit->GetParameter(1, y, yErr);
  gMinuit->GetParameter(2, t, tErr);
  gMinuit->GetParameter(3, p, pErr);
  gMinuit->GetParameter(4, b, bErr);

  double minuiVar[5] = {x, y, t, p, b};
  double minuiVarErr[5] = {xErr, yErr, tErr, pErr, bErr};

  //root file that shall store the graph :
  name += "_MinuitGraph.root";
  TFile *f = new TFile(name.c_str(), "recreate");

  for (int i = 1; i < v.size() + 1; i++) {
    //get the graph and write it to the root file:
    //int i = 0;
    arglist[0] = i;
    arglist[1] = 90;
    arglist[2] = minuiVar[i] - 3 * minuiVarErr[i];
    arglist[3] = minuiVar[i] + 3 * minuiVarErr[i];
    gMinuit->mnexcm("SCAN", arglist, 1, ierflg);
    TGraph *gr = (TGraph *)gMinuit->GetPlot();
    gr->SetTitle(std::to_string(i).c_str());
    gr->Write();
  }
  f->Close();
  //std::cout << "DX: " << x << "+/-" << xErr << "   ";
  //std::cout << "DY: " << y << "+/-" << yErr << "   ";
  std::cout << "DT: " << t << "+/-" << tErr << "   ";
  //std::cout << "DP: " << p << "+/-" << pErr << "   ";
  std::cout << "DB: " << b << "+/-" << bErr << "   ";
  std::cout << "LL: " << ln << "   ";
  std::cout << "SM: " << icstatMigrad << "   ";
  std::cout << "SH: " << icstatHesse << "   ";
  return ln;
}

double Detector::getMinuitVariableAndWParameterMinimization(TMinuit *gMinuit,
                                                            std::vector<double> pH,
                                                            std::vector<double> mean,
                                                            double sig,
                                                            std::shared_ptr<Beam> beam = NULL,
                                                            std::string name = "") {
  gMinuit->SetFCN(minuitFunctionWVariables);
  gMinuit->SetObjectFit((TObject *)this);
  double arglist[10];
  int ierflg = 0;

  //arglist[0] = 1;
  //gMinuit->mnexcm("SET WAR", arglist, 1, ierflg);

  arglist[0] = 0; //Migrad
  gMinuit->mnexcm("SET STR", arglist, 1, ierflg);

  //arglist[0] = 1; //Has to be set at -2log likelihood
  //gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  //arglist[0] = 0.01; //Set Precision of Minuit to 1E-4
  //gMinuit->mnexcm("SET EPS", arglist, 1, ierflg);

  std::vector<double> v;
  v.reserve(5);
  v.push_back(beam->getPosition()->getX());
  v.push_back(beam->getPosition()->getY());
  v.push_back(beam->getDirection()->getTheta());
  v.push_back(beam->getDirection()->getPhi());
  v.push_back(beam->getBeta());

  setMean(mean);
  setPixelHit(pH);
  setVariables(v);
  double vSig[5] = {5E-7, 5E-7, 0.5 / 1000, 7E-7, 0.13 * beam->getMomentum()};
  double particle = beam->getParticleMass();
  vSig[4] = vSig[4] * pow(particle, 2) / pow(particle * particle + beam->getMomentum() * beam->getMomentum(), 1.5);
  //static Double_t step[5] = {, , 1, 1, 1};

  gMinuit->mnparm(0, "x", v[0], 0.01, 0, 0, ierflg);                //0, 2 * 0.097
  gMinuit->mnparm(1, "y", v[1], 0.01, 0, 0, ierflg);                //0, 2 * 0.097
  gMinuit->mnparm(2, "t", v[2], 0.01, 0, 0, ierflg);                //-0.4, 0.4,
  gMinuit->mnparm(3, "p", v[3], 0.01, 0, 0, ierflg);                //0, 2 * TMath::Pi()
  gMinuit->mnparm(4, "bm", v[4], 0.01, 1 / 1.012, 0.99992, ierflg); // 1 / 1.012, 1.001, ierflg); //0, 1

  gMinuit->mnparm(5, "errX", vSig[0], 0, 0, 0, ierflg);
  gMinuit->mnparm(6, "errY", vSig[1], 0, 0, 0, ierflg);
  gMinuit->mnparm(7, "errT", vSig[2], 0, 0, 0, ierflg);
  gMinuit->mnparm(8, "errP", vSig[3], 0, 0, 0, ierflg);
  gMinuit->mnparm(9, "errBM", vSig[4], 0, 0, 0, ierflg);

  //gMinuit->FixParameter(1);
  //gMinuit->FixParameter(2);
  //gMinuit->FixParameter(3);
  //  gMinuit->FixParameter(4);

  gMinuit->FixParameter(5);
  gMinuit->FixParameter(6);
  gMinuit->FixParameter(7);
  gMinuit->FixParameter(8);
  gMinuit->FixParameter(9);

  gMinuit->mnparm(10, "sig", sig, 0, 0, 0, ierflg);
  gMinuit->FixParameter(10);

  double ln, edm, errdef;
  int nvpar, nparx, icstatMigrad, icstatHesse;

  int err_M, err_H;
  arglist[0] = 500;
  gMinuit->mnexcm("MIGRAD", arglist, 1, ierflg);
  err_M = ierflg;
  gMinuit->mnstat(ln, edm, errdef, nvpar, nparx, icstatMigrad);
  gMinuit->mnexcm("HESSE", arglist, 1, ierflg);
  err_H = ierflg;
  gMinuit->mnstat(ln, edm, errdef, nvpar, nparx, icstatHesse);

  //gMinuit->mnparm(4, "bm", v[4], 0.01, 0, 0, ierflg);
  //gMinuit->mnexcm("MIGRAD", arglist, 1, ierflg);
  //gMinuit->mnstat(ln, edm, errdef, nvpar, nparx, icstatMigrad);
  //gMinuit->mnexcm("HESSE", arglist, 1, ierflg);
  //gMinuit->mnstat(ln, edm, errdef, nvpar, nparx, icstatHesse);

  double x,y,t,p,b = 0;
  double xErr, yErr, tErr, pErr, bErr= 0;

  gMinuit->GetParameter(0, x, xErr);
  gMinuit->GetParameter(1, y, yErr);
  gMinuit->GetParameter(2, t, tErr);
  gMinuit->GetParameter(3, p, pErr);
  gMinuit->GetParameter(4, b, bErr);

  double minuiVar[5] = {x, y, t, p, b};
  double minuiVarErr[5] = {xErr, yErr, tErr, pErr, bErr};

  //root file that shall store the graph :
  name += "_MinuitGraph.root";
  TFile *f = new TFile(name.c_str(), "recreate");

  for (int i = 1; i < v.size() + 1; i++) {
    //get the graph and write it to the root file:
    //int i = 0;
    arglist[0] = i;
    arglist[1] = 90;
    arglist[2] = minuiVar[i] - 3 * minuiVarErr[i];
    arglist[3] = minuiVar[i] + 3 * minuiVarErr[i];
    gMinuit->mnexcm("SCAN", arglist, 1, ierflg);
    TGraph *gr = (TGraph *)gMinuit->GetPlot();
    gr->SetTitle(std::to_string(i).c_str());
    gr->Write();
  }
  f->Close();
  //std::cout << "DX: " << x << "+/-" << xErr << "   ";
  //std::cout << "DY: " << y << "+/-" << yErr << "   ";
  std::cout << "DT: " << t << "+/-" << tErr << "   ";
  //std::cout << "DP: " << p << "+/-" << pErr << "   ";
  std::cout << "DB: " << b << "+/-" << bErr << "   ";
  std::cout << "LL: " << ln << "   ";
  std::cout << "WM: " << err_M << "   ";
  std::cout << "WH: " << err_H << "   ";
  std::cout << "SM: " << icstatMigrad << "   ";
  std::cout << "SH: " << icstatHesse << "   ";
  return ln;
}

double Detector::getMinuitVariableAndWParameterMinimizationReturnParameter(std::shared_ptr<TMinuit> gMinuit,
                                                                           std::vector<double> pH,
                                                                           std::vector<double> mean,
                                                                           double sig,
                                                                           std::vector<double> &minuiVar,
                                                                           std::vector<double> &minuiVarErr,
                                                                           std::vector<double> &minuiWarning,
                                                                           std::shared_ptr<Beam> beam = NULL,
                                                                           std::string name = "") {
  gMinuit->SetFCN(minuitFunctionWVariables);
  gMinuit->SetObjectFit((TObject *)this);
  double arglist[10];
  int ierflg = 0;

  //arglist[0] = 1;
  //gMinuit->mnexcm("SET WAR", arglist, 1, ierflg);

  arglist[0] = 0; //Migrad
  gMinuit->mnexcm("SET STR", arglist, 1, ierflg);

  //arglist[0] = 1; //Has to be set at -2log likelihood
  //gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  //arglist[0] = 0.01; //Set Precision of Minuit to 1E-4
  //gMinuit->mnexcm("SET EPS", arglist, 1, ierflg);

  std::vector<double> v;
  v.reserve(5);
  v.push_back(beam->getPosition()->getX());
  v.push_back(beam->getPosition()->getY());
  v.push_back(beam->getDirection()->getTheta());
  v.push_back(beam->getDirection()->getPhi());
  v.push_back(beam->getBeta());

  setMean(mean);
  setPixelHit(pH);
  setVariables(v);
  double vSig[5] = {5E-7, 5E-7, 0.5 / 1000, 7E-7, 0.13 * beam->getMomentum()};
  double particle = beam->getParticleMass();
  vSig[4] = vSig[4] * pow(particle, 2) / pow(particle * particle + beam->getMomentum() * beam->getMomentum(), 1.5);
  //static Double_t step[5] = {, , 1, 1, 1};

  gMinuit->mnparm(0, "x", v[0], 0.01, 0, 0, ierflg);                //0, 2 * 0.097
  gMinuit->mnparm(1, "y", v[1], 0.01, 0, 0, ierflg);                //0, 2 * 0.097
  gMinuit->mnparm(2, "t", v[2], 0.01, -0.4, 0.4, ierflg);                //-0.4, 0.4,
  gMinuit->mnparm(3, "p", v[3], 0.01, 0, 0, ierflg);                //0, 2 * TMath::Pi()
  gMinuit->mnparm(4, "bm", v[4], 0.01, 1 / 1.012, 0.99992, ierflg); // 1 / 1.012, 1.001, ierflg); //0, 1

  gMinuit->mnparm(5, "errX", vSig[0], 0, 0, 0, ierflg);
  gMinuit->mnparm(6, "errY", vSig[1], 0, 0, 0, ierflg);
  gMinuit->mnparm(7, "errT", vSig[2], 0, 0, 0, ierflg);
  gMinuit->mnparm(8, "errP", vSig[3], 0, 0, 0, ierflg);
  gMinuit->mnparm(9, "errBM", vSig[4], 0, 0, 0, ierflg);

  //gMinuit->FixParameter(1);
  //gMinuit->FixParameter(2);
  //gMinuit->FixParameter(3);
  //  gMinuit->FixParameter(4);

  gMinuit->FixParameter(5);
  gMinuit->FixParameter(6);
  gMinuit->FixParameter(7);
  gMinuit->FixParameter(8);
  gMinuit->FixParameter(9);

  gMinuit->mnparm(10, "sig", sig, 0, 0, 0, ierflg);
  gMinuit->FixParameter(10);

  this->setMinuit(gMinuit.get()); 

  double ln, edm, errdef;
  int nvpar, nparx, icstatMigrad, icstatHesse, icstatSimplex = 0;

  int err_M, err_H, err_S = 0;
  arglist[0] = 500;
  gMinuit->mnexcm("MIGRAD", arglist, 1, ierflg);
  err_M = ierflg;
  gMinuit->mnstat(ln, edm, errdef, nvpar, nparx, icstatMigrad);
  gMinuit->mnexcm("HESSE", arglist, 1, ierflg);
  err_H = ierflg;
  gMinuit->mnstat(ln, edm, errdef, nvpar, nparx, icstatHesse);
  
  //gMinuit->mnexcm("SIMPLEX",arglist,0,ierflg);
  //gMinuit->mnstat(ln, edm, errdef, nvpar, nparx, icstatSimplex);
  //err_S = ierflg;

  //gMinuit->mnparm(4, "bm", v[4], 0.01, 0, 0, ierflg);
  //gMinuit->mnexcm("MIGRAD", arglist, 1, ierflg);
  //gMinuit->mnstat(ln, edm, errdef, nvpar, nparx, icstatMigrad);
  //gMinuit->mnexcm("HESSE", arglist, 1, ierflg);
  //gMinuit->mnstat(ln, edm, errdef, nvpar, nparx, icstatHesse);

  double x,y,t,p,b = 0;
  double xErr, yErr, tErr, pErr, bErr= 0;

  gMinuit->GetParameter(0, x, xErr);
  gMinuit->GetParameter(1, y, yErr);
  gMinuit->GetParameter(2, t, tErr);
  gMinuit->GetParameter(3, p, pErr);
  gMinuit->GetParameter(4, b, bErr);

  minuiVar.push_back(x);
  minuiVar.push_back(y);
  minuiVar.push_back(t);
  minuiVar.push_back(p);
  minuiVar.push_back(b);

  minuiVarErr.push_back(xErr);
  minuiVarErr.push_back(yErr);
  minuiVarErr.push_back(tErr);
  minuiVarErr.push_back(pErr);
  minuiVarErr.push_back(bErr);

  /*
  name += "_MinuitGraph.root";
  TFile *f = new TFile(name.c_str(), "recreate");
  int i = 5;
    //get the graph and write it to the root file:
    //int i = 0;
  arglist[0] = i;
  arglist[1] = 90;
  arglist[2] = .9992; //minuiVar[i] - 3 * minuiVarErr[i];
  arglist[3] = .9999; //minuiVar[i] + 3 * minuiVarErr[i];
    //if (arglist[3]>=.99992) {arglist[3]=.99992;}    gMinuit->mnexcm("SCAN", arglist, 4, ierflg);
  gMinuit->mnexcm("SCAN", arglist, 4, ierflg);
  TGraph *gr = (TGraph *)gMinuit->GetPlot();
  gr->SetTitle(std::to_string(i).c_str());
  gr->Write();
  f->Close();
*/
  std::cout << "DT: " << t << "+/-" << tErr << "   ";
  std::cout << "DB: " << b << "+/-" << bErr << "   ";
  std::cout << "LL: " << ln << "   ";
  std::cout << "WM: " << err_M << "   ";
  std::cout << "WH: " << err_H << "   ";
  std::cout << "WS: " << err_S << "    ";
  std::cout << "SM: " << icstatMigrad << "   ";
  std::cout << "SH: " << icstatHesse << "   ";
  std::cout << "SS: " << icstatSimplex << "   ";
  minuiWarning.push_back(err_M);
  minuiWarning.push_back(err_H);
  minuiWarning.push_back(icstatMigrad);
  minuiWarning.push_back(icstatHesse);

  return ln;
}

double Detector::calcLikelihoodForEvent(double pH_event, double mean_pixel) {
  double threshold = 1;
  mean_pixel = mean_pixel + pmt->getDarkNoise();
  // std::cout << mean_pixel << std::endl;
  if (pH_event >= threshold) {
    return log(1 - exp(-mean_pixel));
  } else {
    return -1.0 * mean_pixel;
  }
}

//double Detector::getNuisanceParameter(double event, double mean, double sig, double eventNumb) {

//double threshold = 1;
//mean = mean + pmt->getDarkNoise();
//if (event >= threshold) { //Hit
////Must be performed numerically
//TF1 f("Hit Function", "x + ([0] * exp( -[1] - x))/(exp(-[1]-x)-1)", -0.1, 1);
//f.SetParNames("sig^2", "mean");
//f.SetParameter(0, std::pow(sig, 2));
//f.SetParameter(1, mean);
//ROOT::Math::WrappedTF1 wf1(f);

//// Create the Integrator
//ROOT::Math::BrentRootFinder brf;

//// Set parameters of the method
//brf.SetFunction(wf1, 0, .5);
//brf.Solve();
//double w = brf.Root();
//f.Delete();
//if (mean + w < 0) {
//return 0;
//} else {
//return w;
//}
//} else { //NoHit
//return -pow(sig, 2);
//}
//}

//double Detector::getNuisanceParameterMINUIT(double event, double mean, double sig, double eventNumb) {
//double threshold = 1;
//mean = mean + pmt->getDarkNoise();
//if (event >= threshold) { //Hit
//TMinuit *gMinuit = new TMinuit(2);
//gMinuit->SetPrintLevel(-1);
//gMinuit->SetFCN(FCN);

//double arglist[10];
//int ierflg = 0;

//arglist[0] = 0;
//gMinuit->mnexcm("SET STR", arglist, 1, ierflg);

//arglist[0] = 1;
//gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

//static Double_t vstart[3] = {mean, 0.001, sig};
//static Double_t step[3] = {vstart[0] / 10, vstart[1] / 10, vstart[2] / 10};
//gMinuit->mnparm(0, "mean", vstart[0], step[0], 0, 0, ierflg);
//gMinuit->mnparm(1, "w", vstart[1], step[1], 0, 0, ierflg);
//gMinuit->mnparm(2, "sig", vstart[2], step[2], 0, 0, ierflg);

//gMinuit->FixParameter(0);
//gMinuit->FixParameter(2);

//arglist[0] = 500;
//gMinuit->mnexcm("MIGRAD", arglist, 1, ierflg);
//double w;
//double wErr;
//gMinuit->GetParameter(1, w, wErr);
//if (mean + w < 0) {
//return 0;
//} else {
//return w;
//}
//} else { //NoHit
//return -pow(sig, 2);
//}
//}

std::vector<double> Detector::getLogLikelihood(std::vector<std::vector<double>> pH,
                                               std::vector<double> mean1,
                                               std::vector<double> mean2) {
  std::vector<double> ll_pH;
  for (int rowNo = 0; rowNo < pH.size(); rowNo++) {
    std::vector<double> pH_event = pH.at(rowNo);
    double ll_event = 0.0;
    for (int colNo = 0; colNo < pH_event.size(); colNo++) {
      ll_event = ll_event + -2 * calcLikelihoodForEvent(pH_event.at(colNo), mean1.at(colNo)) +
                 2 * calcLikelihoodForEvent(pH_event.at(colNo), mean2.at(colNo));
    }
    ll_pH.push_back(ll_event);
  }
  return ll_pH;
}

//**********************
//THIS IS JUST A DUMB TEST AND SHOULD BE DELETED
std::vector<double> Detector::getLogLikelihoodStupid(std::vector<std::vector<double>> pH,
                                                     std::vector<double> mean1,
                                                     std::vector<double> mean2) {
  std::vector<double> ll_pH;
  for (int rowNo = 0; rowNo < pH.size(); rowNo++) {
    std::vector<double> pH_event = pH.at(rowNo);
    double ll_event = 0.0;
    for (int colNo = 0; colNo < pH_event.size(); colNo++) {
      if (mean1.at(colNo) < .001) {
        mean1.at(colNo) = 1E-20;
      }
      if (mean2.at(colNo) < .001) {
        mean2.at(colNo) = 1E-20;
      }
      ll_event = ll_event + -2 * calcLikelihoodForEvent(pH_event.at(colNo), mean1.at(colNo)) +
                 2 * calcLikelihoodForEvent(pH_event.at(colNo), mean2.at(colNo));
    }
    ll_pH.push_back(ll_event);
  }
  return ll_pH;
}
//**********************

std::vector<double> Detector::getWeightedLogLikelihood(std::vector<std::vector<double>> pH,
                                                       std::vector<double> mean1,
                                                       std::vector<double> mean2) {
  std::vector<double> weights1 = calcWeight(mean1, pH.size());
  double weightSum1 = calcWeightSum(weights1);
  std::vector<double> weights2 = calcWeight(mean2, pH.size());
  double weightSum2 = calcWeightSum(weights2);

  std::vector<double> ll_pH;
  for (int rowNo = 0; rowNo < pH.size(); rowNo++) {
    std::vector<double> pH_event = pH.at(rowNo);
    double ll_event = 0.0;
    for (int colNo = 0; colNo < pH_event.size(); colNo++) {
      // std::cout << "WI1:" << weights1[colNo] << "    WS1:" << weightSum1 << "
      // W1:" << weights1[colNo]/weightSum1 << "     ";  std::cout << "WI2:" <<
      // weights2[colNo] << "    WS2:" << weightSum2 << "    w2:" <<
      // weights2[colNo]/weightSum2 << std::endl;

      ll_event = ll_event +
                 -2 * weights1[colNo] * calcLikelihoodForEvent(pH_event.at(colNo), mean1[colNo]) /
                     weightSum1 +
                 2 * weights2[colNo] * calcLikelihoodForEvent(pH_event.at(colNo), mean2[colNo]) /
                     weightSum2;
    }
    ll_pH.push_back(ll_event);
  }
  return ll_pH;
}

std::vector<std::vector<double>> Detector::getLogLikelihoodPerPixel(std::vector<std::vector<double>> pH, std::vector<double> mean1,
                                                                    std::vector<double> mean2) {
  std::vector<std::vector<double>> ll_pH_perPixel;
  for (int eventCount = 0; eventCount < pH.size(); eventCount++) {
    std::vector<double> pH_event = pH.at(eventCount);
    double ll_event = 0.0;
    std::vector<double> temp;
    for (int pixelCount = 0; pixelCount < pH_event.size(); pixelCount++) {
      temp.push_back(-2 * calcLikelihoodForEvent(pH_event.at(pixelCount), mean1.at(pixelCount)) +
                     2 * calcLikelihoodForEvent(pH_event.at(pixelCount), mean2.at(pixelCount)));
    }
    ll_pH_perPixel.push_back(temp);
  }
  return ll_pH_perPixel;
}

//std::vector<std::vector<double>> Detector::getLogLikelihoodNuissancePerPixel(std::vector<std::vector<double>> pH, std::vector<double> mean1,
//std::vector<double> mean2) {
//std::vector<std::vector<double>> ll_pH_perPixel;
//for (int rowNo = 0; rowNo < pH.size(); rowNo++) {
//std::vector<double> pH_event = pH.at(rowNo);
//double ll_event = 0.0;
//std::vector<double> temp;
//for (int colNo = 0; colNo < pH_event.size(); colNo++) {
//double sig1 = .00315;
//double w1 = getNuisanceParameter(pH_event.at(colNo), mean1.at(colNo), sig1, pH.size());
//double sig2 = .002851;
//double w2 = getNuisanceParameter(pH_event.at(colNo), mean2.at(colNo), sig2, pH.size());
//temp.push_back(-2 * calcLikelihoodForEvent(pH_event.at(colNo), mean1.at(colNo) + w1) + pow(w1 / sig1, 2) +
//2 * calcLikelihoodForEvent(pH_event.at(colNo), mean2.at(colNo) + w2) - pow(w2 / sig2, 2));
//}
//ll_pH_perPixel.push_back(temp);
//}
//return ll_pH_perPixel;
//}

//std::vector<std::vector<double>> Detector::getWeightedLogLikelihoodPerPixel(std::vector<std::vector<double>> pH,
//std::vector<double> mean1, std::vector<double> mean2) {
//std::vector<double> weights1 = calcWeight(mean1, pH.size());
//double weightSum1 = calcWeightSum(weights1);
//std::vector<double> weights2 = calcWeight(mean2, pH.size());
//double weightSum2 = calcWeightSum(weights2);

//std::vector<std::vector<double>> ll_pH_perPixel;
//for (int rowNo = 0; rowNo < pH.size(); rowNo++) {
//std::vector<double> pH_event = pH.at(rowNo);
//double ll_event = 0.0;
//std::vector<double> temp;
//for (int colNo = 0; colNo < pH_event.size(); colNo++) {
//temp.push_back(-2 * weights1[colNo] *
//calcLikelihoodForEvent(pH_event.at(colNo), mean1.at(colNo)) / weightSum1 +
//2 * weights2[colNo] *
//calcLikelihoodForEvent(pH_event.at(colNo), mean2.at(colNo)) / weightSum2);
//// ll_event = ll_pH_perPixel[colNo][rowNo]+ ll_event;
//// std::cout << colNo << "/" << rowNo << "      ";
//// std::cout << mean1.at(colNo) << "   " << mean2.at(colNo) << "   ";
//// std::cout <<
//// (-2*calcLikelihoodForEvent(pH_event.at(colNo),mean1.at(colNo))+
//// 2*calcLikelihoodForEvent(pH_event.at(colNo),mean2.at(colNo))) <<
//// std::endl;
//}
//ll_pH_perPixel.push_back(temp);
//}
//return ll_pH_perPixel;
//}

std::vector<double> Detector::getLogLikelihoodValuesForSingleParticle(std::vector<std::vector<double>> pH,
                                                                      std::vector<double> mean) {
  std::vector<double> ll_pH;
  for (int rowNo = 0; rowNo < pH.size(); rowNo++) {
    std::vector<double> pH_event = pH.at(rowNo);
    double ll_event = 0.0;
    for (int colNo = 0; colNo < pH_event.size(); colNo++) {
      ll_event = ll_event + -2 * calcLikelihoodForEvent(pH_event.at(colNo), mean.at(colNo));
    }
    ll_pH.push_back(ll_event);
  }
  return ll_pH;
}

double Detector::getLogLikelihoodValuesForSingleEventSingleMean(std::vector<double> pH_event,
                                                                double mean) {
  double ll_event = 0.0;
  for (int colNo = 0; colNo < pH_event.size(); colNo++) {
    ll_event = ll_event + -2 * calcLikelihoodForEvent(pH_event.at(colNo), mean);
  }
  return ll_event;
}

std::vector<double> Detector::getMeanPerPixel(std::vector<std::vector<double>> hit_on_pixels_perEvent) {

  std::vector<double> currentMean(hit_on_pixels_perEvent.at(0).size(),0);
  for (int i = 0; i < hit_on_pixels_perEvent.size(); i++) {
    std::vector<double> pHRow = hit_on_pixels_perEvent.at(i);
    for (int j = 0; j < pHRow.size(); j++) {
      currentMean.at(j) = currentMean.at(j) + pHRow.at(j) / hit_on_pixels_perEvent.size();
    }
  }
  this->mean = currentMean;
  return currentMean;
}

std::vector<double> Detector::calcMeanRelativeUncertainty(std::vector<double> mean,
                                                          int eventCount) {
  std::vector<double> uncertainty;
  for (int i = 0; i < mean.size(); i++) {
    uncertainty.push_back(sqrt(mean[i]) / mean[i]);
  }
  return uncertainty;
}

std::vector<double> Detector::calcWeight(std::vector<double> mean, int eventCount) {
  std::vector<double> weight;
  for (int i = 0; i < mean.size(); i++) {
    weight.push_back(1.0 / std::pow(sqrt(mean[i]) / mean[i], 2));
  }
  return weight;
}

std::vector<double> Detector::getStandardDeviationForPixel(std::vector<std::vector<double>> pH) {
  std::vector<double> std_pixel;
  std::vector<double> eventWeight(pH.size(), 1.0);

  for (int pixelCount = 0; pixelCount < pH[0].size(); pixelCount++) {
    std::vector<double> pixelEvents = getEventsForPixel(pH, pixelCount);
    std::vector<double> pE_doubles(pixelEvents.begin(),pixelEvents.end());
    std::shared_ptr<TH1D> temp = std::make_shared<TH1D>(("pixel" + std::to_string(pixelCount)).c_str(),
                               ("pixel" + std::to_string(pixelCount)).c_str(), 100, 0, 0);
    temp->FillN(pE_doubles.size(), &pE_doubles[0], &eventWeight[0]);
    std_pixel.push_back(temp->GetStdDev());
  }
  return std_pixel;
}

std::vector<double> Detector::getEventsForPixel(std::vector<std::vector<double>> pH, int pixel) {
  std::vector<double> pHPixel;
  for (int eventCount = 0; eventCount < pH.size(); eventCount++) {
    for (int pixelCount = 0; pixelCount < pH[eventCount].size(); pixelCount++) {
      if (pixelCount == pixel) {
        pHPixel.push_back(pH[eventCount][pixelCount]);
      }
    }
  }
  return pHPixel;
}

double Detector::calcWeightSum(std::vector<double> mean, int eventCount) {
  double weightSum = 0;
  for (int i = 0; i < mean.size(); i++) {
    weightSum += (1 / std::pow(sqrt(mean[i]) / sqrt(eventCount), 2));
  }
  return weightSum;
}

double Detector::calcWeightSum(std::vector<double> weights) {
  double weightSum = 0;
  for (int i = 0; i < weights.size(); i++) {
    weightSum += weights[i];
  }
  return weightSum;
}

//void Detector::testSeperationLikelihoodOnVariable(
//int steps, std::shared_ptr<InterpolateLikelihood> linearInter, std::string variableName,
//std::vector<double> inputVar1, std::vector<double> inputVar2,
//std::vector<std::vector<double>> hit_on_pixels_perEvent1,
//std::vector<std::vector<double>> hit_on_pixels_perEvent2, std::vector<double> &delta1,
//std::vector<double> &delta2, std::vector<double> &delta_particle1,
//std::vector<double> &delta_particle2) {
//std::vector<std::vector<double>> inputVariableList1;
//std::vector<std::vector<double>> inputVariableList2;
//std::vector<double> dependentVariable1;
//std::vector<double> dependentVariable2;

//double minV = 0;
//double maxV = 0;
//int stepV = 100;
//std::vector<bool> variableSwitch;
//variableSwitch.push_back(false);
//variableSwitch.push_back(false);
//variableSwitch.push_back(false);
//variableSwitch.push_back(false);
//variableSwitch.push_back(false);

//double inputValue1;
//double inputValue2;

//if (variableName == std::string("x")) {
//minV = linearInter->getMinX();
//maxV = linearInter->getMaxX();
//inputValue1 = inputVar1.at(0);
//inputValue2 = inputVar2.at(0);
//variableSwitch.at(0) = true;
//} else if (variableName == std::string("y")) {
//minV = linearInter->getMinY();
//maxV = linearInter->getMaxY();
//inputValue1 = inputVar1.at(1);
//inputValue2 = inputVar2.at(1);
//variableSwitch.at(1) = true;
//} else if (variableName == std::string("theta")) {
//minV = linearInter->getMinTheta();
//maxV = linearInter->getMaxTheta();
//inputValue1 = inputVar1.at(2);
//inputValue2 = inputVar2.at(2);
//variableSwitch.at(2) = true;
//} else if (variableName == std::string("phi")) {
//minV = linearInter->getMinPhi();
//maxV = linearInter->getMaxPhi();
//inputValue1 = inputVar1.at(3);
//inputValue2 = inputVar2.at(3);
//variableSwitch.at(3) = true;
//} else if (variableName == std::string("beta")) {
//minV = linearInter->getMinBeta();
//maxV = linearInter->getMaxBeta();
//inputValue1 = inputVar1.at(4);
//inputValue2 = inputVar2.at(4);
//variableSwitch.at(4) = true;
//} else {
//return;
//}
//getListOfVariables(steps, inputVar1.size(), minV, maxV, variableSwitch, inputVar1, inputVar2,
//dependentVariable1, dependentVariable2, inputVariableList1,
//inputVariableList2);

//std::vector<double> meanInterp1 = linearInter->interpolateMeanToFirstOrder(inputVar1);
//std::vector<double> meanInterp2 = linearInter->interpolateMeanToFirstOrder(inputVar2);

//// std::vector<double> llInterp1 =
//// getLogLikelihood(hit_on_pixels_perEvent1,meanInterp1,meanInterp2);
//// std::vector<double> llInterp2 =
//// getLogLikelihood(hit_on_pixels_perEvent2,meanInterp1,meanInterp2);

//std::vector<double> min_ll1(hit_on_pixels_perEvent1.size(), 1E9);
//std::vector<double> min_ll2(hit_on_pixels_perEvent1.size(), 1E9);

//std::vector<double> min_ll1_index(hit_on_pixels_perEvent1.size(), 0);
//std::vector<double> min_ll2_index(hit_on_pixels_perEvent1.size(), 0);

//std::vector<double> min_ll1_particle1(hit_on_pixels_perEvent1.size(), 1E9);
//std::vector<double> min_ll2_particle2(hit_on_pixels_perEvent1.size(), 1E9);

//std::vector<double> min_ll1_particle1_index(hit_on_pixels_perEvent1.size(), 0);
//std::vector<double> min_ll2_particle2_index(hit_on_pixels_perEvent1.size(), 0);

//for (unsigned long k = 0; k < inputVariableList1.at(0).size(); k++) {
//std::vector<double> tempVariableList1;
//tempVariableList1.reserve(inputVariableList1.size());
//std::vector<double> tempVariableList2;
//tempVariableList2.reserve(inputVariableList2.size());
//for (unsigned long j = 0; j < inputVariableList1.size(); j++) {
//tempVariableList1.push_back(inputVariableList1[j][k]);
//tempVariableList2.push_back(inputVariableList2[j][k]);
//}
//std::vector<double> tempMean1 = linearInter->interpolateMeanToFirstOrder(tempVariableList1);
//std::vector<double> tempMean2 = linearInter->interpolateMeanToFirstOrder(tempVariableList2);
//std::vector<double> ll1 = getLogLikelihood(hit_on_pixels_perEvent1, tempMean1, meanInterp2);
//std::vector<double> ll2 = getLogLikelihood(hit_on_pixels_perEvent2, meanInterp1, tempMean2);

//std::vector<double> ll1_forParticle1 =
//getLogLikelihoodValuesForSingleParticle(hit_on_pixels_perEvent1, tempMean1);
//std::vector<double> ll2_forParticle2 =
//getLogLikelihoodValuesForSingleParticle(hit_on_pixels_perEvent2, tempMean2);

//std::cout << ll1.size() << "  " << ll2.size() << std::endl;
//for (int v = 0; v < ll1.size(); v++) {
//if (ll1.at(v) < min_ll1.at(v)) {
//min_ll1.at(v) = ll1.at(v);
//min_ll1_index.at(v) = k;
//}
//if (ll2.at(v) < min_ll2.at(v)) {
//min_ll2.at(v) = ll2.at(v);
//min_ll2_index.at(v) = k;
//}

//if (ll1_forParticle1.at(v) < min_ll1_particle1.at(v)) {
//min_ll1_particle1.at(v) = ll1_forParticle1.at(v);
//min_ll1_particle1_index.at(v) = k;
//}
//if (ll2_forParticle2.at(v) < min_ll2_particle2.at(v)) {
//min_ll2_particle2.at(v) = ll2_forParticle2.at(v);
//min_ll2_particle2_index.at(v) = k;
//}

//std::cout << v << "  Particle 1"
//<< "   ";
//std::cout << dependentVariable1.size() << "    ";
//std::cout << inputVariableList1.at(0).size() << "    ";
//std::cout << "In: " << inputValue1 << "    ";
//std::cout << ll1.at(v) << "    ";
//std::cout << min_ll1_index.at(v) << "   ";
//std::cout << dependentVariable1.at(min_ll1_index.at(v)) << "   ";
//std::cout << min_ll1_particle1_index.at(v) << "   ";
//std::cout << dependentVariable1.at(min_ll1_particle1_index.at(v)) << std::endl;
//std::cout << v << "  Particle 2"
//<< "   ";
//std::cout << dependentVariable2.size() << "    ";
//std::cout << inputVariableList2.at(0).size() << "    ";

//std::cout << "In: " << inputValue2 << "    ";
//std::cout << ll2.at(v) << "    ";
//std::cout << min_ll2_index.at(v) << "   ";
//std::cout << dependentVariable2.at(min_ll2_index.at(v)) << "   ";
//std::cout << min_ll2_particle2_index.at(v) << "   ";
//std::cout << dependentVariable2.at(min_ll2_particle2_index.at(v)) << std::endl;
//}
//}
//for (int i = 0; i < min_ll2_index.size(); i++) {
//delta1.push_back(dependentVariable1.at(min_ll1_index.at(i)) - inputValue1);
//delta2.push_back(dependentVariable2.at(min_ll2_index.at(i)) - inputValue2);
//delta_particle1.push_back(dependentVariable1.at(min_ll1_particle1.at(i)) - inputValue1);
//delta_particle2.push_back(dependentVariable2.at(min_ll2_particle2.at(i)) - inputValue2);
//}
//return;
//}
double Detector::calcStepSize(double start, double end, double step) {
  return (end - start) / (step - 1);
}

void Detector::getListOfVariables(int stepV, int varCount, double minV, double maxV,
                                  std::vector<bool> variableSwitch,
                                  std::vector<double> variablesPion,
                                  std::vector<double> variablesKaon,
                                  std::vector<double> &dependentVariablePion,
                                  std::vector<double> &dependentVariableKaon,
                                  std::vector<std::vector<double>> &inputVariableListPion,
                                  std::vector<std::vector<double>> &inputVariableListKaon) {
  double deltaV;
  dependentVariablePion.reserve(stepV - 10);
  dependentVariableKaon.reserve(stepV - 10);
  for (unsigned long k = 0; k < varCount; k++) {
    if (variableSwitch.at(k)) {
      std::vector<double> tempVListPion;
      std::vector<double> tempVListKaon;
      tempVListPion.reserve(stepV - 10);
      tempVListKaon.reserve(stepV - 10);
      deltaV = calcStepSize(minV, maxV, stepV);
      for (double v = minV; v < maxV; v = v + deltaV) {
        tempVListPion.push_back(v);
        tempVListKaon.push_back(v);
      }
      dependentVariablePion = tempVListPion;
      dependentVariableKaon = tempVListKaon;
      inputVariableListPion.push_back(tempVListPion);
      inputVariableListKaon.push_back(tempVListKaon);
    } else {
      std::vector<double> tempVListPion;
      std::vector<double> tempVListKaon;
      tempVListKaon.reserve(stepV);
      tempVListPion.reserve(stepV);
      for (int j = 0; j < stepV; j++) {
        tempVListPion.push_back(variablesPion.at(k));
        tempVListKaon.push_back(variablesKaon.at(k));
      }
      inputVariableListPion.push_back(tempVListPion);
      inputVariableListKaon.push_back(tempVListKaon);
    }
  }
  return;
}
