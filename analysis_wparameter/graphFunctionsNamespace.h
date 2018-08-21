//# CPUPROFILE=./tmp/main.cpp ./main.cpp
//C++
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <numeric>
#include <stdio.h>
#include <stdlib.h> /* srand, rand */
#include <string>
#include <sys/stat.h>
#include <time.h> /* time */
#include <vector>

//Boost
//include "/usr/local/include/boost/algorithm/string.hpp"

//Root
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGFrame.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TPDF.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TRandom2.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVirtualPad.h"

//Classes
#include "aerogelClass.h"
#include "detectedEventClass.h"
#include "detectorClass.h"
#include "directionClass.h"
#include "interpolateLikelihoodClass.h"
#include "pmtClass.h"
#include "positionClass.h"
#include "randomEventClass.h"

namespace graphingAndAnalysisCherenkovPID {

void printArgumentHelperList(){
  std::cout << "1st Argument: filename (extension to name)" << std::endl;
  std::cout << "2nd Argument: number of events" << std::endl;
  std::cout << "3rd Argument: type of interpolation" << std::endl;
  std::cout << "\t 0 = LINEAR \n\t 1 = SPLINE" << std::endl;
  std::cout << "4th Argument: w parameter 0/1 = off/on" << std::endl;
  std::cout << "5th Argument: x value" << std::endl; 
  std::cout << "6th Argument: y value" << std::endl;
  std::cout << "7th Argument: theta value" << std::endl;
  std::cout << "8th Argument: phi value" << std::endl;
  std::cout << "9th Argument: beta bool value" << std::endl;
  std::cout << "10th Argument: beta value" << std::endl;
}

void printArgList(int argc, char **argv){

for (int i = 0; i < argc-1; i++)
{
std::cout << i +1 << "  " << argv[i+1] << std::endl;
}

}

/*
 * saveCanvas = save histograms into seperate  files 
 */
void saveCanvas(std::string fileName, std::shared_ptr<TCanvas> canvas) {
  std::string namePng = fileName + ".png";
  std::string nameRoot = fileName + ".root";
  std::string namePdf = fileName + ".pdf";
  canvas->SaveAs(namePng.c_str());
  canvas->SaveAs(nameRoot.c_str());
  canvas->SaveAs(namePdf.c_str());
  return;
}

/*
 * GetCumulative = root6 getCumulative
 */

 TH1D* GetCumulative(TH1D*h, bool forward=true, const char* suffix="_cumulative")
  {
     const Int_t nbinsx = h->GetNbinsX();
     const Int_t nbinsy = h->GetNbinsY();
     const Int_t nbinsz = h->GetNbinsZ();
     TH1D* hintegrated = (TH1D*) h->Clone(TString(h->GetName())+suffix);
     hintegrated->Reset();
     if (forward) { // Forward computation
        Double_t sum = 0.;
        for (Int_t binz = 1; binz <= nbinsz; ++binz) {
      for (Int_t biny = 1; biny <= nbinsy; ++biny) {
         for (Int_t binx = 1; binx <= nbinsx; ++binx) {
            const Int_t bin = hintegrated->GetBin(binx, biny, binz);
            sum += h->GetBinContent(bin);
            hintegrated->SetBinContent(bin, sum);
         }
      }
        }
     } else { // Backward computation
        Double_t sum = 0.;
        for (Int_t binz = nbinsz; binz >= 1; --binz) {
      for (Int_t biny = nbinsy; biny >= 1; --biny) { 
      for (Int_t binx = nbinsx; binx >= 1; --binx) {
            const Int_t bin = hintegrated->GetBin(binx, biny, binz);
            sum += h->GetBinContent(bin);
            hintegrated->SetBinContent(bin, sum);
         }
      }
        }
     }
    return hintegrated;
  }


/*
 * getSeperationValue = calculates the seperation value from the two particle histograms 
 * The seperation value is considered overlap in the area of the second particle past the 90%
 * area point of the second particle
 */
double getSeperationValue(std::shared_ptr<TH1D> h1, std::shared_ptr<TH1D> h2) {
  double minX1 = h1->GetXaxis()->GetXmin();
  double maxX1 = h1->GetXaxis()->GetXmax();

  double minX2 = h2->GetXaxis()->GetXmin();
  double maxX2 = h2->GetXaxis()->GetXmax();

  double area = 0;
  double x1 = minX1;
  double xCount = 100;
  double dx = (maxX1 - minX1) / (xCount - 1);

  std::shared_ptr<TH1> hC1(GetCumulative(h1.get()));
  std::shared_ptr<TH1> hC2(GetCumulative(h2.get()));

  while (area < .900 && x1 <= maxX1) {
    area = hC1->Interpolate(x1);
    x1 = x1 + dx;
  }
  return 1 - hC2->Interpolate(x1);
}

//##########################################################################################################
/*
 *getFileName functions = takes different inputs and genreates the same fileName
 */
std::string getFileName(bool betaBool, double x, double y, double theta, double phi, double betaMomentum) {
  std::string fileName;
  if (betaBool) {
    fileName = fileName + "beta" + std::to_string(betaMomentum);
  } else {
    fileName = fileName + "momentum" + std::to_string(betaMomentum);
  }

  fileName = fileName + "_xPosition" + std::to_string(x);
  fileName = fileName + "_yPosition" + std::to_string(y);
  fileName = fileName + "_theta" + std::to_string(theta);
  fileName = fileName + "_phi" + std::to_string(phi);
  return fileName;
}

std::string getFileName(bool betaBool, std::vector<double> variable) {
  std::string fileName;
  if (betaBool) {
    fileName = fileName + "beta" + std::to_string(variable[4]);
  } else {
    fileName = fileName + "momentum" + std::to_string(variable[4]);
  }

  fileName = fileName + "_xPosition" + std::to_string(variable[0]);
  fileName = fileName + "_yPosition" + std::to_string(variable[1]);
  fileName = fileName + "_theta" + std::to_string(variable[2]);
  fileName = fileName + "_phi" + std::to_string(variable[3]);
  return fileName;
}

std::string getFileName(bool betaBool, std::shared_ptr<Beam> beam) {
  std::string fileName;
  if (betaBool) {
    fileName = fileName + "beta" + std::to_string(beam->getBeta());
  } else {
    fileName = fileName + "momentum" + std::to_string(beam->getMomentum());
  }

  fileName = fileName + "_xPosition" + std::to_string(beam->getPosition()->getX());
  fileName = fileName + "_yPosition" + std::to_string(beam->getPosition()->getY());
  fileName = fileName + "_theta" + std::to_string(beam->getDirection()->getTheta());
  fileName = fileName + "_phi" + std::to_string(beam->getDirection()->getPhi());
  return fileName;
}
//##########################################################################################################

/*
 * saveAllInFolder = based of the bool saveALL choose whether the 
 * file will be saved in a lump folder for that configurations 
 * or saved with respect to its attribute being a likelihood
 */

std::string saveAllInFolder(std::string type, std::string folderName, bool saveAll) {
  std::string pathName;
  if (!saveAll) {
    pathName = "./Images/" + type + "/";
    mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  } else {
    pathName = "./Images/" + folderName + "/";
    mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }
  return pathName;
}

/*
 * plotBothLLHistograms = plot histogram of each likelihood either in two 
 * seperate pads or in a single graph
 */
void plotBothLLHistograms(std::shared_ptr<TH1D> h1, std::shared_ptr<TH1D> h2,
                          bool showSubGraphs,
                          std::string particle1, std::string particle2,
                          std::string variableNames,
                          std::string appendName = "",
                          bool saveAll = false) {

  double area1 = h1->Integral();
  double area2 = h2->Integral();
  h1->Scale(1 / (area1));
  h2->Scale(1 / (area2));
  h1->SetTitle(particle1.c_str());
  h2->SetTitle(particle2.c_str());
  //h1->SetFillColorAlpha(38, .5);
  //h2->SetFillColorAlpha(46, .5);

  std::string pathName = saveAllInFolder("logLikelihood", variableNames, saveAll);

  //Create Path Name and File names
  mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  std::string fileName = pathName + "ll";
  fileName = fileName + variableNames;
  fileName = fileName + appendName;

  double minX1 = h1->GetXaxis()->GetXmin();
  double maxX1 = h1->GetXaxis()->GetXmax();

  double minX2 = h2->GetXaxis()->GetXmin();
  double maxX2 = h2->GetXaxis()->GetXmax();

  std::string statTest = std::to_string(getSeperationValue(h1, h2));

  //Canvases
  std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>(fileName.c_str(), "", 1200, 1200));

  /*
   * if the minimum of the second particle's histograms is positive then it suggests that there is no overlap
   * between the two histograms so they should be plotted in different graph 
   */
  if (minX2 > 0) {
    //Set Global Style Variables
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    //Create 2 Pads
    c1->SetGrid();
    c1->Divide(2, 0, 0, 0);

    std::shared_ptr<TPaveLabel> temp1(std::make_shared<TPaveLabel>(0.1, 0.96, 0.9, 0.99, ""));
    temp1->Draw();

    //Draw Histogram 1
    c1->cd(1);
    c1->SetGrid();
    h1->GetXaxis()->SetTitle("");
    h1->GetYaxis()->SetTitle("");
    h1->Draw();
    c1->Modified();
    c1->Update();

    //Draw Histogram 2
    c1->cd(2);
    c1->SetGrid();
    gPad->SetTicky(2);
    h2->GetYaxis()->SetLabelOffset(999);
    h2->GetYaxis()->SetLabelSize(0);
    h2->GetYaxis()->SetTitle("");
    h2->Draw();
    h2->GetXaxis()->SetNdivisions(-3);
    c1->Modified();
    c1->Update();

    //Create Legend and Draw
    c1->cd(0);
    std::shared_ptr<TLegend> leg(std::make_shared<TLegend>(0.1, 0.7, 0.48, 0.9));
    leg->SetHeader("Legend"); // option "C" allows to center the header
    leg->AddEntry(h1.get(), particle1.c_str(), "f");
    leg->AddEntry(h2.get(), particle2.c_str(), "f");
    leg->AddEntry((TObject *)0, (std::string("Seperatation: ") + statTest).c_str(), "");
    leg->Draw();

    c1->cd(0);
    //Redraw Canvas
    c1->Modified();
    c1->Update();
    saveCanvas(fileName, c1);

  } else {
    //1 Pads
    gStyle->SetOptTitle(0);
    c1->SetGrid();

    //Histogram Stack
    std::shared_ptr<THStack> hs(std::make_shared<THStack>("1 and 2", ""));
    hs->Add(h1.get());
    hs->Add(h2.get());
    hs->Draw("nostack");
    hs->GetXaxis()->SetTitle("-2 ln(#frac{L_{#pi}}{L_{k}})");

    //Create Legend and Draw
    std::shared_ptr<TLegend> leg(std::make_shared<TLegend>(0.1, 0.7, 0.48, 0.9));
    leg->SetHeader("Legend"); // option "C" allows to center the header
    leg->AddEntry(h1.get(), particle1.c_str(), "f");
    leg->AddEntry(h2.get(), particle2.c_str(), "f");
    leg->AddEntry((TObject *)0, (std::string("Seperatation: ") + statTest).c_str(), "");
    leg->Draw();

    //Pave Label
    std::shared_ptr<TPaveLabel> temp(std::make_shared<TPaveLabel>(0.1, 0.96, 0.9, 0.99, ""));
    temp->Draw();

    //Redraw Canvas
    c1->Modified();
    c1->Update();
    saveCanvas(fileName, c1);
  }

  /*
   * showSubGraphs = Switch to print each individual histogram of the likelihood
   */
  if (showSubGraphs) {
    //PARTICLE 1
    fileName = pathName + "ll" + particle1;
    std::shared_ptr<TCanvas> c2(std::make_shared<TCanvas>(pathName.c_str(), "", 1200, 1200));
    c2->SetGrid();
    h1->Draw();
    saveCanvas(fileName, c2);

    //PARTICLE 2
    fileName = pathName + "ll" + particle2;
    std::shared_ptr<TCanvas> c3(std::make_shared<TCanvas>(pathName.c_str(), "", 1200, 1200));
    c3->SetGrid();
    h2->Draw();
    saveCanvas(fileName, c3);
  }
  return;
}

/*
 * getHistogramOfLL = plot histogram of loglikelihood values 
 * give random name to each histogram
 */
std::shared_ptr<TH1D> getHistogramOfLL(std::vector<double> ll, std::string title = "", bool doDraw = true) {
  /* initialize random seed: */
  //int randTime;
  //std::srand(time(NULL));

  ////To not cause memory leak errors when bulk generating the histograms
  //[> generate secret number between 1 and 10: <]
  //randTime = std::rand() % 1000 + 1;
  std::shared_ptr<TCanvas> c1;
  std::vector<double> weights(ll.size(), 1.0);
  if (doDraw) {
    c1 = std::make_shared<TCanvas>((title).c_str());
    c1->SetGrid();
  }
  std::shared_ptr<TH1D> h1(std::make_shared<TH1D>((title).c_str(),
                                                  "Histogram of loglikelihood ratio",
                                                  100, 0, 0));

  for (int i = 0; i < ll.size(); i++) {
    h1->Fill(ll[i]);
  }
  h1->SetStats();
  h1->SetFillColor(38);
  h1->SetXTitle("-2 ln(#frac{L_{#pi}}{L_{k}})");
  h1->SetYTitle("");
  h1->Draw();
  return h1;
}

void getVariableFromSwitch(std::vector<double> variables, std::vector<bool> variableSwitch, double &value, std::string &variableName) {
  if (variableSwitch.at(0)) {
    variableName = "x";
    value = variables.at(0);
  } else if (variableSwitch.at(1)) {
    variableName = "y";
    value = variables.at(1);
  } else if (variableSwitch.at(2)) {
    variableName = "theta";
    value = variables.at(2);
  } else if (variableSwitch.at(3)) {
    variableName = "phi";
    value = variables.at(3);
  } else if (variableSwitch.at(4)) {
    variableName = "beta";
    value = variables.at(4);
  }
  return;
}

void saveGraphOfMarginalLikelihood(std::vector<double> variableList, std::vector<double> ll_perVariable,
                                   std::vector<bool> variableSwitch,
                                   std::vector<double> variables, bool betaBool, std::string yAxisTitle = "",
                                   bool saveAllInSingleFolder = false) {
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  std::string variableLabel = getFileName(betaBool, variables);
  std::string pathName = saveAllInFolder("Marginals", variableLabel, saveAllInSingleFolder);

  std::string variableName;
  double interpVal = 0;
  getVariableFromSwitch(variables, variableSwitch, interpVal, variableName);

  pathName = pathName + variableName + "/";
  mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  //Make Canvas plot Marginal
  std::string canvasTitle = "Estimate of " + variableName;
  std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>(canvasTitle.c_str()));
  std::shared_ptr<TGraphErrors> g1(std::make_shared<TGraphErrors>(variableList.size(), &variableList[0], &ll_perVariable[0]));
  c1->SetGrid();

  g1->GetXaxis()->SetTitle(variableName.c_str());
  if (yAxisTitle == std::string("")) {
    g1->GetYaxis()->SetTitle("Log LikelihoodValues");
  } else {
    g1->GetYaxis()->SetTitle(yAxisTitle.c_str());
  }
  g1->Draw("AC");

  int minPos = std::distance(ll_perVariable.begin(), min_element(ll_perVariable.begin(), ll_perVariable.end()));

  TPaveText *ps = new TPaveText(.75, .7, .95, .9, "NDC");

  ps->SetName("mystats");
  std::string text = "Min Variable: ";
  text = text + std::to_string(variableList.at(minPos));
  ps->AddText(text.c_str());
  std::string text2 = "Input Variable: ";
  text2 = text2 + std::to_string(interpVal);
  ps->AddText(text2.c_str());
  ps->Draw();

  //Print
  std::string fileName = pathName + variableName + "Marginal";
  fileName = fileName + getFileName(betaBool, variables);
  c1->Modified();
  c1->Update();
  saveCanvas(fileName, c1);
}

/*
 *  getBinEdges = Calculated the binEdges given the PMT value
 *  Input:
 *  Pointer to vector of doubles for binEdgesx
 *  Pointer to vector of doubles for binEdgesy
 */
void getBinEdges(std::shared_ptr<PMT> pmt, std::vector<double> &binEdgesx, std::vector<double> &binEdgesy) {

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
    binEdgesx.push_back(totalX);
    binEdgesy.push_back(totalY);
    totalX = totalX + deltaXD;
    totalY = totalY + deltaYD;
  }
}

/*
 * saveDiffHistogram = generate a histogram of the difference between two variables
 */
void saveDiffHistogram(std::vector<double> delta,
                       std::string particleName, std::string variableNames,
                       std::string appendName = "", bool saveAll = false, std::string xAxisName = "") {
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  std::string pathName = saveAllInFolder("logLikelihood", variableNames, saveAll);
  std::string fileName = pathName + "DifferenceHistograms" + particleName + appendName;
  fileName = fileName + variableNames;

  std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>(("Difference Scatter Histogram" + appendName).c_str()));
  std::shared_ptr<TH1F> h1(std::make_shared<TH1F>(("Difference Scatter Histogram" + appendName).c_str(), "Difference Scatter Histogram", 100, -5, -5));
  std::vector<double> weights1(delta.size(), 1);

  c1->cd();
  c1->SetGrid();
  h1->FillN(delta.size(), &delta[0], &weights1[0]);
  h1->Draw();
  h1->GetYaxis()->SetTitle("Count");
  h1->GetXaxis()->SetTitle(("Difference in the " + xAxisName).c_str());
  h1->SetTitle("Difference Scatter Plot");

  c1->Modified();
  c1->Update();
  saveCanvas(fileName, c1);
}

/*
 * saveXYScatterPlot = creates a scatter plot between the x and y variable
 * this function is used to generate graphs of
 * "true mean vs delta between generated and interpolated"
*/
void saveXYScatterPlot(std::vector<double> x, std::vector<double> y,
                       std::string particleName, std::string variableNames,
                       bool saveAll = false, std::string xAxisName = "",
                       std::string yAxisName = "") {
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  std::string pathName = saveAllInFolder("logLikelihood", variableNames, saveAll);
  std::string fileName = pathName + "DifferenceScatter" + particleName;
  fileName = fileName + variableNames;

  std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>("Difference Scatter Plot"));
  std::shared_ptr<TGraphErrors> g1(std::make_shared<TGraphErrors>(x.size(), &x[0], &y[0]));
  c1->cd();
  g1->Draw("AP");
  c1->SetGrid();
  g1->GetXaxis()->SetTitle((xAxisName).c_str());
  g1->GetYaxis()->SetTitle((yAxisName).c_str());
  g1->SetTitle("Difference Scatter Plot");

  c1->Modified();
  c1->Update();
  saveCanvas(fileName, c1);
}

/*
 * saveMeanInterpScatterPlot = save the mean number of photons as a function of the a
 * selected variable
 */
void saveMeanInterpScatterPlot(std::vector<double> variableList,
                               std::vector<std::vector<double>> meanVariation,
                               std::vector<bool> variableSwitch,
                               std::string variableIdentifiers, bool saveEach = false,
                               bool saveAllInSingleFolder = false) {
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  std::string variableName;
  double variableValue;
  getVariableFromSwitch(variableList, variableSwitch, variableValue, variableName);

  std::string pathName = saveAllInFolder("VariableMeanScatterPlot", variableIdentifiers, saveAllInSingleFolder);
  std::string fileName = "DifferenceScatter" + variableName;
  fileName = fileName + variableIdentifiers;

  std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>("Difference Scatter Plot"));
  std::shared_ptr<TMultiGraph> mg(std::make_shared<TMultiGraph>());
  std::vector<std::shared_ptr<TGraphErrors>> gs;

  std::vector<int> color;
  color.push_back(kRed);
  color.push_back(kBlue);
  color.push_back(kYellow);
  color.push_back(kTeal);
  color.push_back(kMagenta);
  color.push_back(kAzure);
  color.push_back(kCyan);
  color.push_back(kViolet);
  color.push_back(kOrange);
  color.push_back(kGreen);

  int colorCount = 0;
  std::vector<std::string> paths;
  for (unsigned long i = 0; i < meanVariation.size(); i++) {
    for (unsigned long j = 0; j < meanVariation.at(i).size(); j++) {
      if (i == 0) {
        if (saveEach) {
          paths.push_back(pathName + std::to_string(j) + "/");
          mkdir((paths.back()).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        }

        gs.push_back(std::make_shared<TGraphErrors>(meanVariation.size()));
        gs.at(j)->SetTitle(std::to_string(j).c_str());
      }
      gs.at(j)->SetPoint(i, variableList.at(i), meanVariation[i][j]);

      if (i == meanVariation.size() - 1) {
        gs.at(j)->SetMarkerColor(color.at(colorCount));
        gs.at(j)->SetLineColor(color.at(colorCount));
        if (saveEach) {
          std::shared_ptr<TCanvas> ctemp(std::make_shared<TCanvas>(std::to_string(j).c_str()));
          gs.at(j)->SetLineColor(kBlack);
          gs.at(j)->Draw("APL");
          std::string fileNamePerPixel = paths.at(j) + fileName + "_Pixel" + std::to_string(j);

          ctemp->Modified();
          ctemp->Update();
          saveCanvas(fileNamePerPixel, ctemp);
        }

        mg->Add(gs.at(j).get());
        colorCount++;
      }
      if (colorCount == color.size()) {
        colorCount = 0;
      }
    }
  }

  mg->Draw("APL");
  mg->GetXaxis()->SetTitle(variableName.c_str());
  mg->GetYaxis()->SetTitle("MeanValue");
  mg->SetTitle("Difference Scatter Plot");

  c1->Modified();
  c1->Update();
  saveCanvas(pathName + fileName, c1);
}

/*
 * organizeMeans = puts generated means in the same order as the interpolated means
 */
std::vector<double> organizeMeans(std::vector<double> mean, std::shared_ptr<PMT> pmt) {
  int pixelN = pmt->getPixelCount();
  std::vector<double> binEdgesx;
  binEdgesx.reserve(sqrt(pixelN));
  std::vector<double> binEdgesy;
  binEdgesy.reserve(sqrt(pixelN));
  getBinEdges(pmt, binEdgesx, binEdgesy);

  std::shared_ptr<TH2F> h2(std::make_shared<TH2F>("Test Mean",
                                                  "Test Mean",
                                                  binEdgesx.size() - 1,
                                                  &binEdgesx[0],
                                                  binEdgesy.size() - 1,
                                                  &binEdgesy[0]));

  int meanCount = 0;
  for (unsigned long j = 1; j < binEdgesy.size(); j++) {
    for (unsigned long k = 1; k < binEdgesx.size(); k++) {
      h2->SetBinContent(k, j, mean.at(meanCount));
      meanCount++;
    }
  }

  std::vector<double> tempMean;
  tempMean.reserve(pixelN);
  //meanCount = 0;
  for (unsigned long j = binEdgesy.size() - 1; j > 0; j--) {
    for (unsigned long k = 1; k < binEdgesx.size(); k++) {
      tempMean.push_back(h2->GetBinContent(k, j));
    }
  }
  return tempMean;
}

/*
 * saveHistogramOfMean = save 2D histogram as a mean 
 */
void saveHistogramOfMean(std::vector<double> mean,
                         std::shared_ptr<PMT> pmt,
                         std::string variableName,
                         std::string AppendName = "",
                         bool saveAllInSingleFolder = false,
                         bool printNumbers = false) {
  std::string pathName = saveAllInFolder("Mean2DHistogram", variableName, saveAllInSingleFolder);
  std::string fileName = pathName + "Mean2DHistogram_" + AppendName + "_";
  fileName = fileName + variableName;

  std::cout << fileName << std::endl;

  int pixelN = pmt->getPixelCount();
  int pixelSide = sqrt(pixelN);
  std::vector<double> binEdgesx;
  binEdgesx.reserve(sqrt(pixelN));
  std::vector<double> binEdgesy;
  binEdgesy.reserve(sqrt(pixelN));
  getBinEdges(pmt, binEdgesx, binEdgesy);
  std::cout << pixelN << "    ";
  std::cout << binEdgesx.size() << "    ";
  std::cout << binEdgesy.size() << "    ";
  std::cout << mean.size() << std::endl;
  std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>(fileName.c_str(), "", 640, 480));
  std::shared_ptr<TH2F> h2(std::make_shared<TH2F>("2D Histogram",
                                                  "2D Histogram of Detector Face (Each bin is Sum)",
                                                  binEdgesx.size() - 1,
                                                  &binEdgesx[0],
                                                  binEdgesy.size() - 1,
                                                  &binEdgesy[0]));

  std::vector<double> meanCountList;
  int meanCount = 0;

  for (int j = 1; j < binEdgesy.size(); j++) {
    for (int k = 1; k < binEdgesx.size(); k++) {
      double flippedPosition = (pixelN - pixelSide) + 2 * (meanCount - pixelSide * TMath::Floor(meanCount / pixelSide)) - meanCount;
      h2->SetBinContent(k, j, mean.at(flippedPosition));
      meanCountList.push_back(flippedPosition);
      meanCount++;
    }
  }

  h2->GetXaxis()->SetRange(binEdgesx.front(), binEdgesx.back());
  h2->GetYaxis()->SetRange(binEdgesy.front(), binEdgesy.back());
  h2->SetStats();
  h2->GetXaxis()->SetNdivisions(-16);
  h2->GetYaxis()->SetNdivisions(-16);
  h2->GetXaxis()->SetLabelSize(0.015);
  h2->GetYaxis()->SetLabelSize(0.015);
  h2->GetXaxis()->SetTitle("X Position");
  h2->GetYaxis()->SetTitle("Y Position");
  c1->SetGrid();
  h2->Draw("COLZ2");

  meanCount = 0;
  //meanCountList = organizeMeans(meanCountList, pmt);
  if (printNumbers) {
    for (int j = 1; j < binEdgesy.size(); j++) {
      double y = h2->GetYaxis()->GetBinCenter(j);
      for (int k = 1; k < binEdgesx.size(); k++) {
        double x = h2->GetXaxis()->GetBinCenter(k);
        //std::cout << "x:" << x << "   y:" << y << std::endl;
        std::shared_ptr<TText> t = std::make_shared<TText>();
        c1->cd();
        t->SetTextAngle(0);
        t->SetTextSize(0.02);
        t->SetTextAlign(33);
        std::stringstream ss;
        ss << meanCountList.at(meanCount);
        std::string meanCountString = ss.str();
        t->DrawText(x, y, meanCountString.c_str());
        meanCount++;
      }
    }
  }
  c1->Update();
  c1->Modified();
  saveCanvas(fileName, c1);
}
void plotSingleParticleValues(std::shared_ptr<Detector> detector,
                              std::shared_ptr<InterpolateLikelihood> interp,
                              std::vector<double> variables1,
                              std::vector<std::vector<double>> pixelHits1,
                              std::vector<double> &meanInterp1,
                              std::string variableNames,
                              std::string particle,
                              int interpSwitch = 0,
                              std::string Title = "",
                              std::vector<bool> variableSwitch = std::vector<bool>()) {
  bool saveAllInSingleFolder = true;

  meanInterp1 = interp->interpolationSwitch(variables1, interpSwitch, variableSwitch,
                                            detector->getPMT()->getWidthPixelCount(), detector->getTotalDistance(),
                                            detector->getMaxWidth(), detector->getMaxLength(),
                                            detector->getPMT()->getWidthPixel(), detector->getPMT()->getLengthPixel());

  if (meanInterp1.empty()) {
    return;
  }
  std::string name = "_Interpolated_" + particle + "_Means_" + Title;
  saveHistogramOfMean(meanInterp1, detector->getPMT(), variableNames, name, saveAllInSingleFolder);
  std::vector<double> meanGenerated1 = detector->getMeanPerPixel(pixelHits1);
  std::vector<double> mean1Delta;
  mean1Delta.reserve(meanGenerated1.size());
  for (int pixel = 0; pixel < meanGenerated1.size(); pixel++) {
    mean1Delta.push_back(meanGenerated1[pixel] - meanInterp1[pixel]);
  }
  saveXYScatterPlot(meanGenerated1, mean1Delta, particle + "_GenereatedMean_And_" + Title, variableNames, saveAllInSingleFolder, "Mean_{True}", "Mean_{True} - Mean_{Interp}");
  saveDiffHistogram(mean1Delta, particle, variableNames, particle + "_GenereatedMean_And_" + Title, saveAllInSingleFolder, "Mean_{True} - Mean_{Interp}");
}

/*
 * plotLikelihoodValues = plots all possible comparison graphs between the generated means and the 
 * the interpolated means
 */
void plotComparisonBetweenGeneratedAndInterpolated(std::vector<double> variables1,
                                                   std::vector<double> variables2,
                                                   std::vector<double> generatedLL1,
                                                   std::vector<double> generatedLL2,
                                                   std::shared_ptr<TH1D> generatedLL1Plot,
                                                   std::shared_ptr<TH1D> generatedLL2Plot,
                                                   std::vector<std::vector<double>> pixelHits1,
                                                   std::vector<std::vector<double>> pixelHits2,
                                                   std::shared_ptr<Detector> detector,
                                                   std::shared_ptr<InterpolateLikelihood> interp,
                                                   std::string variableNames,
                                                   int interpSwitch = 0,
                                                   std::string Title = "",
                                                   bool nuissance = false,
                                                   std::vector<bool> variableSwitch = std::vector<bool>()) {
  bool showSubGraphs = false;
  bool saveAllInSingleFolder = true;

  std::vector<double> meanInterp1;
  plotSingleParticleValues(detector, interp, variables1, pixelHits1, meanInterp1, variableNames, "Pion", interpSwitch, Title, variableSwitch);

  std::vector<double> meanInterp2;
  std::vector<double> meanGenerated2;
  plotSingleParticleValues(detector, interp, variables2, pixelHits2, meanInterp2, variableNames, "Kaon", interpSwitch, Title, variableSwitch);

  if (meanInterp1.empty() || meanInterp2.empty()) {
    return;
  }

  std::vector<double> ll1Interp;
  std::vector<double> ll2Interp;
  if (nuissance == false) {
    ll1Interp = detector->getLogLikelihood(pixelHits1, meanInterp1, meanInterp2);
    ll2Interp = detector->getLogLikelihood(pixelHits2, meanInterp1, meanInterp2);
  } else {
    ll1Interp = detector->getLogLikelihoodWithNuisanceMinuit(pixelHits1, meanInterp1, meanInterp2);
    ll2Interp = detector->getLogLikelihoodWithNuisanceMinuit(pixelHits2, meanInterp1, meanInterp2);
  }

  std::shared_ptr<TH1D> ll_1_hist_Interp = getHistogramOfLL(ll1Interp);
  std::shared_ptr<TH1D> ll_2_hist_Interp = getHistogramOfLL(ll2Interp);

  std::vector<double> ll_delta_1;
  ll_delta_1.reserve(generatedLL1.size());
  std::vector<double> ll_delta_2;
  ll_delta_2.reserve(generatedLL2.size());
  for (int event = 0; event < generatedLL1.size(); event++) {
    ll_delta_1.push_back(generatedLL1[event] - ll1Interp[event]);
    ll_delta_2.push_back(generatedLL2[event] - ll2Interp[event]);
  }

  //Both Histograms
  plotBothLLHistograms(ll_1_hist_Interp, ll_2_hist_Interp,
                       showSubGraphs, "Pion", "Kaon",
                       variableNames, Title,
                       saveAllInSingleFolder);

  //Comparing 1 histogramsgeneratedLL1Plot
  plotBothLLHistograms(generatedLL1Plot, ll_1_hist_Interp,
                       showSubGraphs, "Pion", "Pion Interpolated " + Title,
                       variableNames, "_Comparing_Pions_",
                       saveAllInSingleFolder);

  //Comparing 2 Histograms
  plotBothLLHistograms(generatedLL2Plot, ll_2_hist_Interp,
                       showSubGraphs, "Kaon", "Kaon Interpolated " + Title,
                       variableNames, "_Comparing_Kaons_",
                       saveAllInSingleFolder);

  saveDiffHistogram(ll_delta_1, "Pion", variableNames, "Events_Interp_" + Title, saveAllInSingleFolder, "-2 ln(#frac{L_{#pi}}{L_{k}}) values for Pion Events");
  saveDiffHistogram(ll_delta_2, "Kaon", variableNames, "Events_Interp_" + Title, saveAllInSingleFolder, "-2 ln(#frac{L_{#pi}}{L_{k}}) values for Kaon Events");
  return;
}

void SaveMeanQuick(bool betaBool, std::shared_ptr<Beam> beam, std::vector<double> mean) {
  int VARIABLE_SIZE = 5;

  std::string pathName = "./Textfiles/";
  mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  std::string fileName = pathName;
  fileName += getFileName(betaBool, beam);
  fileName += ".bin";
  FILE *f = fopen(fileName.c_str(), "wb");
  assert(f); //does f exist

  double v[VARIABLE_SIZE];
  if (betaBool) {
    v[0] = beam->getBeta();
  } else {
    v[0] = beam->getMomentum();
  }
  v[1] = beam->getPosition()->getX();
  v[2] = beam->getPosition()->getY();
  v[3] = beam->getDirection()->getTheta();
  v[4] = beam->getDirection()->getPhi();
  fwrite(v, sizeof(v[0]), VARIABLE_SIZE, f);
  fwrite("\n", sizeof(char), 1, f);
  //for (int pixel = 0; pixel < mean.size(); pixel++) {
  fwrite(&mean[0], sizeof(mean), mean.size(), f);
  //fwrite("\n", sizeof(char), 1, f);
  //}
}

/*
 * saveAllPixelsAndSums = plot the likelihood of each event for each pixel
 */
void saveAllPixelsAndSums(std::vector<std::vector<double>> llPerPixel, std::vector<double> mean1, std::vector<double> mean2, std::string name,
                          std::vector<double> weight1 = std::vector<double>(), std::vector<double> weight2 = std::vector<double>()) {
  bool first = true;
  std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>(name.c_str(), name.c_str(), 1200, 1200));
  std::vector<std::shared_ptr<TH1D>> hs;
  std::shared_ptr<TH1D> sumHisto = std::make_shared<TH1D>((name + "Sum").c_str(), (name + "Sum").c_str(), 100, 0, 0);
  std::vector<double> sumPixels;
  for (int i = 0; i < llPerPixel.size(); i++) {
    double sumPixels = 0;
    for (int pixel = 0; pixel < llPerPixel[i].size(); pixel++) {
      if (i == 0) {
        std::string Histname = name + "perPixel";
        Histname = Histname + std::to_string(pixel);
        hs.push_back(std::make_shared<TH1D>(Histname.c_str(), Histname.c_str(), 50, 0, 0));
      }
      sumPixels += llPerPixel[i][pixel];
      hs[pixel]->Fill(llPerPixel[i][pixel]);
      if (i == llPerPixel.size() - 1 && hs[pixel]->GetMean() != 0 && hs[pixel]->GetStdDev() != 0) {
        hs[pixel]->Draw();
        std::string titleName1;
        titleName1 += "1: " + std::to_string(mean1[pixel]) + " +/- " + std::to_string(sqrt(mean1[pixel]) / sqrt(1E4));
        titleName1 += "   LLVal Hit:" + std::to_string(-2 * log(1 - exp(-mean1[pixel])));
        titleName1 += "   LLVal NoHit:" + std::to_string(-2 * log(exp(-mean1[pixel])));
        std::string titleName2;
        titleName2 += "2: " + std::to_string(mean2[pixel]) + " +/- " + std::to_string(sqrt(mean2[pixel]) / sqrt(1E4));
        titleName2 += "   LLVal Hit:" + std::to_string(-2 * log(1 - exp(-mean2[pixel])));
        titleName2 += "   LLVal noHit:" + std::to_string(-2 * log(exp(-mean2[pixel])));

        if (!weight1.empty() && !weight2.empty()) {
          titleName1 += "    W:" + std::to_string(weight1[pixel]);
          titleName2 += "    W:" + std::to_string(weight2[pixel]);
        }

        TPaveText *ps = new TPaveText(.1, .85, .9, 1, "NDC");
        ps->SetName("mystats");
        ps->AddText((titleName1 + "\n\r").c_str());
        ps->AddText(titleName2.c_str());
        ps->Draw();

        c1->Update();
        c1->Modified();
        std::string fileName = "./Images/" + name;
        fileName = fileName + ".pdf";
        if (first) {
          fileName = fileName + "(";
          first = false;
        }
        c1->Print(fileName.c_str());
      }
    }
    sumHisto->Fill(sumPixels);
  }
  sumHisto->Draw();
  std::string fileName = "./Images/" + name;
  fileName = fileName + ".pdf)";
  c1->Print(fileName.c_str());
}

//Organize Pixel Events
std::vector<std::vector<double>> quickOrg(std::vector<std::vector<double>> pixelHits) {
  std::vector<std::vector<double>> tempPH;
  tempPH.reserve(pixelHits.size());
  for (unsigned long i = 0; i < pixelHits.size(); i++) {
    std::vector<double> tempPHevent(pixelHits.at(i).size(), 0.0);
    for (unsigned long j = 0; j < sqrt(pixelHits.at(i).size()); j++) {
      for (unsigned long k = 0; k < sqrt(pixelHits.at(i).size()); k++) {
        tempPHevent.at((sqrt(pixelHits.at(i).size()) - j) * (sqrt(pixelHits.at(i).size())) - (sqrt(pixelHits.at(i).size()) - k)) = pixelHits[i][j * sqrt(pixelHits.at(i).size()) + k];
      }
    }
    tempPH.push_back(tempPHevent);
  }
  return tempPH;
}

//Seperation Code wrapper
double getSeperationValueFromPixelHitsAndMean(std::shared_ptr<Detector> detector,
                                              std::vector<std::vector<double>> pH1, std::vector<std::vector<double>> pH2,
                                              std::vector<double> mean1, std::vector<double> mean2,
                                              std::string title, bool nuisance = false) {

  std::vector<double> ll1;
  std::vector<double> ll2;

  if (nuisance) {
    //std::cout << "ROOT FINDING" << std::endl;
    //ll1 = detector->getLogLikelihoodWithNuisance(pH1, mean1, mean2, false);
    //ll2 = detector->getLogLikelihoodWithNuisance(pH2, mean1, mean2, false);
    //std::cout << "ll1: " << ll1[0] << "    ";
    //std::cout << "ll2: " << ll2[0] << std::endl;

    std::cout << "Minimization" << std::endl;
    ll1 = detector->getLogLikelihoodWithNuisanceMinuit(pH1, mean1, mean2);
    ll2 = detector->getLogLikelihoodWithNuisanceMinuit(pH2, mean1, mean2);

  } else {
    ll1 = detector->getLogLikelihood(pH1, mean1, mean2);
    ll2 = detector->getLogLikelihood(pH2, mean1, mean2);
  }

  std::shared_ptr<TH1D> ll_1_hist = getHistogramOfLL(ll1, "1" + title);
  std::shared_ptr<TH1D> ll_2_hist = getHistogramOfLL(ll2, "2" + title);

  double area1 = ll_1_hist->Integral();
  double area2 = ll_2_hist->Integral();
  ll_1_hist->Scale(1 / (area1));
  ll_2_hist->Scale(1 / (area2));
  double sep = getSeperationValue(ll_1_hist, ll_2_hist);

  return sep;
}

double calcStepSize(double start, double end, double step) {
  return (end - start) / (step - 1);
}

void setVariablesFromJobFiles(std::string fileName,
                              std::vector<double> *xVariables,
                              std::vector<double> *yVariables,
                              std::vector<double> *thetaVariables,
                              std::vector<double> *phiVariables,
                              std::vector<double> *betaMomentumVariables) {
  //std::cout << fileName << std::endl;
  std::ifstream infile(fileName.c_str());
  if (infile.is_open()) {
    std::cout << "TEST" << std::endl;
    while (infile.good()) {
      std::string line;
      getline(infile, line);
      //std::cout << line << std::endl;
      if (line.find("-") == std::string::npos) {
        std::istringstream iss(line);
        std::string variableName;
        double min, max, step;
        iss >> variableName >> min >> max >> step;

        if ("X" == variableName) {
          //std::cout << "X" << std::endl;
          xVariables->at(0) = min;
          xVariables->at(1) = max;
          xVariables->at(2) = calcStepSize(min, max, step);
          std::cout << variableName << "   ";
          std::cout << min << "   ";
          std::cout << max << "   ";
          std::cout << calcStepSize(min, max, step) << std::endl;
        } else if ("Y" == variableName) {
          std::cout << "Y" << std::endl;
          yVariables->at(0) = min;
          yVariables->at(1) = max;
          yVariables->at(2) = calcStepSize(min, max, step);
          std::cout << variableName << "   ";
          std::cout << min << "   ";
          std::cout << max << "   ";
          std::cout << calcStepSize(min, max, step) << std::endl;

        } else if ("Theta" == variableName) {
          std::cout << "Theta" << std::endl;
          thetaVariables->at(0) = min;
          thetaVariables->at(1) = max;
          thetaVariables->at(2) = calcStepSize(min, max, step);
          std::cout << variableName << "   ";
          std::cout << min << "   ";
          std::cout << max << "   ";
          std::cout << calcStepSize(min, max, step) << std::endl;
        } else if ("Phi" == variableName) {
          std::cout << "Phi" << std::endl;
          phiVariables->at(0) = min;
          phiVariables->at(1) = max;
          phiVariables->at(2) = calcStepSize(min, max, step);
          std::cout << variableName << "   ";
          std::cout << min << "   ";
          std::cout << max << "   ";
          std::cout << calcStepSize(min, max, step) << std::endl;
        } else if ("Beta" == variableName) {
          std::cout << "Beta" << std::endl;
          betaMomentumVariables->at(0) = min;
          betaMomentumVariables->at(1) = max;
          betaMomentumVariables->at(2) = calcStepSize(min, max, step);
          std::cout << variableName << "   ";
          std::cout << min << "   ";
          std::cout << max << "   ";
          std::cout << calcStepSize(min, max, step) << std::endl;
        }
	else if ("Momentum" == variableName) {
          std::cout << "Momentum" << std::endl;
          betaMomentumVariables->at(0) = min;
          betaMomentumVariables->at(1) = max;
          betaMomentumVariables->at(2) = calcStepSize(min, max, step);
          std::cout << variableName << "   ";
          std::cout << min << "   ";
          std::cout << max << "   ";
          std::cout << calcStepSize(min, max, step) << std::endl;
	
	}
      }
    }
  }
}

void readRootFile(std::string fileName,
                  std::vector<double> *xVariables,
                  std::vector<double> *yVariables,
                  std::vector<double> *thetaVariables,
                  std::vector<double> *phiVariables,
                  std::vector<double> *betaMomentumVariables,
                  std::vector<double> *particleType) {
  std::shared_ptr<TFile> f = std::make_shared<TFile>(fileName.c_str());
  std::shared_ptr<TTree> t((TTree *)f->Get("h1000"));
  float mom, x, y, z, vx, vy, vz;
  float pos[3];
  float vec[3];
  int gpid[50];
  int ng;
  t->SetBranchAddress("mom", &mom);
  t->SetBranchAddress("pos", pos);
  t->SetBranchAddress("vec", vec);
  t->SetBranchAddress("gpid", gpid);
  t->SetBranchAddress("ng", &ng);
  Int_t nentries = (Int_t)t->GetEntries();
  for (Int_t i = 0; i < nentries; i++) {
    t->GetEntry(i);
    //std::cout << mom << "  ";
    //std::cout << pos[0] << "  ";
    //std::cout << pos[1] << "  ";
    //std::cout << pos[2] << "  ";
    //std::cout << vec[0] << "  ";
    //std::cout << vec[1] << "  ";
    //std::cout << vec[2] << std::endl;
    pos[0] = pos[0] / 100;
    pos[1] = pos[1] / 100;
    pos[2] = pos[2] / 100;

    vec[0] = vec[0] / 100;
    vec[1] = vec[1] / 100;
    vec[2] = vec[2] / 100;

    if (vec[2] >= 0) {
      betaMomentumVariables->push_back(mom * 1000);                         //MeV
      xVariables->push_back(pos[0]);                                        //m
      yVariables->push_back(pos[1]);                                        //m
      double r = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]); //m^3
      thetaVariables->push_back(acos(vec[2] / r));                          //rad
      phiVariables->push_back(atan(vec[1] / vec[0]));                       //rad
      particleType->push_back(gpid[ng]);
    }

    //std::cout << betaMomentumVariables->back() << "   ";
    //std::cout << xVariables->back() << "   ";
    //std::cout << yVariables->back() << "   ";
    //std::cout << thetaVariables->back() << "   ";
    //std::cout << phiVariables->back() << std::endl;
  }
}

void saveSeperationParameter(std::string filename,
                             bool betaBool, std::shared_ptr<Beam> beam,
                             double genSep, double genNucSep, double intSep, double intNucSep) {
  int VARIABLE_SIZE = 5;
  int SEPERATION_SIZE = 4;
  std::string pathName = "./Textfiles/";
  mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  std::string pathFileName = pathName;
  pathFileName += filename;
  pathFileName += ".txt";
  FILE *f = fopen(pathFileName.c_str(), "at");
  assert(f); //does f exist

  double v[VARIABLE_SIZE];
  v[0] = beam->getMomentum();
  v[1] = beam->getPosition()->getX();
  v[2] = beam->getPosition()->getY();
  v[3] = beam->getDirection()->getTheta();
  v[4] = beam->getDirection()->getPhi();
  //fwrite(v, sizeof(v[0]), VARIABLE_SIZE, f);
  for (int variable = 0; variable < VARIABLE_SIZE; variable++) {
    //std::cout << v[variable] << " ";
    fprintf(f, "%g ", v[variable]);
  }
  double s[SEPERATION_SIZE];
  s[0] = genSep;
  s[1] = genNucSep;
  s[2] = intSep;
  s[3] = intNucSep;
  for (int sperationVariable = 0; sperationVariable < SEPERATION_SIZE; sperationVariable++) {
    //std::cout << v[variable] << " ";
    fprintf(f, "%g ", s[sperationVariable]);
  }
  fwrite("\n", sizeof(char), 1, f);
  fclose(f);
}

void SaveSepQuick(std::string fileName, bool betaBool, std::shared_ptr<Beam> beam, std::vector<double> mean) {
  int VARIABLE_SIZE = 5;

  std::string pathName = "./Textfiles/";
  mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  std::string pathFileName = pathName;
  pathFileName += fileName;
  pathFileName += ".txt";
  FILE *f = fopen(pathFileName.c_str(), "a");
  assert(f); //does f exist

  double v[VARIABLE_SIZE];
  if (betaBool) {
    v[0] = beam->getBeta();
  } else {
    v[0] = beam->getMomentum();
  }
  v[1] = beam->getPosition()->getX();
  v[2] = beam->getPosition()->getY();
  v[3] = beam->getDirection()->getTheta();
  v[4] = beam->getDirection()->getPhi();
  fwrite(v, sizeof(v[0]), VARIABLE_SIZE, f);
  fwrite("\n", sizeof(char), 1, f);
  fwrite(&mean[0], sizeof(mean), mean.size(), f);
}

void readParticleTable(std::vector<double> *particleMasses, std::vector<double> *wSigmas) {
  std::ifstream infile("../ParticleTable.txt");
  while (infile.good()) {
    std::string line;
    getline(infile, line);
    if (line.find("#") == std::string::npos) {
      std::istringstream iss(line);
      std::string particleName;
      int pID, option, charge;
      double mass, lifetime, sig;
      iss >> particleName >> pID >> option >> mass >> charge >> lifetime >> sig;
      particleMasses->push_back(1000 * mass); //GeV/c^2
      if (sig == 0) {
        sig = sig + 1;
      }
      wSigmas->push_back(sig);
    }
  }
}

void saveEachPixelMeanVsParameter(std::string filename,
                                  std::vector<std::vector<double>> mean,
                                  std::vector<std::vector<double>> delta,
                                  std::vector<double> p) {

  std::string pathName = "./Textfiles/";
  mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  std::string pathFileName = pathName;
  pathFileName += filename;
  pathFileName += ".txt";
  FILE *f = fopen(pathFileName.c_str(), "at");
  assert(f); //does f exist
  for (int parameterCount = 0; parameterCount < mean.size(); parameterCount++) {
    fprintf(f, "%g ", p[parameterCount]);
    for (int pixel = 0; pixel < mean[parameterCount].size(); pixel++) {
      fprintf(f, "%g ", mean[parameterCount][pixel]);
    }
    for (int pixel = 0; pixel < delta[parameterCount].size(); pixel++) {
      fprintf(f, "%g ", delta[parameterCount][pixel]);
    }
    fwrite("\n", sizeof(char), 1, f);
  }
  fclose(f);
}


void saveParameterAndLikelihood(std::string filename,
                                bool betaBool, std::shared_ptr<Beam> beam,
                                std::vector<std::vector<double>> vMin,
                                std::vector<std::vector<double>> vMin_err,
                                std::vector<std::vector<double>> vWar,
                                std::vector<double> ll) {
  int VARIABLE_SIZE = 5;
  int MINIMIZATION_SIZE = 11;
  std::string pathName = "./Textfiles/";
  mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  std::string pathFileName = pathName;
  pathFileName += filename;
  pathFileName += ".txt";
  FILE *f = fopen(pathFileName.c_str(), "at");
  assert(f); //does f exist

  double v[VARIABLE_SIZE];
  v[0] = beam->getMomentum();
  v[1] = beam->getPosition()->getX();
  v[2] = beam->getPosition()->getY();
  v[3] = beam->getDirection()->getTheta();
  v[4] = beam->getDirection()->getPhi();
  //fwrite(v, sizeof(v[0]), VARIABLE_SIZE, f);
  for (int variable = 0; variable < VARIABLE_SIZE; variable++) {
    //std::cout << v[variable] << " ";
    fprintf(f, "%g ", v[variable]);
  }
  fwrite("\n", sizeof(char), 1, f);

  for (int event = 0; event < vMin.size(); event++) {
    fprintf(f, "%d ", event);
    for (int parameter = 0; parameter < vMin[event].size(); parameter++) {
      fprintf(f, "%g ", vMin[event][parameter]);
      fprintf(f, "%g ", vMin_err[event][parameter]);
    }
    for (int warnings = 0; warnings < vWar[event].size(); warnings++) {
      fprintf(f, "%g ", vWar[event][warnings]);
    }
    fprintf(f, "%g ", ll[event]);
    fwrite("\n", sizeof(char), 1, f);
  }
  fclose(f);
}

} // namespace graphingAndAnalysisCherenkovPID
