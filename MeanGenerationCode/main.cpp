
// C++
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <numeric>
#include <regex>
#include <string>
#include <sys/stat.h>
#include <vector>
// Boost
// include "/usr/local/include/boost/algorithm/string.hpp"

// Root
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TEllipse.h"
#include "TGFrame.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TVirtualPad.h"

// Classes
#include "aerogelClass.h"
#include "detectedEventClass.h"
#include "detectorClass.h"
#include "directionClass.h"
#include "pmtClass.h"
#include "positionClass.h"
#include "randomEventClass.h"

void checkRandomOutput(
		std::vector<std::shared_ptr<RandomEvent>> eventsForAerogels) {
	// Checks the Randomly Generated Events by outputing their variables to the
	// terminal
	for (int i = 0; i < eventsForAerogels.size(); i++) {
		std::shared_ptr<RandomEvent> events = eventsForAerogels.at(i);
		// std::cout << "Aerogel:  " << i << std::endl;
		std::vector<double> radius = events->getRadius();
		std::vector<double> theta = events->getPhi();
		std::vector<double> wvl = events->getWavelength();
		std::vector<double> energy = events->getEnergy();
		std::vector<double> qE = events->getQuantumEfficiency();
		// std::cout << "j"    << "        ";
		// std::cout << "rad"  << "        ";
		// std::cout << "theta"<< "        ";
		// std::cout << "wvl"  << "        ";
		// std::cout << "efficieny" << "     ";
		// std::cout << "energy" << std::endl;
		for (int j = 0; j < events->getN(); j++) {
			// std::cout << j << "    ";
			// std::cout << radius.at(j) << "    ";
			// std::cout << theta.at(j)  << "    ";
			// std::cout << wvl.at(j)    << "    ";
			// std::cout << qE.at(j)     << "    ";
			// std::cout << energy.at(j) << std::endl;
		}
	}
}

void checkDetectedOutput(
		std::vector<std::shared_ptr<RandomEvent>> eventsForAerogels) {
	// Checkes the Detected Generated Events by outputing their variables to the
	// terminal
	for (int i = 0; i < eventsForAerogels.size(); i++) {
		std::shared_ptr<RandomEvent> events = eventsForAerogels.at(i);
		std::shared_ptr<DetectedEvent> DetectedEvent = events->getDetectedHits();
		// std::cout << "Aerogel:  " << i << std::endl;
		std::vector<double> radius = DetectedEvent->getRadius();
		std::vector<double> theta = DetectedEvent->getPhi();
		std::vector<double> wvl = DetectedEvent->getWavelength();
		std::vector<double> energy = DetectedEvent->getEnergy();
		// std::vector<double> qE = DetectedEvent->getQuantumEfficiency();
		// std::cout << "j"    << "        ";
		// std::cout << "rad"  << "        ";
		// std::cout << "theta"<< "        ";
		// std::cout << "wvl"  << "        ";
		////std::cout << "efficieny" << "     ";
		// std::cout << "energy" << std::endl;
		for (int j = 0; j < DetectedEvent->getN(); j++) {
			// std::cout << j << "    ";
			// std::cout << radius.at(j) << "    ";
			// std::cout << theta.at(j)  << "    ";
			// std::cout << wvl.at(j)    << "    ";
			////std::cout << qE.at(j)     << "    ";
			// std::cout << energy.at(j) << std::endl;
		}
	}
}

void plotHistogramOfN(std::vector<std::vector<std::shared_ptr<DetectedEvent>>>
		eventsPerRunPerAerogels,
		std::string particle) {
	// Plot the Histogram of the detected Events on the face of the PMT
	// This is not per pixel this is the total detected photons
	// Input: 2D Vector of Detected Events

	// Create Path Name and File names
	std::string pathName = "./Images/TotalDetectedPhotons" + particle;
	std::string namePng = pathName + ".png";
	std::string nameRoot = pathName + ".root";
	std::string namePdf = pathName + ".pdf";

	std::vector< double> detectedPhotons;
	detectedPhotons.reserve(eventsPerRunPerAerogels.size());
	std::vector< double> detectedPhotonsAerogel1;
	detectedPhotonsAerogel1.reserve(eventsPerRunPerAerogels.size());
	std::vector< double> detectedPhotonsAerogel2;
	detectedPhotonsAerogel2.reserve(eventsPerRunPerAerogels.size());

	for (int i = 0; i < eventsPerRunPerAerogels.size(); i++) {
		std::vector<std::shared_ptr<DetectedEvent>> eventsForAerogels =
			eventsPerRunPerAerogels.at(i);
		double totalN = 0;
		for (int j = 0; j < eventsForAerogels.size(); j++) {
			std::shared_ptr<DetectedEvent> dE = eventsForAerogels.at(j);
			if (j == 0) {
				detectedPhotonsAerogel1.push_back(dE->getN());
			} else if (j == 1) {
				detectedPhotonsAerogel2.push_back(dE->getN());
			} else {
				std::cout << "SOMETHING HAS GONE TERRIBLY WRONG" << std::endl;
			}
			totalN = totalN + dE->getN();
		}
		detectedPhotons.push_back(totalN);
	}

	// Note Weights are generated of TH1 FILLN Function
	std::vector< double> weights(detectedPhotons.size(), 1.0);

	// Histogram of Detected Events on the face of the PMT
	std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>(pathName.c_str()));
	std::shared_ptr<TH1D> h1(std::make_shared<TH1D>("h1", "Histogram", 10, 0, 0));
	c1->SetGrid();
	h1->FillN(detectedPhotons.size(), &(detectedPhotons[0]), &(weights[0]));
	h1->SetStats();
	h1->SetFillColor(38);
	h1->SetXTitle("Number of Detected Photons");
	h1->SetYTitle("Counts");
	h1->Scale(1);
	h1->Draw();
	c1->Modified();
	c1->SaveAs(namePng.c_str());
	c1->SaveAs(nameRoot.c_str());
	c1->SaveAs(namePdf.c_str());

	// Histogram of Detected Events for the first  aerogel
	// Create Path Name and File names
	std::string pathNameAero1 = "./Images/TotalDetectedPhotonsAero1" + particle;
	namePng = pathNameAero1 + ".png";
	nameRoot = pathNameAero1 + ".root";
	namePdf = pathNameAero1 + ".pdf";

	std::shared_ptr<TCanvas> c2(std::make_shared<TCanvas>(pathNameAero1.c_str()));
	std::shared_ptr<TH1D> h2(
			std::make_shared<TH1D>("h2", "Histogram of Aerogel1", 10, 0, 0));
	c2->SetGrid();
	h2->FillN(detectedPhotonsAerogel1.size(), &(detectedPhotonsAerogel1[0]),
			&(weights[0]));
	h2->SetStats();
	h2->SetFillColor(38);
	h2->SetXTitle("Number of Detected Photons");
	h2->SetYTitle("Counts");
	h2->Scale(1);
	h2->Draw();
	c2->Modified();
	c2->SaveAs(namePng.c_str());
	c2->SaveAs(nameRoot.c_str());
	c2->SaveAs(namePdf.c_str());

	// Histogram of Detected Events for the second aerogel
	// Create Path Name and File names
	std::string pathNameAero2 = "./Images/TotalDetectedPhotonsAero2" + particle;
	namePng = pathNameAero2 + ".png";
	nameRoot = pathNameAero2 + ".root";
	namePdf = pathNameAero2 + ".pdf";

	std::shared_ptr<TCanvas> c3(std::make_shared<TCanvas>(pathNameAero2.c_str()));
	std::shared_ptr<TH1D> h3(
			std::make_shared<TH1D>("h3", "Histogram of Aerogel1", 10, 0, 0));
	c3->SetGrid();
	h3->FillN(detectedPhotonsAerogel2.size(), &(detectedPhotonsAerogel2[0]),
			&(weights[0]));
	h3->SetStats();
	h3->SetFillColor(38);
	h3->SetXTitle("Number of Detected Photons");
	h3->SetYTitle("Counts");
	h3->Scale(1);
	h3->Draw();
	c3->Modified();
	c3->SaveAs(namePng.c_str());
	c3->SaveAs(nameRoot.c_str());
	c3->SaveAs(namePdf.c_str());
}

void getBinEdges(std::shared_ptr<PMT> pmt, std::vector< double> *binEdgesx,
		std::vector< double> *binEdgesy) {
	// Calculated the binEdges given the PMT value
	// Input:
	// Pointer to vector of doubles for binEdgesx
	// Pointer to vector of doubles for binEdgesy
	double dW = pmt->getWidth();
	double dL = pmt->getLength();

	int pixelN = pmt->getPixelCount();
	int binEdgesN = pmt->getLengthPixelCount() +
		1; // where getLengthPixelCount = getWidthPixelCount

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

std::vector<std::vector<int>>
getPixelHitsPerEvents(std::vector<std::vector<std::shared_ptr<DetectedEvent>>>
		eventsPerRunPerAerogels,
		std::shared_ptr<PMT> pmt) {
	// Calculates the 2D vector of Histogram hits per pixel for each event
	// INPUT:
	// 2D Vector of Detected Hits
	// PMT Object
	double pixelN = pmt->getPixelCount();
	std::vector<double> binEdgesx;
	std::vector<double> binEdgesy;
	getBinEdges(pmt, &binEdgesx, &binEdgesy);

	// double dW = pmt->getWidth();
	// double dL = pmt->getLength();

	// int pixelN = pmt->getPixelCount();
	// int binEdgesN = pmt->getLengthPixelCount()+1; // where getLengthPixelCount
	// = getWidthPixelCount

	// double xD = pmt->getPosition()->getX();
	// double yD = pmt->getPosition()->getY();

	// double deltaXD = dW/(binEdgesN-1);//((xD+dW)-(xD))/N;
	// double deltaYD = dL/(binEdgesN-1);//((yD+dL)-(yD))/N;

	// double totalX = xD;
	// double totalY = yD;
	//////std::cout << xD << "   " << xD+dW << std::endl;
	//////std::cout << yD << "   " << yD+dL << std::endl;

	// for (int i = 0; i < binEdgesN; i++)
	//{
	//////std::cout << i << "    ";
	//////std::cout << totalX << "     ";
	//////std::cout << totalY << std::endl;

	// binEdgesx.push_back(totalX);
	// binEdgesy.push_back(totalY);
	// totalX = totalX + deltaXD;
	// totalY = totalY + deltaYD;
	//}

	// Loops through all the the events (for each aerogel)
	// Then "merges" all the x values/y values for each aerogel (per an event)
	// Then (for that event) calculates the histogram values of x and y
	std::vector<std::vector<int>> hit_on_pixels_perEvent;
	for (int i = 0; i < eventsPerRunPerAerogels.size(); i++) {
		std::vector<std::shared_ptr<DetectedEvent>> eventsForAerogels =
			eventsPerRunPerAerogels.at(i);
		double totalN = 0;
		std::vector<double> x_perEvent;
		std::vector<double> y_perEvent;
		for (int j = 0; j < eventsForAerogels.size(); j++) {
			std::shared_ptr<DetectedEvent> dE = eventsForAerogels.at(j);
			std::vector<double> tempX = dE->getX();
			std::vector<double> tempY = dE->getY();
			x_perEvent.insert(x_perEvent.end(), tempX.begin(), tempX.end());
			y_perEvent.insert(y_perEvent.end(), tempY.begin(), tempY.end());
		}
		std::shared_ptr<TH2F> h2(std::make_shared<TH2F>(
					std::to_string(i).c_str(), std::to_string(i).c_str(),
					binEdgesx.size() - 1, &binEdgesx[0], binEdgesy.size() - 1,
					&binEdgesy[0]));
		std::vector< double> weights(x_perEvent.size(), 1.0);
		h2->FillN(x_perEvent.size(), &(x_perEvent[0]), &(y_perEvent[0]),
				&(weights[0]));
		h2->Draw();

		// Creates a new vector out of each hit per event
		std::vector<int> pixelHits;
		pixelHits.reserve(pixelN);
		for (int j = 1; j < binEdgesx.size(); j++) {
			for (int k = 1; k < binEdgesy.size(); k++) {
				pixelHits.push_back(h2->GetBinContent(j, k));
			}
		}
		hit_on_pixels_perEvent.push_back(pixelHits);
	}
	return hit_on_pixels_perEvent;
}

std::vector<double>
getSumOfHitsonPixel(std::vector<std::vector<double>> hit_on_pixels_perEvent) {
	std::vector<double> sums(hit_on_pixels_perEvent.at(0).size());
	for (int i = 0; i < hit_on_pixels_perEvent.size(); i++) {
		std::vector<double> pHRow = hit_on_pixels_perEvent.at(i);
		for (int j = 0; j < pHRow.size(); j++) {
			sums.at(j) = sums.at(j) + pHRow.at(j);;
		}
	}
	return sums;
}

std::vector<double> organizeMeans(std::vector<double> mean,
		std::shared_ptr<PMT> pmt) {
	int pixelN = pmt->getPixelCount();
	std::vector< double> binEdgesx;
	binEdgesx.reserve(sqrt(pixelN));
	std::vector< double> binEdgesy;
	binEdgesy.reserve(sqrt(pixelN));
	getBinEdges(pmt, &binEdgesx, &binEdgesy);


	std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>("test"));
	std::shared_ptr<TH2D> h2(std::make_shared<TH2D>(
				"Test Mean", "Test Mean", binEdgesx.size() - 1, &binEdgesx[0],
				binEdgesy.size() - 1, &binEdgesy[0]));
	int meanCount = 0;
	for (int j = 1; j < binEdgesy.size(); j++) {
		for (int k = 1; k < binEdgesx.size(); k++) {
			h2->SetBinContent(k, j, mean.at(meanCount));
			meanCount++;
		}
	}


	h2->Draw("colz");
	c1->SaveAs("test.pdf");
	std::vector<double> tempMean;
	tempMean.reserve(pixelN);
	// meanCount = 0;
	for (int j = binEdgesy.size() - 1; j > 0; j--) {
		for (int k = 1; k < binEdgesx.size(); k++) {
			// std::cout << meanCount << "   " <<  h2->GetBinContent(k,j) <<
			// std::endl;  meanCount++;
			tempMean.push_back(h2->GetBinContent(k, j));
		}
	}
	return tempMean;
}

void plot2DHistogramForHits(
		std::vector<std::vector<double>> hit_on_pixels_perEvent,
		std::shared_ptr<PMT> pmt, std::string particle, double betaMomentum,
		double theta = 0, double phi = 0, double x = 0, double y = 0,
		bool betaBool =
		false) // std::vector<std::vector<std::shared_ptr<DetectedEvent>>>
// eventsPerRunPerAerogels,std::shared_ptr<PMT> pmt)
{
	// Plots a 2D Histogram For the TOTAL Hits
	// INPUT:
	// 2D Vector of Int of Hits on pixels
	// PMT object

	// Create Path Name and File names
	std::string pathName = "./Images/Histogram2D/";

	mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	std::string fileName = pathName + "2DHist";
	fileName = fileName + particle;
	if (betaBool) {
		fileName = fileName + "beta" + std::to_string(betaMomentum);
	} else {
		fileName = fileName + "momentum" + std::to_string(betaMomentum);
	}
	fileName = fileName + "_Theta" + std::to_string(theta);
	fileName = fileName + "_Phi" + std::to_string(phi);
	fileName = fileName + "_x" + std::to_string(x);
	fileName = fileName + "_y" + std::to_string(y);

	// gStyle->SetPalette(kRainBow);

	std::string namePng = fileName + ".png";
	std::string nameRoot = fileName + ".root";
	std::string namePdf = fileName + ".pdf";

	// double dW = pmt->getWidth();
	// double dL = pmt->getLength();

	int pixelN = pmt->getPixelCount();
	std::vector< double> binEdgesx;
	binEdgesx.reserve(sqrt(pixelN));
	std::vector< double> binEdgesy;
	binEdgesy.reserve(sqrt(pixelN));
	std::vector< double> sums = getSumOfHitsonPixel(hit_on_pixels_perEvent);

	getBinEdges(pmt, &binEdgesx, &binEdgesy);
	std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>(
				fileName.c_str(), (particle + " 5 GeV").c_str(), 640, 480));
	std::shared_ptr<TH2F> h2(std::make_shared<TH2F>(
				"2D Histogram", (particle + " 5 GeV").c_str(), binEdgesx.size() - 1,
				&binEdgesx[0], binEdgesy.size() - 1, &binEdgesy[0]));

	int sumCount = 0;
	for (int j = 1; j < binEdgesy.size(); j++) {
		for (int k = 1; k < binEdgesx.size(); k++) {
			h2->SetBinContent(k, j, sums.at(sumCount));
			sumCount++;
		}
	}
	gStyle->SetOptStat(0);
	h2->GetXaxis()->SetRange(binEdgesx.front(), binEdgesx.back());
	h2->GetYaxis()->SetRange(binEdgesy.front(), binEdgesy.back());
	h2->SetStats();
	h2->Scale();
	// h2->GetXaxis()->SetNdivisions(003);
	// h2->GetYaxis()->SetNdivisions(003);
	h2->GetXaxis()->SetLabelSize(0.07);
	h2->GetYaxis()->SetLabelSize(0.07);
	c1->SetFixedAspectRatio();
	h2->GetXaxis()->SetTitle("X");
	h2->GetYaxis()->SetTitle("Y");
	c1->SetGrid();
	h2->Draw("COLZ2");
	// h2->GetXaxis()->SetRangeUser(0.04, 0.08);
	// h2->GetYaxis()->SetRangeUser(0.04, 0.08);
	gStyle->SetTitleFontSize(0.1);
	// TBox *b1 = new TBox(0.06, 0.06, 0.0660625, 0.0660625);
	// b1->SetFillColor(6);
	// b1->SetLineColor(6);
	// b1->Draw();

	// TEllipse *el1 = new TEllipse(x, y, .0005, .0005);
	// el1->SetFillColor(1);
	// el1->SetLineColor(1);
	// el1->Draw();

	c1->Modified();
	c1->SaveAs(namePng.c_str());
	c1->SaveAs(nameRoot.c_str());
	// c1->SaveAs(namePdf.c_str());
}

void plotHistogramForEachPixel(std::vector<std::vector<int>> pH) {
	// Plots a histogram for each pixel on the distribution of hits
	// INPUT
	// 2D vector of ints for the pixel Hits
	std::string pathName = "./Images/pixelHits/";

	for (size_t colNo = 0; colNo < pH.at(0).size(); ++colNo) {
		mkdir((pathName + std::to_string(colNo)).c_str(),
				S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		std::string fileName =
			pathName + std::to_string(colNo) + "/" + std::to_string(colNo) + "Hist";
		std::string namePng = fileName + ".png";
		std::string nameRoot = fileName + ".root";
		std::string namePdf = pathName + ".pdf";

		std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>(pathName.c_str()));
		std::shared_ptr<TH1D> h1(std::make_shared<TH1D>(
					pathName.c_str(),
					("Histogram of pixel" + std::to_string(colNo)).c_str(), 10, 0, 10));
		c1->SetGrid();

		for (const auto &row : pH) {
			h1->Fill(row.at(colNo));
		}
		h1->SetStats();
		h1->SetFillColor(38);
		h1->SetXTitle("Number of Detected Photons");
		h1->SetYTitle("Counts");
		h1->Draw();
		c1->Modified();
		c1->SaveAs(namePng.c_str());
		c1->SaveAs(nameRoot.c_str());
		c1->SaveAs(namePdf.c_str());
	}
}

void plotHistogramForEachPixel(std::vector<std::shared_ptr<TH1D>> hists,
		std::string particle, int momentum) {
	// Plots a histogram for each pixel on the distribution of hits
	// INPUT
	// 2D vector of TH1D histograms for the pixel Hits
	std::string pathName = "./Images/pixelHits/";
	////std::cout << hists.size() << std::endl;
	for (size_t colNo = 0; colNo < hists.size(); ++colNo) {
		mkdir((pathName + std::to_string(colNo)).c_str(),
				S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		std::string fileName = pathName + std::to_string(colNo) + "/" +
			std::to_string(colNo) + "Hist" + particle +
			std::to_string(momentum) + "GeV";
		std::string namePng = fileName + ".png";
		std::string nameRoot = fileName + ".root";
		std::string namePdf = pathName + ".pdf";

		std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>(pathName.c_str()));
		hists.at(colNo)->Draw();
		c1->Modified();
		c1->SaveAs(namePng.c_str());
		c1->SaveAs(nameRoot.c_str());
		c1->SaveAs(namePdf.c_str());
	}
}

void plotHistogramForEachPixel(std::vector<std::shared_ptr<TH1D>> hists,
		std::vector<std::shared_ptr<TGraph>> graphs,
		std::string particle, int momentum) {
	// Plots a histogram for each pixel on the distribution of hits
	// INPUT
	// 2D vector of TH1D histograms for the pixel Hits
	std::string pathName = "./Images/pixelHits/";
	mkdir(pathName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	////std::cout << hists.size() << std::endl;
	for (size_t colNo = 0; colNo < hists.size(); ++colNo) {
		mkdir((pathName + std::to_string(colNo)).c_str(),
				S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		std::string fileName = pathName + std::to_string(colNo) + "/" +
			std::to_string(colNo) + "Hist" + particle +
			std::to_string(momentum) + "GeV";
		std::string namePng = fileName + ".png";
		std::string nameRoot = fileName + ".root";
		std::string namePdf = fileName + ".pdf";

		std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>(pathName.c_str()));
		hists.at(colNo)->Draw();
		graphs.at(colNo)->Draw("SAME");
		c1->Update();
		c1->Modified();
		c1->SaveAs(namePng.c_str());
		c1->SaveAs(nameRoot.c_str());
		c1->SaveAs(namePdf.c_str());
	}
}

void plotAllPixelsOnSingleCanvas(std::vector<std::vector<int>> pH,
		std::string particle) {
	// Plots the distribution of hits on each pixel and plots it on an
	// (PixelNxPixelNy) array  Input:  A 2D vector of ints for the pixelHits

	std::string pathName = "./Images/MultiHistogramCanvas/";

	// Create Path Name and File names
	mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	std::string fileName = pathName + particle + "mulitHist";
	std::string namePng = fileName + ".png";
	std::string nameRoot = fileName + ".root";
	std::string namePdf = fileName + ".pdf";

	// Transpose the 2D vector
	std::vector<std::vector<int>> T_pH(pH[0].size(), std::vector<int>(pH.size()));
	for (int i = 0; i < pH.size(); i++) {
		for (int j = 0; j < pH.at(0).size(); j++) {
			T_pH[j][i] = pH[i][j];
		}
	}

	// Loop through each pixel (which is now the row of the T_PH 2D Vector)
	// and pass the entire row to the FillN function for the Histogram
	// The canvas must be divided into an 8by8 array
	std::vector<std::shared_ptr<TH1D>> hMulti(T_pH.size());
	std::shared_ptr<TCanvas> c1(
			std::make_shared<TCanvas>(pathName.c_str(), "", 640, 480));
	c1->DivideSquare(pH.at(0).size());

	for (int pixelCount = 0; pixelCount < T_pH.size(); pixelCount++) {
		// Have to convert the vector of ints to vector of doubles to plot
		std::vector< double> pixelHistogramValues(T_pH.at(pixelCount).begin(),
				T_pH.at(pixelCount).end());
		std::vector< double> weights(pixelHistogramValues.size(), 1.0);

		hMulti.push_back(std::make_shared<TH1D>(
					(pathName + std::to_string(pixelCount)).c_str(),
					("Histogram of pixel" + std::to_string(pixelCount)).c_str(), 10, 0,
					10));
		hMulti.back()->FillN(pixelHistogramValues.size(), &pixelHistogramValues[0],
				&weights[0]);
		hMulti.back()->SetStats();
		hMulti.back()->SetFillColor(38);
		hMulti.back()->SetXTitle("Number of Detected Photons");
		hMulti.back()->SetYTitle("Counts");
		// hMulti.back()->Scale();
		c1->cd(pixelCount + 1);
		hMulti.back()->Draw();
		gPad->Modified();
		gPad->Update();
	}
	c1->cd(0);
	c1->SaveAs(namePng.c_str());
	c1->SaveAs(nameRoot.c_str());
	c1->SaveAs(namePdf.c_str());
}

void plotAllPixelsOnSingleCanvas(std::vector<std::shared_ptr<TH1D>> hists,
		std::vector<std::shared_ptr<TGraph>> graphs,
		std::string particle, int p) {
	// Plots the distribution of hits on each pixel and plots it on an
	// (PixelNxPixelNy) array  Input:  A 2D vector of ints for the pixelHits

	std::string pathName = "./Images/MultiHistogramCanvas/";

	// Create Path Name and File names
	mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	std::string fileName =
		pathName + particle + "mulitHist" + std::to_string(p) + "GeV";
	std::string namePng = fileName + ".png";
	std::string nameRoot = fileName + ".root";
	std::string namePdf = fileName + ".pdf";

	// Loop through each pixel (which is now the row of the T_PH 2D Vector)
	// and pass the entire row to the FillN function for the Histogram
	// The canvas must be divided into an 8by8 array
	std::shared_ptr<TCanvas> c1(
			std::make_shared<TCanvas>(pathName.c_str(), "", 1200, 1200));
	c1->DivideSquare(hists.size(), 0.001, 0.001);

	for (int pixelCount = 0; pixelCount < hists.size(); pixelCount++) {
		c1->cd(pixelCount + 1);
		hists.at(pixelCount)->Draw();
		graphs.at(pixelCount)->Draw("SAME");
		gPad->Modified();
		gPad->Update();
	}
	c1->cd(0);
	c1->SaveAs(namePng.c_str());
	c1->SaveAs(nameRoot.c_str());
	c1->SaveAs(namePdf.c_str());
}

void checkPixelHits(std::vector<std::vector<int>> pH) {
	// Outputs the array of pixel hits to the terminal
	// Input: 2D vector of pixelHits
	for (int i = 0; i < pH.size(); i++) {
		std::vector<int> pHEvent = pH.at(i);
		std::cout << "E: " << i << "   ";
		for (int j = 0; j < pHEvent.size(); j++) {
			if (j % 5 == 0) {
				std::cout << j << ": " << pHEvent.at(j) << "    ";
			}
		}
		std::cout << std::endl;
	}
}

std::shared_ptr<TH1D> getHistogramOfLL(std::vector<double> ll) {
	std::vector< double> weights(ll.size(), 1.0);
	std::vector<double> ll_double(ll.begin(),ll.end());
	std::shared_ptr<TCanvas> c1(
			std::make_shared<TCanvas>(std::to_string(ll_double.at(0)).c_str()));
	std::shared_ptr<TH1D> h1(
			std::make_shared<TH1D>(std::to_string(ll_double.at(0)).c_str(),
				"Histogram of loglikelihood ratio", 10, 0, 0));
	c1->SetGrid();
	h1->FillN(ll_double.size(), &ll_double[0], &weights[0]);
	h1->SetStats();
	h1->SetFillColor(38);
	h1->SetXTitle("-2 ln(#frac{L_{#pi}}{L_{k}})");
	h1->SetYTitle("");
	h1->Draw();
	return h1;
}

double getSeperationValue(std::shared_ptr<TH1D> h1, std::shared_ptr<TH1D> h2) {
	return 1;
}

void plotBothLLHistograms(std::shared_ptr<TH1D> h1, std::shared_ptr<TH1D> h2,
		bool showSubGraphs, std::string particle1,
		std::string particle2, std::string title,
		bool betaBool, double betaMomentum, double theta,
		double phi, double x, double y) {
	double area1 = h1->Integral();
	double area2 = h2->Integral();
	h1->Scale(1 / (area1));
	h2->Scale(1 / (area2));
	h1->SetTitle(particle1.c_str());
	h1->SetTitle(particle1.c_str());
	h2->SetTitle(particle2.c_str());
	// h1->SetFillColorAlpha(38, .5);
	// h2->SetFillColorAlpha(46, .5);

	std::string pathName = "./Images/LogLikelihood/";

	// Create Path Name and File names
	mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	std::string fileName = pathName + "ll";
	if (betaBool) {
		fileName = fileName + "beta" + std::to_string(betaMomentum);
	} else {
		fileName = fileName + "momentum" + std::to_string(betaMomentum);
	}

	fileName = fileName + "_xPosition" + std::to_string(x);
	fileName = fileName + "_yPosition" + std::to_string(y);
	fileName = fileName + "_theta" + std::to_string(theta);
	fileName = fileName + "_phi" + std::to_string(phi);

	std::string namePng = fileName + ".png";
	std::string nameRoot = fileName + ".root";
	std::string namePdf = fileName + ".pdf";

	double minX1 = h1->GetXaxis()->GetXmin();
	double maxX1 = h1->GetXaxis()->GetXmax();
	// //std::cout << "H1:" << minX1 << "   " << maxX1 << std::endl;

	double minX2 = h2->GetXaxis()->GetXmin();
	double maxX2 = h2->GetXaxis()->GetXmax();
	////std::cout << "H2:" << minX2 << "   " << maxX2 << std::endl;

	// double area=0;
	// double x1 = minX1;
	// double xCount = 100;
	// double dx = (maxX1-minX1)/(xCount-1);

	// std::shared_ptr<TH1> hC1(h1->GetCumulative());
	// std::shared_ptr<TH1> hC2(h2->GetCumulative());
	std::string statTest = "NAN";

	// while(area<.900 && x1<=maxX1)
	//{
	// area = hC1->Interpolate(x1);
	//////std::cout << "x: " << x1 << "    ";
	//////std::cout << area << std::endl;
	// x1=x1+dx;
	//}
	////std::cout << "x1: " << x1<< "   ";
	////std::cout << "minX2: " << minX2 << std::endl;
	statTest = std::to_string(getSeperationValue(h1, h2));

	if (minX2 > 0) {
		std::shared_ptr<TCanvas> c1(
				std::make_shared<TCanvas>(fileName.c_str(), title.c_str(), 1200, 1200));
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		std::shared_ptr<TPaveLabel> temp1(
				std::make_shared<TPaveLabel>(0.1, 0.96, 0.9, 0.99, title.c_str()));
		temp1->Draw();
		// std::shared_ptr<TPaveLabel> temp2(std::make_shared<TPaveLabel>
		// (0.75,0.55,0.95,0.75,title.c_str()));  temp2->Draw();
		c1->SetGrid();
		c1->Divide(2, 0, 0, 0);
		c1->cd(1);
		c1->SetGrid();
		h1->GetXaxis()->SetTitle("");
		h1->GetYaxis()->SetTitle("");
		h1->Draw();
		c1->cd(2);
		c1->SetGrid();
		gPad->SetTicky(2);
		h2->GetYaxis()->SetLabelOffset(999);
		h2->GetYaxis()->SetLabelSize(0);
		h2->GetYaxis()->SetTitle("");
		h2->Draw();
		h2->GetXaxis()->SetNdivisions(-3);
		c1->cd(0);
		std::shared_ptr<TLegend> leg(
				std::make_shared<TLegend>(0.1, 0.7, 0.48, 0.9));
		leg->SetHeader("Legend"); // option "C" allows to center the header
		leg->AddEntry(h1.get(), particle1.c_str(), "f");
		leg->AddEntry(h2.get(), particle2.c_str(), "f");
		leg->AddEntry((TObject *)0,
				(std::string("Seperatation: ") + statTest).c_str(), "");
		leg->Draw();
		c1->Modified();
		c1->Update();
		// gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
		c1->SaveAs(namePng.c_str());
		c1->SaveAs(nameRoot.c_str());
		c1->SaveAs(namePdf.c_str());
	} else {
		std::shared_ptr<TCanvas> c1(
				std::make_shared<TCanvas>(fileName.c_str(), "", 1200, 1200));
		gStyle->SetOptTitle(0);
		c1->SetGrid();
		std::shared_ptr<THStack> hs(std::make_shared<THStack>(title.c_str(), ""));
		hs->Add(h1.get());
		hs->Add(h2.get());
		hs->Draw("nostack");
		hs->GetXaxis()->SetTitle("-2 ln(#frac{L_{#pi}}{L_{k}})");
		// gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

		std::shared_ptr<TLegend> leg(
				std::make_shared<TLegend>(0.1, 0.7, 0.48, 0.9));
		leg->SetHeader("Legend"); // option "C" allows to center the header
		leg->AddEntry(h1.get(), particle1.c_str(), "f");
		leg->AddEntry(h2.get(), particle2.c_str(), "f");
		leg->AddEntry((TObject *)0,
				(std::string("Seperatation: ") + statTest).c_str(), "");
		leg->Draw();
		std::shared_ptr<TPaveLabel> temp(
				std::make_shared<TPaveLabel>(0.1, 0.96, 0.9, 0.99, title.c_str()));
		temp->Draw();
		c1->Modified();
		c1->Update();
		// gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
		c1->SaveAs(namePng.c_str());
		c1->SaveAs(nameRoot.c_str());
		c1->SaveAs(namePdf.c_str());
	}

	if (showSubGraphs) {
		fileName = pathName + "ll" + particle1;
		namePng = fileName + ".png";
		nameRoot = fileName + ".root";
		namePdf = fileName + ".pdf";

		std::shared_ptr<TCanvas> c2(
				std::make_shared<TCanvas>(pathName.c_str(), "", 1200, 1200));
		c2->SetGrid();
		// h1->Scale(1/(area1));
		h1->Draw();
		c2->SaveAs(namePng.c_str());
		c2->SaveAs(nameRoot.c_str());
		c2->SaveAs(namePdf.c_str());

		fileName = pathName + "ll" + particle2;
		namePng = fileName + ".png";
		nameRoot = fileName + ".root";
		namePdf = fileName + ".pdf";

		std::shared_ptr<TCanvas> c3(
				std::make_shared<TCanvas>(pathName.c_str(), "", 1200, 1200));
		c3->SetGrid();
		h2->Draw();
		c3->SaveAs(namePng.c_str());
		c3->SaveAs(nameRoot.c_str());
		c3->SaveAs(namePdf.c_str());
	}
}

std::vector<double> getStdDevPerPixel(std::vector<std::shared_ptr<TH1D>> h) {
	std::vector<double> stdDev;
	for (int i = 0; i < h.size(); i++) {
		stdDev.push_back(h.at(i)->GetStdDev());
	}
	return stdDev;
}

void plotHistogramSep(std::vector<double> sepTest, std::string title, int p) {

	std::vector<double> sepTest_double(sepTest.begin(),sepTest.end());
	std::string pathName = "./Images/LogLikelihood/err/";

	// Create Path Name and File names
	mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	std::string fileName =
		pathName + "seperationUncertainty" + std::to_string(p) + "GeV";
	std::string namePng = fileName + ".png";
	std::string nameRoot = fileName + ".root";
	std::string namePdf = fileName + ".pdf";

	std::vector<double> weights(sepTest_double.size(), 1.0);
	std::shared_ptr<TCanvas> c1(
			std::make_shared<TCanvas>(std::to_string(sepTest_double.at(0)).c_str()));
	std::shared_ptr<TH1D> h1(
			std::make_shared<TH1D>(std::to_string(sepTest_double.at(0)).c_str(),
				"Histogram of loglikelihood ratio", 10, 0, 0));
	c1->SetGrid();
	h1->FillN(sepTest_double.size(), &sepTest_double[0], &weights[0]);
	h1->SetStats();
	h1->SetFillColor(38);
	h1->SetXTitle("SepValue");
	h1->SetYTitle("");
	h1->Draw();
	c1->Modified();
	c1->Update();
	// gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	c1->SaveAs(namePng.c_str());
	c1->SaveAs(nameRoot.c_str());
	c1->SaveAs(namePdf.c_str());
}

void SaveMeanAsTextFile(double p, double beta, double x, double y,
		std::vector<double> mean1, std::vector<double> mean2,
		std::string particle1, std::string particle2) {
	std::string pathName = "./Textfiles/";
	mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	std::string fileName = pathName + "momentum" + std::to_string(double(p));
	fileName = fileName + "_xPosition" + std::to_string(x);
	fileName = fileName + "_yPosition" + std::to_string(y);
	fileName = fileName + ".txt";
	int w = 20;
	std::ofstream ofile(fileName.c_str());
	if (ofile.is_open()) {
		ofile << "*"
			<< "Momentum:" << std::to_string(p).c_str() << "\n";
		ofile << "*"
			<< "Beta:" << std::to_string(beta).c_str() << "\n";
		ofile << "*"
			<< "xPosition:" << std::to_string(x).c_str() << "\n";
		ofile << "*"
			<< "yPosition:" << std::to_string(y).c_str() << "\n";
		ofile << "#i" << std::setw(w) << (particle1.c_str()) << std::setw(w)
			<< (particle2.c_str()) << "\n";
		for (int i = 0; i < mean1.size(); i++) {
			ofile << i << std::setw(w) << (mean1.at(i)) << std::setw(w)
				<< (mean2.at(i)) << "\n";
			// ofile << i << std::setw(w) << i << std::setw(w) << i << "\n";
		}
		ofile.close();
	}
}

void SaveMeanBetaThetaPhiAsTextFile(bool betaBool, double betaMomentum,
		double x, double y, double theta,
		double phi, std::vector<double> mean1,
		std::vector<double> mean2,
		std::string particle1,
		std::string particle2) {
	std::string pathName = "./Textfiles/";
	mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	std::string fileName = pathName;
	if (betaBool) {
		fileName = fileName + "beta" + std::to_string(betaMomentum);
	} else {
		fileName = fileName + "momentum" + std::to_string(betaMomentum);
	}
	fileName = fileName + "_xPosition" + std::to_string(x);
	fileName = fileName + "_yPosition" + std::to_string(y);
	fileName = fileName + "_theta" + std::to_string(theta);
	fileName = fileName + "_phi" + std::to_string(phi);
	fileName = fileName + ".txt";
	int w = 20;
	std::ofstream ofile(fileName.c_str(), std::ios::out | std::ios::binary);
	if (ofile.is_open()) {
		if (betaBool) {
			ofile << "*"
				<< "Beta:" << std::to_string(betaMomentum).c_str() << "\n";
		} else {
			ofile << "*"
				<< "Momentum:" << std::to_string(betaMomentum).c_str() << "\n";
		}
		ofile << "*"
			<< "xPosition:" << std::to_string(x).c_str() << "\n";
		ofile << "*"
			<< "yPosition:" << std::to_string(y).c_str() << "\n";
		ofile << "*"
			<< "Theta:" << std::to_string(theta).c_str() << "\n";
		ofile << "*"
			<< "Phi:" << std::to_string(phi).c_str() << "\n";
		ofile << "#i" << std::setw(w) << (particle1.c_str());
		if (!particle2.empty()) {
			ofile << std::setw(w) << (particle2.c_str());
		}
		ofile << "\n";
		for (int i = 0; i < mean1.size(); i++) {
			ofile << i << std::setw(w) << (mean1.at(i));
			if (!mean2.empty()) {
				ofile << std::setw(w) << (mean2.at(i));
			}
			ofile << "\n";
		}
		ofile.close();
	}
}

std::string getFileName(bool betaBool, std::shared_ptr<Beam> beam) {
	std::string fileName;
	if (betaBool) {
		fileName = fileName + "beta" + std::to_string(beam->getBeta());
	} else {
		fileName = fileName + "momentum" + std::to_string(beam->getMomentum());
	}

	fileName =
		fileName + "_xPosition" + std::to_string(beam->getPosition()->getX());
	fileName =
		fileName + "_yPosition" + std::to_string(beam->getPosition()->getY());
	fileName =
		fileName + "_theta" + std::to_string(beam->getDirection()->getTheta());
	fileName = fileName + "_phi" + std::to_string(beam->getDirection()->getPhi());
	return fileName;
}


std::string getFileNameBeta(bool betaBool, std::shared_ptr<Beam> beam) {
	std::string fileName;

	fileName =
		fileName + "xPosition" + std::to_string(beam->getPosition()->getX());
	fileName =
		fileName + "_yPosition" + std::to_string(beam->getPosition()->getY());
	fileName =
		fileName + "_theta" + std::to_string(beam->getDirection()->getTheta());
	fileName = fileName + "_phi" + std::to_string(beam->getDirection()->getPhi());
	return fileName;
}


void SaveMultipleMeanQuick(bool betaBool, std::vector<std::shared_ptr<Beam>> beamList,
		std::vector<std::vector<double>> mean, std::string path, std::string file = "") {
	int VARIABLE_SIZE = 5;

	std::cout <<"CALLED TO PRINT" << std::endl;

	std::string pathName;
	if (path.empty()) {
		pathName = "./Textfiles/";
	} else {
		pathName = path;
	}
	std::cout << pathName << std::endl;
	mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	std::string fileName = "";//pathName;
	if (file==""){
		fileName += pathName;
		fileName += getFileNameBeta(betaBool, beamList.back());
		fileName += ".txt";
		std::cout << fileName << std::endl;
	}
	else{
		fileName +=  file;
		fileName += ".txt";
	}
	FILE *f = fopen(fileName.c_str(), "wt");
	// assert(f); //does f exist

	double v[VARIABLE_SIZE];
	// std::cout << std::endl;

	std::cout << mean.size() << std::endl;
	for(int loop = 0; loop < mean.size(); loop++)
	{
		std::cout << loop << "   " << beamList[loop]->getBeta() << "  " <<  mean[loop].size()<< std::endl;
		std::shared_ptr<Beam> beam = beamList[loop];
		if (betaBool) {
			v[0] = beam->getBeta();
		} else {
			v[0] = beam->getMomentum();
		}
		v[1] = beam->getPosition()->getX();
		v[2] = beam->getPosition()->getY();
		v[3] = beam->getDirection()->getTheta();
		v[4] = beam->getDirection()->getPhi();
		for (int variable = 0; variable < VARIABLE_SIZE+ mean[loop].size(); variable++) {
			// std::cout << v[variable] << " ";
			if (variable < VARIABLE_SIZE){
				fprintf(f, "%g ", double( v[variable]));
			}
			else{
				fprintf(f, "%g ", double(mean[loop][variable - VARIABLE_SIZE]));
			}
			// std::cout << mean[pixel] << std::endl;
		}
		fprintf(f,"\n");
	}
	fprintf(f, "\n");
	// fwrite("\n", sizeof(char), 1, f);
	//}
	fclose(f);
	}



void SaveMeanQuickNoLines(bool betaBool, std::shared_ptr<Beam> beam,
		std::vector<double> mean, std::string path) {
	int VARIABLE_SIZE = 5;

	std::string pathName;
	if (path.empty()) {
		pathName = "./Textfiles/";
	} else {
		pathName = path;
	}
	mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	std::string fileName = pathName;
	fileName += getFileName(betaBool, beam);
	fileName += ".txt";
	FILE *f = fopen(fileName.c_str(), "wt");
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
	for (int variable = 0; variable < VARIABLE_SIZE + mean.size(); variable++) {
		// std::cout << v[variable] << " ";
		if (variable < VARIABLE_SIZE){
			fprintf(f, "%g ", double(v[variable]));
		}
		else{
			fprintf(f, "%g ",double( mean[variable-VARIABLE_SIZE]));
		}
	}
	fprintf(f, "\n");
	fclose(f);
}

void SaveMeanQuick(bool betaBool, std::shared_ptr<Beam> beam,
		std::vector<double> mean, std::string path) {
	int VARIABLE_SIZE = 5;

	std::string pathName;
	if (path.empty()) {
		pathName = "./Textfiles/";
	} else {
		pathName = path;
	}
	mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	std::string fileName = pathName;
	fileName += getFileName(betaBool, beam);
	fileName += ".txt";
	FILE *f = fopen(fileName.c_str(), "wt");
	// assert(f); //does f exist

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
	for (int variable = 0; variable < VARIABLE_SIZE; variable++) {
		// std::cout << v[variable] << " ";
		fprintf(f, "%g ", double(v[variable]));
	}
	// std::cout << std::endl;
	fprintf(f, "\n");
	for (int pixel = 0; pixel < mean.size(); pixel++) {
		fprintf(f, "%g\n",double( mean[pixel]));
		// std::cout << mean[pixel] << std::endl;
	}
	fprintf(f, "\n");
	// fwrite("\n", sizeof(char), 1, f);
	//}
	fclose(f);
	}

double calcStepSize(double start, double end, double step) {
	return (end - start) / (step - 1);
}
std::vector<std::vector<double>> getJobLimits(std::string fileName) {
	// std::regex pattern1("Jobs/job");
	// std::regex pattern2(".txt");
	// jobName = std::regex_replace(jobName,pattern1,"");
	// jobName = std::regex_replace(jobName,pattern2,"");
	// std::cout << jobName << std::endl;
	// int jobnum = std::stoi(jobName);
	std::vector<std::vector<double>> limitVariables;
	std::ifstream infile(fileName.c_str());
	while (infile.good()) {
		std::vector<double> temp;
		std::string line;
		getline(infile, line);
		std::istringstream iss(line);
		std::string variableName;
		int job, stupid;
		double v1, v2, v3, v4, v5;
		iss >> job >> variableName >> stupid >> v1 >> v2 >> v3 >> v4 >> v5;
		// if (job == jobnum){
		// std::cout << v1 << " " << v2 << " " << v3 << " " << v4 << " " << v5 <<
		// std::endl;
		temp.push_back(v1);
		temp.push_back(v2);
		temp.push_back(v3);
		temp.push_back(v4);
		temp.push_back(v5);
		limitVariables.push_back(temp);
		// return limitVariables;
		//}
	}

	return limitVariables;
}


std::vector<double> linspace(double begin, double end, double steps)
{

	std::vector<double> v;
	v.reserve(steps);

	double delta = (end-begin)/(steps-1);

	double value = begin;
	for (int i = 0; i < steps; i++)
	{
		v.push_back(value);
		value+=delta;
	}
	return v;
}

void setVariablesFromJobFilesBeta(std::string fileName,
		std::vector<double> &xVariables,
		std::vector<double> &yVariables,
		std::vector<double> &thetaVariables,
		std::vector<double> &phiVariables,
		std::vector<double> &betaMomentumVariables) {
	// std::cout << fileName << std::endl;
	std::ifstream infile(fileName.c_str());
	if (infile.is_open()) {
		std::cout << "TEST" << std::endl;
		while (infile.good()) {
			std::string line;
			getline(infile, line);
			std::cout << line << std::endl;
			std::istringstream iss(line);
			double minX, minY, minTheta,minPhi=0;
			iss >> minX >> minY >> minTheta >> minPhi;

			// std::cout << "X" << std::endl;
			xVariables.at(0) = minX;
			xVariables.at(1) = minX+0.1;
			xVariables.at(2) = 1;//calcStepSize(min, max, step);
			for (int i = 0; i < xVariables.size();i++)
			{
				std::cout << xVariables.at(i) << "  ";
			}					
			std::cout << std::endl;
			yVariables.at(0) = minY;
			yVariables.at(1) = minY+0.1;
			yVariables.at(2) = 1.0; //calcStepSize(min, max, step);
			for (int i = 0; i < yVariables.size();i++)
			{
				std::cout << yVariables.at(i) << "  ";
			}					
			std::cout << std::endl;

			thetaVariables.at(0) = minTheta;
			thetaVariables.at(1) = minTheta+0.1;
			thetaVariables.at(2) = 1;//calcStepSize(min, max, step);
			for (int i = 0; i < thetaVariables.size();i++)
			{
				std::cout << thetaVariables.at(i) << "  ";
			}					
			std::cout << std::endl;
			phiVariables.at(0) = minPhi;
			phiVariables.at(1) = minPhi+0.1;
			phiVariables.at(2) = 1;//calcStepSize(min, max, step);
			for (int i = 0; i < phiVariables.size();i++)
			{
				std::cout << phiVariables.at(i) << "  ";
			}					
			std::cout << std::endl;
			std::vector<double> betaValues1 = linspace(0.988,0.999,20);
			std::vector<double> betaValues2 = linspace(0.999,1.0,15);
			betaValues2.erase(betaValues2.begin()+0);
			betaValues1.insert(betaValues1.end(),betaValues2.begin(),betaValues2.end());
			betaMomentumVariables = betaValues1;
			for (int i = 0; i < betaMomentumVariables.size();i++)
			{
				std::cout << betaMomentumVariables.at(i) << "  ";
			}					
			std::cout << std::endl;

		}
	}

}

void setVariablesFromJobFiles(std::string fileName,
		std::vector<double> *xVariables,
		std::vector<double> *yVariables,
		std::vector<double> *thetaVariables,
		std::vector<double> *phiVariables,
		std::vector<double> *betaMomentumVariables) {

	std::ifstream infile(fileName.c_str());
	if (infile.is_open()) {
		std::cout << "TEST" << std::endl;
		while (infile.good()) {
			std::string line;
			getline(infile, line);
			// std::cout << line << std::endl;
			if (line.find("-") == std::string::npos) {
				std::istringstream iss(line);
				std::string variableName;
				double min, max, step;
				iss >> variableName >> min >> max >> step;

				if ("X" == variableName) {
					// std::cout << "X" << std::endl;
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
			}
		}
	}
}

std::vector<double> quickMean(std::vector<std::vector<double>> pixelHits1) {
	std::vector<double> mean(pixelHits1.at(0).size(), 0.0);
	for (int i = 0; i < pixelHits1.size(); i++) {
		std::vector<double> eventsPerPxiel = pixelHits1.at(i);
		double eventsMeans = 0;
		// std::cout << eventsPerPxiel.size() << std::endl;
		for (int j = 0; j < eventsPerPxiel.size(); j++) {
			mean.at(j) = eventsPerPxiel.at(j) / pixelHits1.size() + mean.at(j);
		}
	}
	return mean;
}

std::vector<std::vector<double>>
quickOrg(std::vector<std::vector<double>> pixelHits) {
	std::vector<std::vector<double>> tempPH;
	for (int i = 0; i < pixelHits.size(); i++) {
		std::vector<double> tempPHevent(pixelHits.at(i).size(), 0.0);
		for (int j = 0; j < sqrt(pixelHits.at(i).size()); j++) {
			for (int k = 0; k < sqrt(pixelHits.at(i).size()); k++) {
				std::cout << (sqrt(pixelHits.at(i).size()) - j) *
					(sqrt(pixelHits.at(i).size())) -
					(sqrt(pixelHits.at(i).size()) - k)
					<< std::endl;
				tempPHevent.at((sqrt(pixelHits.at(i).size()) - j) *
						(sqrt(pixelHits.at(i).size())) -
						(sqrt(pixelHits.at(i).size()) - k)) =
					pixelHits[i][j * sqrt(pixelHits.at(i).size()) + k];
			}
		}
		tempPH.push_back(tempPHevent);
	}
	return tempPH;
}

bool exists_test1(std::string name) {
	if (FILE *file = fopen(name.c_str(), "r")) {
		fclose(file);
		return true;
	} else {
		return false;
	}
}

bool isLimitVariables(double beta, double x, double y, double t, double p,
		std::vector<std::vector<double>> variableLimits) {
	for (int i = 0; i < variableLimits.size(); i++) {
		if (variableLimits[i][0] == beta) {
			bool limited = beta <= variableLimits[i][0];    // beta
			limited = limited && x <= variableLimits[i][1]; // x
			limited = limited && y <= variableLimits[i][2]; // y
			limited = limited && t <= variableLimits[i][3]; // t
			limited = limited && p <= variableLimits[i][4]; // p
			std::cout << limited << std::endl;
			return limited;
		}
	}
}

std::vector<double>
getLimitVariables(double beta, double x, double y, double t, double p,
		std::vector<std::vector<double>> variableLimits) {
	double a[5] = {beta, x, y, t, p};
	int savedIndex = 0;
	double savedCounter =10000;
	for (int i = 0; i < 5; i++){std::cout << a[i] << "***";}
	std::cout << std::endl;
	std::vector<std::vector<double> > saveBeta;
	for (int i = 0; i < variableLimits.size(); i++) {
		std::cout << variableLimits[i][0] << "    " << beta << std::endl;	
		if (std::fabs(variableLimits[i][0]- beta) < 0.00001 ){saveBeta.push_back(variableLimits[i]);}
	}
	for (int i = 0 ; i < saveBeta.size(); i++)
	{
		double counter =  std::fabs(saveBeta[i][3]-t);
		if (counter < savedCounter && saveBeta[i][3]<= t){
			savedIndex = i;
			savedCounter = counter;
		}
	}
	for (int i = 0; i < 5; i++){std::cout << saveBeta[savedIndex][i] << "***";}
	std::cout << std::endl;
	return saveBeta[savedIndex];
}

int main(int argc, char *argv[]) {


	std::cout << "Arg 1: Job File" << std::endl;
	std::cout << "Default: Beta:1000, Theta:0, Phi:0, beamX:0.097m, beamY:0.097m, beamZ:0.00m" << std::endl;
	std::cout << "Arg 2: Save Path" << std::endl;
	std::cout << "Default: $PWD" << std::endl;

	double aeroWidth = .12;     // m
	double aeroLength = .12;    // m
	double aeroThickness = .02; // m

	// Refractive Index
	double n1 = 1.010;
	double n2 = 1.012;

	/*Aerogel Parameters*/
	std::shared_ptr<Aerogel> aero1(
			std::make_shared<Aerogel>(aeroWidth, aeroLength, aeroThickness, n1));
	std::shared_ptr<Aerogel> aero2(
			std::make_shared<Aerogel>(aeroWidth, aeroLength, aeroThickness, n2));
	std::vector<std::shared_ptr<Aerogel>> aerogels;
	aerogels.push_back(aero1);
	aerogels.push_back(aero2);

	/*PMT Parameters*/
	double multiplePMTs = 4; // THIS IS A HACK IDEALLY IT SHOULD BE THAT YOU MAKE
				 // MULTIPLE PMT Objects
	double pmtWidth = multiplePMTs * 48.5 / 1000;  // m
	double pmtLength = multiplePMTs * 48.5 / 1000; // m
	int sidePixelCount = 8;
	int pixelN = std::pow(multiplePMTs * sidePixelCount, 2);
	std::shared_ptr<PMT> pmt(std::make_shared<PMT>(pmtWidth, pmtLength, pixelN));

	std::string path;
	int eventCount = 100000;

	// Momentum
	double startMomentum = 1000;//MeV
	double endMomentum = 10000;//MeV
	double momentumCount = 1;
	double stepMomentum = calcStepSize(startMomentum, endMomentum, momentumCount);

	// Beta
	double startBeta = 0.99999; // FORCES TO BE SHOWN
	double endBeta = 1.0;
	double betaCount = 1;
	double stepBeta = calcStepSize(startBeta, endBeta, betaCount);
	bool betaBool = true;

	// Phi
	double startPhi = 0;
	double endPhi = TMath::Pi() / 4; // 2*TMath::Pi()/4;//radians => max = 2*PI
	double phiCount = 1;
	double stepPhi = calcStepSize(startPhi, endPhi, phiCount);
	
	// Theta
	double startTheta = 0; //-TMath::Pi()/4;//radians
	double endTheta = .4;    // TMath::Pi()/4;//radians
	double thetaCount = 1;
	double stepTheta = calcStepSize(startTheta, endTheta, thetaCount);
	
	// Mass
	double pionMass = 137.2735; // MeV/c^2
	double kaonMass = 495.665;  // MeV/c^2
	// Position
	double beamx = 0.097;               // center m //#CHANGE=.06
	double beamy = 0.097;               // center m //#CHANGE=.06
	double pixelDist = 48.5 / (8000.0); // dimensions of a pixel
	
	int beamCountX = 1;
	int beamCountY = 1;
	double beamStepX = calcStepSize(beamx, beamx + pixelDist, beamCountX);
	double beamStepY = calcStepSize(beamy, beamy + pixelDist, beamCountY);
	double beamz = 0.00; // m
 
	/*Detector Parameters*/
	std::shared_ptr<Detector> detector(
			std::make_shared<Detector>(aerogels, pmt, .09));
 
	double startI = startMomentum;
	double endI = endMomentum;
	double stepI = stepMomentum;

	if (betaBool) {
		startI = startBeta;
		endI = endBeta;
		stepI = stepBeta;
	}
	gROOT->SetBatch(kTRUE);

	std::string filename;

	std::vector<double> xVariables;
	std::vector<double> yVariables;
	std::vector<double> thetaVariables;
	std::vector<double> phiVariables;
	std::vector<double> betaMomentumVariables;

	xVariables.push_back(beamx);
	xVariables.push_back(beamx + pixelDist);
	xVariables.push_back(beamStepX);

	yVariables.push_back(beamy);
	yVariables.push_back(beamy + pixelDist);
	yVariables.push_back(beamStepY);

	thetaVariables.push_back(startTheta);
	thetaVariables.push_back(endTheta);
	thetaVariables.push_back(stepTheta);

	phiVariables.push_back(startPhi);
	phiVariables.push_back(endPhi);
	phiVariables.push_back(stepPhi);

	if (betaBool) {
		betaMomentumVariables.push_back(startBeta);
		betaMomentumVariables.push_back(endBeta);
		betaMomentumVariables.push_back(stepBeta);
	} else {
		betaMomentumVariables.push_back(startMomentum);
		betaMomentumVariables.push_back(endMomentum);
		betaMomentumVariables.push_back(stepMomentum);
	}

	if (argc > 1) {
		std::string filename(argv[1]);
		std::cout << filename << std::endl;
		//This file assumes beta ranges from  0.988-0.999 (20) 0.999-1.0 (15)
		setVariablesFromJobFilesBeta(filename, xVariables, yVariables,
				thetaVariables, phiVariables,
				betaMomentumVariables);			
	}
	
	std::string pathAllfiles;
	if (argc >2) {
		pathAllfiles = std::string(argv[2]);
		std::cout << pathAllfiles<< std::endl;
	}

	int loops = 1;
	std::vector<std::shared_ptr<Beam>> loopBeam;
	loopBeam.reserve(betaMomentumVariables.size());

	std::vector<std::vector<double>> loopMeans;//(betaMomentumVariables.size(),std::vector<double>(std::pow(multiplePMTs * sidePixelCount+2, 2)));
	

	for(int lo = 0; lo < loops; lo++){
		for (double x = xVariables.at(0); x <= xVariables.at(1);
				x = x + xVariables.at(2)) {
			for (double y = yVariables.at(0); y <= yVariables.at(1);
					y = y + yVariables.at(2)) {
				for (double t = thetaVariables.at(0); t <= 0.01 + thetaVariables.at(1);
						t = t + thetaVariables.at(2)) {
					for (double p = phiVariables.at(0); p <= 0.01 + phiVariables.at(1);
							p = p + phiVariables.at(2)) {
							std::shared_ptr<Position> pos(
									std::make_shared<Position>(x, y, beamz));
							std::shared_ptr<Direction> dir(std::make_shared<Direction>(t, p));
						for (int i = 0;
								i < betaMomentumVariables.size();
								i++){ // I index can either be momentum or beta

							/*
							   std::cout << i << "   ";
							   std::cout << x << "   ";
							   std::cout << y << "   ";
							   std::cout << t << "   ";
							   std::cout << p << std::endl;;
							   */
							/*		if (x < trueLimits[1]) {
									continue;
									} else if (y < trueLimits[2]) {
									continue;
									} else if (i < trueLimits[0]) {
									continue;
									} else if (t < trueLimits[3]) {
									continue;
									} else if (p < trueLimits[4]) {
									continue;
									}
									*/

							std::shared_ptr<Beam> beamPion;
							//std::cout << betaMomentumVariables[i] << std::endl;
							if (betaBool) {
								beamPion = std::make_shared<Beam>(betaMomentumVariables[i], dir, pos);
								beamPion->setParticleMass(pionMass);
							} else {
								beamPion = std::make_shared<Beam>(betaMomentumVariables[i], pionMass, dir, pos);
							}
							std::string fileName = pathAllfiles;
							fileName += getFileNameBeta(betaBool, beamPion);
							fileName += ".txt";
							std::cout << "x: " << x << "    ";
							std::cout << "y: " << y << "    ";
							std::cout << "pionP: " << beamPion->getMomentum() << "     ";
							std::cout << "pionB: " << beamPion->getBeta() << "     ";
							std::cout << "theta: " << t << "    ";
							std::cout << "phi: " << p << "    ";
							std::cout << "\n";
							std::cout << fileName << std::endl;
							std::clock_t begin_single = std::clock();
							std::vector<std::vector<std::shared_ptr<DetectedEvent>>> dEPion =
								detector->getDetectedEvents(beamPion, eventCount);
							std::clock_t end_single = std::clock();
								double elapsed_secs_single =  double(end_single- begin_single) / CLOCKS_PER_SEC;
							std::cout << "TIME THE CODE TOOK SINGLE:" << elapsed_secs_single << std::endl;
				
							std::vector<std::vector<double>> pixelHitsPion =
							detector->getPixelHitsPerEvents(dEPion, false);
							
							std::vector<double> meanPion =
								detector->getMeanPerPixel(pixelHitsPion);


							//if(loops==1) {SaveMeanQuick(betaBool, beamPion, meanPion, path+"unorganized");}
							//meanPion = organizeMeans(meanPion, pmt);
							/*	
								for (int l = 0; l < meanPion.size(); l++)
								{
								std::cout << l << "   " << meanPion[l] << std::endl;
								}*/
							loopBeam.push_back(beamPion);
							loopMeans.push_back(meanPion);
							//std::clock_t startLine = std::clock();
							//if(loops==1) {SaveMeanQuick(betaBool, beamPion, meanPion, path+"lines");}
							//std::clock_t endLine = std::clock();
							//std::cout << "Time of Line: " << (endLine - startLine) / CLOCKS_PER_SEC;
							
							
							//No Lines
							//std::clock_t startNoLine = std::clock();
							//if(loops==1) {SaveMeanQuickNoLines(betaBool,beamPion, meanPion, path+"no_lines");}
							//std::clock_t endNoLine = std::clock();
							//std::cout << "Time of Line: " << (endLine - startLine) / CLOCKS_PER_SEC;
						}
					}
				}
			}
		}
	}

	SaveMultipleMeanQuick(betaBool,loopBeam,loopMeans,pathAllfiles);
}

