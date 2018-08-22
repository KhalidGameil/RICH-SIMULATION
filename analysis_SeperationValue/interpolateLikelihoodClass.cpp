#include "interpolateLikelihoodClass.h"

InterpolateLikelihood::InterpolateLikelihood() {
}

//InterpolateLikelihood::InterpolateLikelihood(std::vector<double> variables) {
//readVariableLists();
//testIndexFinder(variables);
//std::cout << "Completed Reading datafile" << std::endl;
//}

InterpolateLikelihood::InterpolateLikelihood(std::string datafile) {
  std::cout << datafile.c_str() << std::endl;
  readDataBaseFileFast(datafile);
  std::cout << "Completed Reading datafile" << std::endl;

  //sortVariables();
}

//InterpolateLikelihood::InterpolateLikelihood(std::string dataFile,
//int xCombinations, int yCombinations,
//int thetaCombinations, int phiCombinations,
//int betaCombinations) {
//xList.reserve(xCombinations);
//yList.reserve(yCombinations);
//thetaList.reserve(thetaCombinations);
//phiList.reserve(phiCombinations);
//betaList.reserve(betaCombinations);
//measurements.reserve(xCombinations * yCombinations * thetaCombinations * phiCombinations * betaCombinations);
//readDataBaseFile(dataFile);
//std::cout << "Completed Reading datafile" << std::endl;
//}

void InterpolateLikelihood::addMeasurement(double x, double y, double theta, double phi, double beta, std::vector<double> means) {
  measuredBatch newBatch;
  newBatch.x = x;
  newBatch.y = y;
  newBatch.theta = theta;
  newBatch.phi = phi;
  newBatch.beta = beta;
  newBatch.mean.push_back(means);
  measurements.push_back(newBatch);
}

void InterpolateLikelihood::readDataBaseFileFast(std::string filename) {
  std::ifstream infile(filename.c_str());
  bool newConfig = false;

  if (infile.is_open()) {
    measuredBatch currentBatch;

    std::vector<double> means;
    int lineCounter = 1;
    int beamChunk = 0;
    int chunkCount = 0;
    int valueP = 0;
    int pixelCount = 0;
    std::string line;
    while (infile.good()) {
      std::string line;
      std::getline(infile, line);
      double value0, value1, value2, value3, value4 = 0;
      std::istringstream iss(line);

      if (lineCounter == 1) {
        std::cout << lineCounter << std::endl;
        iss >> value0 >> value1 >> value2 >> value3 >> value4 >> valueP;
        std::cout << value0 << "    " << value1 << "    " << value2 << "     " << value3 << "     " << value4 << "      " << valueP << std::endl;
        betaList.reserve(value0 * value1 * value2 * value3 * value4);
        xList.reserve(value0 * value1 * value2 * value3 * value4);
        yList.reserve(value0 * value1 * value2 * value3 * value4);
        thetaList.reserve(value0 * value1 * value2 * value3 * value4);
        phiList.reserve(value0 * value1 * value2 * value3 * value4);
        measurements.reserve(value0 * value1 * value2 * value3 * value4);
        means.reserve(valueP);
        beamChunk = valueP + 2; //one additional count for the space and another for the configuration
      } else if (lineCounter == 2 + beamChunk * chunkCount) {
        //std::cout << lineCounter << std::endl;
        iss >> value0 >> value1 >> value2 >> value3 >> value4;
        betaList.push_back(value0);
        currentBatch.beta = value0;

        xList.push_back(value1);
        currentBatch.x = value1;

        yList.push_back(value2);
        currentBatch.y = value2;

        thetaList.push_back(value3);
        currentBatch.theta = value3;

        phiList.push_back(value4);
        currentBatch.phi = value4;

      } else if (lineCounter == beamChunk + 1 + beamChunk * chunkCount) {
        pixelCount = 0;
        //std::cout << lineCounter << "   ";
        //std::cout << means.size() << std::endl;
        currentBatch.mean.push_back(means);
        measurements.push_back(currentBatch);
        currentBatch = measuredBatch();
        means.clear();
        means.reserve(valueP);
        chunkCount++;
      } else {
        //std::cout << lineCounter - beamChunk * valueP << "   ";
        iss >> value0;
        //std::cout << value0 << std::endl;
        means.push_back(value0);
        pixelCount++;
      }
      lineCounter++;
    }
  } else {
    std::cout << "Problem reading file" << std::endl;
  }
  return;
}

void InterpolateLikelihood::readDataBaseFileFastC(std::string filename) {
  //std::clrscr();
  FILE *fin;

  fin = fopen(filename.c_str(), "r");

  if (fin == NULL) {
    std::cout << "NO FILE" << std::endl;
    return;
  }

  measuredBatch currentBatch;

  std::vector<double> means;
  int lineCounter = 1;
  int configChunk = 0;
  int chunkCount = 0;
  int valueP = 0;
  int pixelCount = 0;
  while (!feof(fin)) {
    double value0, value1, value2, value3, value4 = 0;

    if (lineCounter == 1) {
      fscanf(fin, "%lf %lf %lf %lf %lf %d",
             &value0, &value1, &value2, &value3, &value4, &valueP);
      betaList.reserve(value0);
      xList.reserve(value1);
      yList.reserve(value2);
      phiList.reserve(value3);
      thetaList.reserve(value4);
      measurements.reserve(value0 * value1 * value2 * value3 * value4);
      means.reserve(valueP);
      configChunk = valueP + 2; //one additional count for the space and another for the configuration
    } else if (lineCounter == 2 + configChunk * chunkCount) {
      fscanf(fin, "%lf %lf %lf %lf %lf", &value0, &value1, &value2, &value3, &value4);

      betaList.push_back(value0);
      currentBatch.beta = value0;

      xList.push_back(value1);
      currentBatch.x = value1;

      yList.push_back(value2);
      currentBatch.y = value2;

      thetaList.push_back(value3);
      currentBatch.theta = value3;

      phiList.push_back(value4);
      currentBatch.phi = value4;

    } else if (lineCounter == configChunk + 1 + configChunk * chunkCount) {
      currentBatch.mean.push_back(means);
      measurements.push_back(currentBatch);
      currentBatch = measuredBatch();
      means.clear();
      means.reserve(valueP);

      chunkCount++;
    } else {
      fscanf(fin, "%lf", &value0);
      means[pixelCount] = value0;
      pixelCount++;
    }
    lineCounter++;
  }
  return;
}

double InterpolateLikelihood::getMinX() { return *std::min_element(xList.begin(), xList.end()); }
double InterpolateLikelihood::getMaxX() { return *std::max_element(xList.begin(), xList.end()); }

double InterpolateLikelihood::getMinY() { return *std::min_element(yList.begin(), yList.end()); }
double InterpolateLikelihood::getMaxY() { return *std::max_element(yList.begin(), yList.end()); }

double InterpolateLikelihood::getMinTheta() { return *std::min_element(thetaList.begin(), thetaList.end()); }
double InterpolateLikelihood::getMaxTheta() { return *std::max_element(thetaList.begin(), thetaList.end()); }

double InterpolateLikelihood::getMinPhi() { return *std::min_element(phiList.begin(), phiList.end()); }
double InterpolateLikelihood::getMaxPhi() { return *std::max_element(phiList.begin(), phiList.end()); }

double InterpolateLikelihood::getMinBeta() { return *std::min_element(betaList.begin(), betaList.end()); }
double InterpolateLikelihood::getMaxBeta() { return *std::max_element(betaList.begin(), betaList.end()); }

void InterpolateLikelihood::readDataBaseFile(std::string filename) {
  std::ifstream infile(filename.c_str());
  bool newConfig = false;

  if (infile.is_open()) {
    measuredBatch currentBatch;

    std::vector<double> means;
    means.reserve(256);
    while (infile.good()) {
      double m1 = 0;
      double m2 = 0;
      int i = 0;
      std::string line;
      getline(infile, line);

      if (line.find("#") == std::string::npos &&
          line.find("*") == std::string::npos && !line.empty()) {
        std::istringstream iss(line);
        iss >> i >> m1 >> m2;
        //std::cout << i << std::endl;
        means.push_back((m1 + m2) / 2);
        if (i == 255) {
          newConfig = true;
        }

      } else if (line.find("*") != std::string::npos) {
        std::string delimeter = ":";
        int pos = line.find(delimeter);
        std::string identifier = line.substr(0, pos);
        line.erase(0, pos + delimeter.size());
        if (identifier.find("Momentum") != std::string::npos) {
          currentBatch.p = stod(line);
        } else if (identifier.find("xPosition") != std::string::npos) {
          double x = stod(line);
          if (x == 0.066062) {
            x = 0.0660625;
          }
          currentBatch.x = x;
          xList.push_back(x);
        } else if (identifier.find("yPosition") != std::string::npos) {
          double y = stod(line);
          if (y == 0.066062) {
            y = 0.0660625;
          }
          currentBatch.y = y;
          yList.push_back(y);
        } else if (identifier.find("Theta") != std::string::npos) {
          currentBatch.theta = stod(line);
          thetaList.push_back(stod(line));
        } else if (identifier.find("Phi") != std::string::npos) {
          currentBatch.phi = stod(line);
          phiList.push_back(stod(line));
        } else if (identifier.find("Beta") != std::string::npos) {
          currentBatch.beta = stod(line);
          betaList.push_back(stod(line));
          currentBatch.p = 0;
        }
      }
      if (newConfig) {
        currentBatch.mean.push_back(means);
        measurements.push_back(currentBatch);
        currentBatch = measuredBatch();
        means.clear();
        means.reserve(256);
        newConfig = false;
      }
    }
  }
}

std::string InterpolateLikelihood::getFolderLocation() {
  return folderLocation;
}

void InterpolateLikelihood::setFolderLocation(std::string fL) {
  this->folderLocation = fL;
}

void InterpolateLikelihood::readTextFile(std::string fileName) {
  measuredBatch currentBatch;
  std::ifstream infile(fileName.c_str());
  if (infile.is_open()) {
    double i;
    std::vector<double> means1;
    std::vector<double> means2;
    std::vector<double> means3;

    double m1;
    double m2;
    double m3;

    int particleCounts = 0;
    while (infile.good()) {
      std::string line;
      getline(infile, line);
      ////std::cout << line.c_str() << std::endl;
      ////std::cout << line.c_str() << std::endl;
      //std::cout << "TEST" << std::endl;
      //
      if (line.find("#") != std::string::npos) {
        char *pch;
        std::vector<char> writable(line.begin(), line.end());
        writable.push_back('\0');
        pch = strtok(&writable[0], " ");
        while (pch != NULL) {
          printf("%s\n", pch);
          if (strstr(pch, "#") == NULL) {
            currentBatch.particle.push_back(pch);
            if (strstr(pch, "kaon") != NULL) {
              currentBatch.mass.push_back(kaonMass);
            } else if (strstr(pch, "pion") != NULL) {
              currentBatch.mass.push_back(pionMass);
            } else if (strstr(pch, "proton") != NULL) {
              currentBatch.mass.push_back(protonMass);
            }

            particleCounts++;
          }
          pch = strtok(NULL, " ");
        }
        ////std::cout << particleCounts << std::endl;
      }
      //std::cout << "TEST" << std::endl;
      if (line.find("#") == std::string::npos &&
          line.find("*") == std::string::npos && !line.empty()) {
        std::istringstream iss(line);
        ////std::cout << particleCounts<< "   " <<  line.c_str() << std::endl;
        if (particleCounts == 2) //if there are two particles in the file
        {
          iss >> i >> m1 >> m2;
          means1.push_back(m1);
          means2.push_back(m2);
          //std::cout << i << "   ";
          //std::cout << m1 << "   ";
          //std::cout << m2 << std::endl;
        } else if (particleCounts == 3) // if there are three
        {
          iss >> i >> m1 >> m2 >> m3;
          means1.push_back(m1);
          means2.push_back(m2);
          means3.push_back(m3);
        }
      } else if (line.find("*") != std::string::npos) {
        ////std::cout << line.c_str() << std::endl;
        std::string delimeter = ":";
        int pos = line.find(delimeter);
        std::string identifier = line.substr(0, pos);
        line.erase(0, pos + delimeter.size());
        ////std::cout << identifier.c_str() << "   ";
        ////std::cout << line.c_str() << std::endl;
        if (identifier.find("Momentum") != std::string::npos) {
          currentBatch.p = stod(line);
        } else if (identifier.find("xPosition") != std::string::npos) {
          currentBatch.x = stod(line);
          xList.push_back(stod(line));
        } else if (identifier.find("yPosition") != std::string::npos) {
          currentBatch.y = stod(line);
          yList.push_back(stod(line));
        } else if (identifier.find("Theta") != std::string::npos) {
          currentBatch.theta = stod(line);
          thetaList.push_back(stod(line));
        } else if (identifier.find("Phi") != std::string::npos) {
          currentBatch.phi = stod(line);
          phiList.push_back(stod(line));
        } else if (identifier.find("Beta") != std::string::npos) {
          currentBatch.beta = stod(line);
          betaList.push_back(stod(line));
          currentBatch.p = 0;
        }
      }
    }
    currentBatch.mean.push_back(means1);
    currentBatch.mean.push_back(means2);
    if (particleCounts == 3) {
      currentBatch.mean.push_back(means3);
    }
    infile.close();
  } else {
  } //std::cout << "ERROR NO QUANTUM EFFICIENCY FILE" << std::endl;}
  measurements.push_back(currentBatch);
}

void InterpolateLikelihood::checkMeasurements() {
  ////std::cout << "##############" << std::endl;
  ////std::cout << "TEST" << std::endl;
  ////std::cout << "##############" << std::endl;
  ////std::cout << "Particle #:" << measurements.at(0).mean.size() << std::endl;
  ////std::cout << "Pixel #:" << measurements.at(0).mean.at(0).size() << std::endl;
  for (int i = 0; i < measurements.size(); i++) {
    std::cout << "Beta:" << measurements.at(i).beta << "  ";
    std::cout << "x:" << measurements.at(i).x << "   ";
    std::cout << "y:" << measurements.at(i).y << "   ";
    std::cout << "theta:" << measurements.at(i).theta << "   ";
    std::cout << "phi:" << measurements.at(i).phi << "   ";
    //for (int j = 0; j < measurements.at(i).mean.at(0).size(); j++) {
    //for (int k = 0; k < measurements.at(i).mean.size(); k++) {
    ////std::cout << measurements.at(i).mean[k][j] << "    ";
    //}
    std::cout << std::endl;
  }
}

void InterpolateLikelihood::readFolder() {
  TSystemDirectory dir(folderLocation.c_str(), folderLocation.c_str());
  TList *files = dir.GetListOfFiles();

  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file = (TSystemFile *)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith("txt")) {
        ////std::cout << fname.Data() << std::endl;
        this->readTextFile(folderLocation + fname.Data());
      }
    }
  }
}

double InterpolateLikelihood::calcBeta(double momentum, double particle) {
  momentum = momentum * 1000;
  double beta = momentum / sqrt(pow(momentum, 2) + pow(particle, 2));
  return beta;
}

double InterpolateLikelihood::calcMomentum(double beta, double particle) {
  double momentum = particle * beta / sqrt(1 - pow(beta, 2));
  return momentum;
}

std::vector<std::vector<std::shared_ptr<TGraphErrors>>> InterpolateLikelihood::getMeanVariableGraph(std::string variableString) {
  std::transform(variableString.begin(), variableString.end(), variableString.begin(), ::tolower);
  //Generates a graph of the mean number of photons for each pixel with resepct to a chosen variable
  //Input: Variable is the string that selects the interpolation varaible
  //Output
  std::vector<std::vector<std::shared_ptr<TGraphErrors>>> allGraphs;
  allGraphs.reserve(measurements.at(0).mean.size());

  //std::cout << "pTotal: " << measurements.at(0).mean.size() <<std::endl;
  //std::cout << "iTotal: " << measurements.size() << std::endl;
  //std::cout << "jTotal: " << measurements.at(0).mean.at(0).size() << std::endl;

  for (int p = 0; p < measurements.at(0).mean.size(); p++) {
    std::vector<std::shared_ptr<TGraphErrors>> pGraphs;
    pGraphs.reserve(measurements.at(0).mean.at(p).size());
    //std::cout << "####################################" << std::endl;
    //std::cout << "PARTICLE: " << measurements.at(0).particle.at(p) << std::endl;
    //std::cout << "####################################" << std::endl;

    for (int i = 0; i < measurements.size(); i++) //LIST OF TEXTFILES
    {
      double variableValue;

      if (variableString == std::string("beta")) {
        if (measurements.at(i).p != 0) {
          variableValue = calcBeta(measurements.at(i).p, measurements.at(i).mass.at(p));
        } else {
          variableValue = measurements.at(i).beta;
        }
      } else if (variableString == std::string("x")) {
        variableValue = measurements.at(i).x;
      } else if (variableString == std::string("y")) {
        variableValue = measurements.at(i).y;
      } else if (variableString == std::string("r")) {
        variableValue = sqrt(pow(measurements.at(i).y, 2) + pow(measurements.at(i).x, 2));
      } else if (variableString == std::string("theta")) {
        variableValue = measurements.at(i).theta;
      } else if (variableString == std::string("phi")) {
        variableValue = measurements.at(i).phi;
      }

      for (int j = 0; j < measurements.at(i).mean.at(p).size(); j++) //Number of Pixels
      {
        //std::cout << p << "   " << i << "    " << j << "    ";
        //std::cout << measurements.at(i).p << "    ";
        //std::cout << measurements.at(i).mean[p][j] << "   +/-    ";
        //std::cout << sqrt(measurements.at(i).mean[p][j]) << "     ";
        if (i == 0) {
          pGraphs.push_back(std::make_shared<TGraphErrors>((measurements.size())));
        }
        pGraphs.at(j)->SetPoint(i, variableValue, measurements.at(i).mean[p][j]);
        if (i == measurements.size() - 1) {
          pGraphs.at(j)->Sort();
          pGraphs.at(j)->Draw();
          std::string title = "XPosition " + std::to_string(measurements.at(i).x);
          title = title + " and YPosition " + std::to_string(measurements.at(i).y);
          title = title + " and Pixels " + std::to_string(j);
          title = title + " for a " + measurements.at(i).particle.at(p);
          pGraphs.at(j)->SetTitle(title.c_str());
        }
        //std::cout << std::endl;
      }
    }
    allGraphs.push_back(pGraphs);
  }
  return allGraphs;
}

std::vector<std::vector<std::shared_ptr<TGraph2D>>> InterpolateLikelihood::getMean2VariableGraph(std::string variableString) {
  std::transform(variableString.begin(), variableString.end(), variableString.begin(), ::tolower);
  std::vector<std::vector<std::shared_ptr<TGraph2D>>> allGraphs;
  allGraphs.reserve(measurements.at(0).mean.size());

  //std::cout << "pTotal: " << measurements.at(0).mean.size() <<std::endl;
  //std::cout << "iTotal: " << measurements.size() << std::endl;
  //std::cout << "jTotal: " << measurements.at(0).mean.at(0).size() << std::endl;

  for (int p = 0; p < measurements.at(0).mean.size(); p++) {
    std::vector<std::shared_ptr<TGraph2D>> pGraphs;
    pGraphs.reserve(measurements.at(0).mean.at(p).size());
    //std::cout << "####################################" << std::endl;
    //std::cout << "PARTICLE: " << measurements.at(0).particle.at(p) << std::endl;
    //std::cout << "####################################" << std::endl;
    double level = 0;
    for (int i = 0; i < measurements.size(); i++) {
      for (int j = 0; j < measurements.at(i).mean.at(p).size(); j++) {
        if (i == 0) {
          pGraphs.push_back(std::make_shared<TGraph2D>((measurements.size())));
          pGraphs.back()->SetName(std::to_string(j).c_str());
        }

        std::string variable1String;
        std::string variable2String;
        if (variableString == std::string("xy")) {
          pGraphs.at(j)->SetPoint(i, measurements.at(i).x, measurements.at(i).y, measurements.at(i).mean[p][j]);
          variable1String = std::string("XPosition");
          variable2String = std::string("yPosition");
        } else if (variableString == std::string("thetaphi")) {
          pGraphs.at(j)->SetPoint(i, cos(measurements.at(i).phi) * sin(measurements.at(i).theta), sin(measurements.at(i).phi) * sin(measurements.at(i).theta), measurements.at(i).mean[p][j]);
          variable1String = std::string("Theta");
          variable2String = std::string("Phi");
        }
        if (i == measurements.size() - 1) {
          //pGraphs.at(j)->Sort();
          pGraphs.at(j)->Draw();
          std::string title = "Measuring " + variable1String;
          title = title + " and " + variable2String;
          title = title + " and Pixels " + std::to_string(j);
          title = title + " for a " + measurements.at(i).particle.at(p);
          pGraphs.at(j)->SetTitle(title.c_str());
        }
      }
    }
    allGraphs.push_back(pGraphs);
  }
  return allGraphs;
}

void InterpolateLikelihood::sortVariables() {
  std::sort(xList.begin(), xList.end());
  std::sort(yList.begin(), yList.end());
  std::sort(thetaList.begin(), thetaList.end());
  std::sort(phiList.begin(), phiList.end());
  std::sort(betaList.begin(), betaList.end());
}

void InterpolateLikelihood::checkVector(std::vector<int> v) {
  for (int i = 0; i < v.size(); i++) {
    std::cout << i << "    " << v.at(i) << std::endl;
  }
}

void InterpolateLikelihood::checkVector(std::vector<double> v) {
  for (int i = 0; i < v.size(); i++) {
    std::cout << i << "    " << v.at(i) << std::endl;
  }
}

std::vector<std::vector<double>> InterpolateLikelihood::getAllPoints(std::vector<double> v1, std::vector<double> v2) {
  std::vector<std::vector<double>> allPoints;
  allPoints.push_back(v2);

  int counter = 0;
  std::vector<int> permSwp(v1.size(), 1);
  for (int i = 0; i < v1.size(); i++) {
    std::sort(permSwp.begin(), permSwp.end());
    permSwp.at(i) = 0;
    do {
      //std::cout << counter << "   ";
      //for (int j = 0; j < permSwp.size(); j++) {
      //std::cout << permSwp[j] << "   ";
      //}
      //std::cout << std::endl;
      counter++;
      allPoints.push_back(returnSwap(v1, v2, permSwp));

    } while (std::next_permutation(permSwp.begin(), permSwp.end()));
  }
  removeDuplicates(allPoints);
  //for (int i = 0; i < allPoints.size(); i++) {
  //std::cout << i << "   ";
  //for (int j = 0; j < allPoints[i].size(); j++) {
  //std::cout << allPoints[i][j] << "   ";
  //}
  //std::cout << std::endl;
  //}
  return allPoints;
}

std::vector<double> InterpolateLikelihood::returnSwap(std::vector<double> v1, std::vector<double> v2, std::vector<int> permSwp) {
  for (int i = 0; i < v1.size(); i++) {
    if (permSwp.at(i) == 1) {
      v1.at(i) = v2.at(i);
    }
  }
  return v1;
}

std::vector<int> InterpolateLikelihood::getIndex(std::vector<std::vector<double>> vAll) {
  std::vector<int> iList;
  for (int i = 0; i < vAll.size(); i++) {
    iList.push_back(findMeanPhotonList(vAll.at(i)));
  }
  return iList;
}

std::vector<double> InterpolateLikelihood::getDifferenceBetween2Vectors(std::vector<double> v1, std::vector<double> v2) {
  std::vector<double> diff;
  for (int i = 0; i < v1.size(); i++) {
    diff.push_back(v1.at(i) - v2.at(i));
  }
  return diff;
}

void InterpolateLikelihood::removeDuplicates(std::vector<std::vector<double>> &vec) {
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

std::vector<double> InterpolateLikelihood::interpolateMeanLinearly(std::vector<double> variables) {
  //std::cout << "Check Variables" << std::endl;
  //std::cout << "Input Variables" << std::endl;
  //checkVector(variables);

  //std::cout << "##################" << std::endl;
  //std::cout << "LINEARLY" << std::endl;
  //std::cout << "##################" << std::endl;

  std::vector<double> meanInt;
  meanInt.reserve(measurements.at(0).mean.at(0).size());

  int upperIndex = 0;
  int lowerIndex = 0;
  getVariableBoundariesList(variables, upperIndex, lowerIndex);
  std::vector<double> vUpp = getVariablesFromIndex(upperIndex);
  std::vector<double> vLow = getVariablesFromIndex(lowerIndex);
  //checkVector(vUpp);
  //checkVector(vLow);

  std::vector<double> vDiff = getDifferenceBetween2Vectors(vUpp, vLow);
  for (int i = 0; i < vDiff.size(); i++) {
    if (vDiff.at(i) == 0) {
      vDiff.at(i) = 1;
    }
  }

  std::vector<std::vector<double>> vAll = getAllPoints(vUpp, vLow);

  //std::cout << "ALL VALUES LINEARLY: " << vAll.size() << std::endl;

  std::vector<int> indexList = getIndex(vAll);
  std::vector<double> means(measurements.at(0).mean[0].size(), 0.0);

  //std::cout << indexList.size() << std::endl;

  double deltaWeight = 1;
  for (int i = 0; i < vDiff.size(); i++) {
    deltaWeight = vDiff.at(i) * deltaWeight;
  }
  for (int i = 0; i < indexList.size(); i++) {
    int fI = indexList.at(findFurthest(i, vAll));

    std::vector<double> vDelta;
    vDelta.push_back(fabs(measurements.at(fI).x - variables.at(0)));
    vDelta.push_back(fabs(measurements.at(fI).y - variables.at(1)));
    vDelta.push_back(fabs(measurements.at(fI).theta - variables.at(2)));
    vDelta.push_back(fabs(measurements.at(fI).phi - variables.at(3)));
    vDelta.push_back(fabs(measurements.at(fI).beta - variables.at(4)));

    double vDeltaProd = 1;
    for (int d = 0; d < vDelta.size(); d++) {
      if (vDelta.at(d) == 0) {
        vDelta.at(d) = 1;
      }
      vDeltaProd = vDeltaProd * vDelta.at(d);
    }
    double weight = vDeltaProd / deltaWeight;
    for (int pixel = 0; pixel < measurements.at(0).mean[0].size(); pixel++) {
      means.at(pixel) = means.at(pixel) + measurements.at(indexList[i]).mean[0][pixel] * weight;
    }
  }
  return means;
}

std::vector<std::vector<double>> InterpolateLikelihood::getSplineValuesAlongAxis(std::vector<std::vector<double>> vAll, std::vector<double> v, std::vector<bool> vSwitch) {
  int splineVariable;
  for (int i = 0; i < vSwitch.size(); i++) {
    if (vSwitch[i]) {
      splineVariable = i;
    }
  }
  //std::cout << "ALL VALUES: " << vAll.size() << std::endl;
  std::vector<std::vector<double>> meansOfInterpolatedSpline;
  for (int i = 0; i < vAll.size(); i++) {
    std::vector<int> meanIndex;
    std::vector<double> variableList;
    std::vector<std::shared_ptr<TSpline3>> tempSplines;
    getListOfVariableIndex(vSwitch, vAll[i], meanIndex, variableList);
    //std::cout << "MEAN INDEX:" << meanIndex.size() << std::endl;
    //std::cout << "VARIABLE INDEX:" << variableList.size() << std::endl;
    tempSplines = getSplines(variableList, meanIndex, "test");
    //std::cout << "SPLINE INDEX:" << tempSplines.size() << std::endl;
    std::vector<double> tempMeans;
    //std::cout << "EVAL SPLINE: " << std::endl;
    for (int k = 0; k < tempSplines.size(); k++) {
      double tMean = tempSplines[k]->Eval(v.at(splineVariable));
      //std::cout << "PIXEL: " << k << "   " << tMean << std::endl;
      if (tMean < 0) {
        tMean = 1E-20;
      }
      tempMeans.push_back(tMean);
    }
    meansOfInterpolatedSpline.push_back(tempMeans);
  }
  return meansOfInterpolatedSpline;
}

double InterpolateLikelihood::calcGaussian(double x, double mean, double sigma, double a0) {

  return a0 * exp(-std::pow(x - mean, 2) / (2 * std::pow(sigma, 2)));
}

std::vector<std::vector<double>> InterpolateLikelihood::getGaussianValuesAlongAxis(std::vector<std::vector<double>> vAll, std::vector<double> v, std::vector<bool> vSwitch) {
  int splineVariable;
  for (int i = 0; i < vSwitch.size(); i++) {
    if (vSwitch[i]) {
      splineVariable = i;
    }
  }
  //std::cout << "ALL VALUES: " << vAll.size() << std::endl;
  std::vector<std::vector<double>> meansOfInterpolatedSpline;
  for (int i = 0; i < vAll.size(); i++) {
    std::vector<int> meanIndex;
    std::vector<double> variableList;
    std::vector<TFitResultPtr> tempGaussians;
    getListOfVariableIndex(vSwitch, vAll[i], meanIndex, variableList);
    //std::cout << "MEAN INDEX:" << meanIndex.size() << std::endl;
    //std::cout << "VARIABLE INDEX:" << variableList.size() << std::endl;
    tempGaussians = getGaussian(variableList, meanIndex, "test");
    //std::cout << "SPLINE INDEX:" << tempSplines.size() << std::endl;
    std::vector<double> tempMeans;
    //std::cout << "EVAL SPLINE: " << std::endl;
    for (int k = 0; k < tempGaussians.size(); k++) {
      double p0 = tempGaussians[k]->Value(0);
      double p1 = tempGaussians[k]->Value(1);
      double p2 = tempGaussians[k]->Value(2);
      double tMean = calcGaussian(v.at(splineVariable), p1, p2, p0);
      std::cout << "PIXEL: " << k << ":   "
                << "v:" << v.at(splineVariable) << "  " << p0 << "  " << p1 << "  " << p2 << "  " << tMean << std::endl;
      if (tMean < 0) {
        tMean = 1E-20;
      }
      tempMeans.push_back(tMean);
    }
    meansOfInterpolatedSpline.push_back(tempMeans);
  }
  return meansOfInterpolatedSpline;
}

std::vector<std::vector<double>> InterpolateLikelihood::removeVariable(std::vector<std::vector<double>> vAll, std::vector<bool> vSwitch) {
  for (int i = 0; i < vAll.size(); i++) {
    for (int j = 0; j < vAll[i].size(); j++) {
      if (vSwitch[j]) {
        vAll[i][j] = 0;
      }
    }
  }
  removeDuplicates(vAll);
  return vAll;
}

std::vector<double> InterpolateLikelihood::interpolateMeansLinearlyWithSpline(std::vector<double> variables, std::vector<bool> vSwitch) {

  //std::cout << "##################" << std::endl;
  //std::cout << "SPLINE" << std::endl;
  //std::cout << "##################" << std::endl;

  //std::cout << "Check Variables" << std::endl;
  //std::cout << "Input Variables" << std::endl;
  //checkVector(variables);

  int upperIndex = 0;
  int lowerIndex = 0;
  getVariableBoundariesList(variables, upperIndex, lowerIndex);
  std::vector<double> vUpp = getVariablesFromIndex(upperIndex);
  std::vector<double> vLow = getVariablesFromIndex(lowerIndex);
  //std::cout << "Upper Variables" << std::endl;
  //checkVector(vUpp);
  //std::cout << "Lower Variables" << std::endl;
  //checkVector(vLow);

  std::vector<double> vDiff = getDifferenceBetween2Vectors(vUpp, vLow);
  for (int i = 0; i < vDiff.size(); i++) {
    if (vDiff.at(i) == 0) {
      vDiff.at(i) = 1;
    }
  }

  std::vector<std::vector<double>> vAll = getAllPoints(vUpp, vLow);
  vAll = removeVariable(vAll, vSwitch);
  for (int i = 0; i < vAll.size(); i++) {
    //std::cout << "VARIABLE ALL: " << i << std::endl;
    checkVector(vAll[i]);
  }
  std::vector<std::vector<double>> meansAlongAxis = getSplineValuesAlongAxis(vAll, variables, vSwitch);
  std::vector<double> means(measurements.at(0).mean[0].size(), 0.0);
  //std::cout << "MEANS ALONG AXIS" << std::endl;
  for (int i = 0; i < meansAlongAxis.size(); i++) {
    std::cout << i << "   ";
    for (int j = 0; j < meansAlongAxis.size(); j++) {
      std::cout << meansAlongAxis[i][j] << "   ";
    }
    std::cout << std::endl;
  }

  double deltaWeight = 1;
  for (int i = 0; i < vDiff.size(); i++) {
    if (!vSwitch[i]) {
      deltaWeight = vDiff.at(i) * deltaWeight;
    }
  }

  for (int i = 0; i < meansAlongAxis.size(); i++) {
    int fI = findFurthest(i, vAll);
    std::vector<double> vDelta;
    vDelta.push_back(fabs(vAll[fI][0] - variables.at(0)));
    vDelta.push_back(fabs(vAll[fI][1] - variables.at(1)));
    vDelta.push_back(fabs(vAll[fI][2] - variables.at(2)));
    vDelta.push_back(fabs(vAll[fI][3] - variables.at(3)));
    vDelta.push_back(fabs(vAll[fI][4] - variables.at(4)));

    double vDeltaProd = 1;
    for (int d = 0; d < vDelta.size(); d++) {
      if (vSwitch[d] || vDelta[d] == 0) {
        vDelta.at(d) = 1;
      }
      vDeltaProd = vDeltaProd * vDelta.at(d);
    }
    double weight = vDeltaProd / deltaWeight;

    for (int pixel = 0; pixel < meansAlongAxis[i].size(); pixel++) {
      means.at(pixel) = means.at(pixel) + meansAlongAxis[i][pixel] * weight;
    }
  }
  return means;
}

std::vector<double> InterpolateLikelihood::interpolateMeansLinearlyWithGaussianFit(std::vector<double> variables, std::vector<bool> vSwitch) {

  //std ::cout << "##################" << std::endl;
  //std::cout << "SPLINE" << std::endl;
  //std::cout << "##################" << std::endl;

  //std::cout << "Check Variables" << std::endl;
  //std::cout << "Input Variables" << std::endl;
  //checkVector(variables);

  int upperIndex = 0;
  int lowerIndex = 0;
  getVariableBoundariesList(variables, upperIndex, lowerIndex);
  std::vector<double> vUpp = getVariablesFromIndex(upperIndex);
  std::vector<double> vLow = getVariablesFromIndex(lowerIndex);
  //std::cout << "Upper Variables" << std::endl;
  //checkVector(vUpp);
  //std::cout << "Lower Variables" << std::endl;
  //checkVector(vLow);

  std::vector<double> vDiff = getDifferenceBetween2Vectors(vUpp, vLow);
  for (int i = 0; i < vDiff.size(); i++) {
    if (vDiff.at(i) == 0) {
      vDiff.at(i) = 1;
    }
  }

  std::vector<std::vector<double>> vAll = getAllPoints(vUpp, vLow);
  vAll = removeVariable(vAll, vSwitch);
  //for (int i = 0; i < vAll.size(); i++) {
  ////std::cout << "VARIABLE ALL: " << i << std::endl;
  //checkVector(vAll[i]);
  //}
  std::vector<std::vector<double>> meansAlongAxis = getGaussianValuesAlongAxis(vAll, variables, vSwitch);
  std::vector<double> means(measurements.at(0).mean[0].size(), 0.0);
  //std::cout << "MEANS ALONG AXIS" << std::endl;
  //for (int i = 0; i < meansAlongAxis.size(); i++) {
  //std::cout << i << "   ";
  //for (int j = 0; j < meansAlongAxis.size(); j++) {
  //std::cout << meansAlongAxis[i][j] << "   ";
  //}
  //std::cout << std::endl;
  //}

  double deltaWeight = 1;
  for (int i = 0; i < vDiff.size(); i++) {
    if (!vSwitch[i]) {
      deltaWeight = vDiff.at(i) * deltaWeight;
    }
  }

  for (int i = 0; i < meansAlongAxis.size(); i++) {
    int fI = findFurthest(i, vAll);
    std::vector<double> vDelta;
    vDelta.push_back(fabs(vAll[fI][0] - variables.at(0)));
    vDelta.push_back(fabs(vAll[fI][1] - variables.at(1)));
    vDelta.push_back(fabs(vAll[fI][2] - variables.at(2)));
    vDelta.push_back(fabs(vAll[fI][3] - variables.at(3)));
    vDelta.push_back(fabs(vAll[fI][4] - variables.at(4)));

    double vDeltaProd = 1;
    for (int d = 0; d < vDelta.size(); d++) {
      if (vSwitch[d] || vDelta[d] == 0) {
        vDelta.at(d) = 1;
      }
      vDeltaProd = vDeltaProd * vDelta.at(d);
    }
    double weight = vDeltaProd / deltaWeight;

    for (int pixel = 0; pixel < meansAlongAxis[i].size(); pixel++) {
      means.at(pixel) = means.at(pixel) + meansAlongAxis[i][pixel] * weight;
    }
  }
  return means;
}
std::vector<std::vector<double>> InterpolateLikelihood::getPointsAlongAxis(std::vector<bool> variableSwitch, std::vector<std::vector<double>> vAll) {
  std::vector<std::vector<double>> temp;
  for (int i = 0; i < vAll.size(); i++) {
    int dupCounter = 0;
    for (int j = 0; j < vAll.size(); j++) {
      for (int k = 0; k < vAll[j].size(); k++) {
        if (i != j && !variableSwitch[k] && vAll[j][k] == vAll[i][k]) {
          dupCounter++;
        }
      }
    }
    if (dupCounter == 2) {
      temp.push_back(vAll[i]);
    }
  }
  return temp;
}

int InterpolateLikelihood::findFurthest(int vI, std::vector<std::vector<double>> vAll) {
  int furthestIndex;
  double r = 0;
  for (int i = 0; i < vAll.size(); i++) {
    if (r < calcDistance(vAll.at(vI), vAll.at(i))) {
      r = calcDistance(vAll.at(vI), vAll.at(i));
      furthestIndex = i;
    }
  }
  return furthestIndex;
}

double InterpolateLikelihood::calcDistance(std::vector<double> v, std::vector<double> vp) {
  double r = 0;
  for (int i = 0; i < v.size(); i++) {
    r = r + std::pow((v.at(i) - vp.at(i)), 2);
  }
  return std::sqrt(r);
}

std::vector<double> InterpolateLikelihood::interpolateMeanToFirstOrder(std::vector<double> variables) {
  std::cout << "Check Variables" << std::endl;
  std::cout << "Input Variables" << std::endl;
  checkVector(variables);

  std::vector<double> meanInt;
  meanInt.reserve(measurements.at(0).mean.at(0).size());

  int upperIndex = 0;
  int lowerIndex = 0;

  getVariableBoundariesList(variables, upperIndex, lowerIndex);

  std::vector<double> vUpp = getVariablesFromIndex(upperIndex);
  std::vector<double> vLow = getVariablesFromIndex(lowerIndex);

  bool isEqual = true;
  for (int i = 0; i < variables.size(); i++) {
    if (vUpp.at(i) == variables.at(i) && vLow.at(i) == variables.at(i)) {
      isEqual = true && isEqual;
    } else {
      isEqual = false && isEqual;
    }
  }
  if (isEqual) {
    return measurements.at(upperIndex).mean.at(0);
  }

  std::vector<double> interpolatedMeans;

  for (int i = 0; i < measurements.at(upperIndex).mean.at(0).size(); i++) {
    double mean = 0;
    int countMeans = 0;
    for (int j = 0; j < variables.size(); j++) {
      if ((vUpp.at(j) - vLow.at(j)) != 0) {
        mean = measurements.at(upperIndex).mean.at(0).at(i) * (variables.at(j) - vLow.at(j)) / (vUpp.at(j) - vLow.at(j)) +
               measurements.at(lowerIndex).mean.at(0).at(i) * (vUpp.at(j) - variables.at(j)) / (vUpp.at(j) - vLow.at(j)) + mean;
        countMeans++;
      }
    }
    mean = mean / countMeans;
    meanInt.push_back(mean);
  }
  return meanInt;
}

void InterpolateLikelihood::getListOfVariableIndex(std::vector<bool> variableSwitch,
                                                   std::vector<double> inputVariables,
                                                   std::vector<int> &meanIndexList,
                                                   std::vector<double> &listOfVariables) {
  /*   Find the list of measurement indexes given 4 out of 5 of the variables also returns
 *   a list of variables that corresepond to the missing variable.
 *   For example: if the fixed variables are (0.06,0.06,0.30,0.56) and the variable switch is 
 *   set to x [i.e variableSwitch.at(0) = true] then it will return the indices which correspond to
 *   x from 0.06 to 0.06+side of pixel
 *
 *   Input: 
 *   VariableSwitch: bool vector that chooses the variable
 *   fixedVariables: vector of 5 doubles that remain fixed in the loop
 *   meanIndexList: pbr of the returned list of means
 *   listOfVariables: pbr of the returned list of variables
 */
  for (int i = 0; i < measurements.size(); i++) {
    std::vector<double> tempMeasurementVariables;
    tempMeasurementVariables.push_back(measurements.at(i).x);
    tempMeasurementVariables.push_back(measurements.at(i).y);
    tempMeasurementVariables.push_back(measurements.at(i).theta);
    tempMeasurementVariables.push_back(measurements.at(i).phi);
    tempMeasurementVariables.push_back(measurements.at(i).beta);
    bool isTrue = true;
    double variable = 0.0;
    for (int j = 0; j < inputVariables.size(); j++) {
      if (!variableSwitch.at(j)) {
        isTrue = isTrue && tempMeasurementVariables.at(j) == inputVariables.at(j);
      } else {
        variable = tempMeasurementVariables.at(j);
      }
    }
    if (isTrue) {
      meanIndexList.push_back(i);
      listOfVariables.push_back(variable);
    }
  }
  return;
}

std::vector<std::shared_ptr<TSpline3>> InterpolateLikelihood::getSplines(std::vector<double> independent, std::vector<int> dependent, std::string fileName = "") {
  std::string pathName;
  pathName = "./Images/";

  std::vector<std::shared_ptr<TGraph>> gTempPlot;
  std::vector<std::shared_ptr<TSpline3>> meanSplines;
  int colorCount = 0;

  //std::cout << measurements.at(0).mean.at(0).size() << std::endl;
  //std::cout << independent.size() << std::endl;
  //std::cout << "GET POINTS FOR SPLINES" << std::endl;
  for (int pixel = 0; pixel < measurements.at(0).mean.at(0).size(); pixel++) {
    std::shared_ptr<TGraph> gTemp(std::make_shared<TGraph>(independent.size()));
    //std::cout << "PIXEL:  " << pixel << std::endl;
    for (int i = 0; i < independent.size(); i++) {
      //std::cout << independent.at(i) << "   " << measurements.at(dependent.at(i)).mean.at(0).at(pixel) << std::endl;
      gTemp->SetPoint(i, independent.at(i), measurements.at(dependent.at(i)).mean.at(0).at(pixel));
    }
    meanSplines.push_back(std::make_shared<TSpline3>(std::to_string(pixel).c_str(), gTemp.get()));
  }
  return meanSplines;
}

std::vector<TFitResultPtr> InterpolateLikelihood::getGaussian(std::vector<double> independent, std::vector<int> dependent, std::string fileName = "") {
  std::string pathName;
  pathName = "./Images/";

  std::vector<std::shared_ptr<TGraph>> gTempPlot;
  std::vector<TFitResultPtr> meanGaussian;
  int colorCount = 0;

  //std::cout << measurements.at(0).mean.at(0).size() << std::endl;
  //std::cout << independent.size() << std::endl;
  //std::cout << "GET POINTS FOR Gaussian" << std::endl;
  for (int pixel = 0; pixel < measurements.at(0).mean.at(0).size(); pixel++) {
    std::shared_ptr<TGraph> gTemp(std::make_shared<TGraph>(independent.size()));
    //std::cout << "PIXEL:  " << pixel << std::endl;
    for (int i = 0; i < independent.size(); i++) {
      //std::cout << independent.at(i) << "   " << measurements.at(dependent.at(i)).mean.at(0).at(pixel) << std::endl;
      gTemp->SetPoint(i, independent.at(i), measurements.at(dependent.at(i)).mean.at(0).at(pixel));
    }
    meanGaussian.push_back(gTemp->Fit("gaus", "S"));
  }
  return meanGaussian;
}

std::vector<double> InterpolateLikelihood::getVariablesFromIndex(std::vector<int> index) {
  std::vector<double> vIn;
  vIn.push_back(xList.at(index.at(0)));
  vIn.push_back(yList.at(index.at(1)));
  vIn.push_back(thetaList.at(index.at(2)));
  vIn.push_back(phiList.at(index.at(3)));
  vIn.push_back(betaList.at(index.at(4)));
  return vIn;
}

std::vector<double> InterpolateLikelihood::getVariablesFromIndex(int index) {
  std::vector<double> vIn;
  vIn.push_back(measurements.at(index).x);
  vIn.push_back(measurements.at(index).y);
  vIn.push_back(measurements.at(index).theta);
  vIn.push_back(measurements.at(index).phi);
  vIn.push_back(measurements.at(index).beta);
  return vIn;
}

int InterpolateLikelihood::findMeanPhotonList(std::vector<double> variables) {
  for (int i = 0; i < measurements.size(); i++) {
    if (measurements.at(i).x == variables.at(0) &&
        measurements.at(i).y == variables.at(1) &&
        measurements.at(i).theta == variables.at(2) &&
        measurements.at(i).phi == variables.at(3) &&
        measurements.at(i).beta == variables.at(4)) {
      return i;
    }
  }
  return -1;
}

void InterpolateLikelihood::getVariableBoundariesList(std::vector<double> variables,
                                                      int &indexUpper,
                                                      int &indexLower) {
  bool firstUpp = true;
  bool firstLow = true;
  //std::cout << "UPPER" << std::endl;
  for (int i = 0; i < measurements.size(); i++) {
    double upX = measurements.at(i).x - variables.at(0);
    double upY = measurements.at(i).y - variables.at(1);
    double upTheta = measurements.at(i).theta - variables.at(2);
    double upPhi = measurements.at(i).phi - variables.at(3);
    double upBeta = measurements.at(i).beta - variables.at(4);
    //std::cout << i << "  ";
    //std::cout << upX << "   ";
    //std::cout << upY << "   ";
    //std::cout << upTheta << "   ";
    //std::cout << upPhi << "   ";
    //std::cout << upBeta << std::endl;
    if (upX >= 0 && upY >= 0 && upTheta >= 0 && upPhi >= 0 && upBeta >= 0) {
      indexUpper = i;
      break;
    }
  }

  //std::cout << "LOWER" << std::endl;
  for (int i = measurements.size() - 1; i > 0; i--) {
    double lowX = measurements.at(i).x - variables.at(0);
    double lowY = measurements.at(i).y - variables.at(1);
    double lowTheta = measurements.at(i).theta - variables.at(2);
    double lowPhi = measurements.at(i).phi - variables.at(3);
    double lowBeta = measurements.at(i).beta - variables.at(4);

    //std::cout << i << "  ";
    //std::cout << lowX << "   ";
    //std::cout << lowY << "   ";
    //std::cout << lowTheta << "   ";
    //std::cout << lowPhi << "   ";
    //std::cout << lowBeta << std::endl;

    if (lowX <= 0 && lowY <= 0 && lowTheta <= 0 && lowPhi <= 0 && lowBeta <= 0) {
      indexLower = i;
      break;
    }
  }

  return;
}

bool InterpolateLikelihood::isVectorPositive(std::vector<int> v) {
  for (int i = 0; i < v.size(); i++) {
    if (v.at(i) < 0) {
      return false;
    }
  }
  return true;
}

bool InterpolateLikelihood::isPositive(int v) { return (v >= 0); }

double InterpolateLikelihood::calcDelta(std::vector<double> grid) {
  //std::sort(grid.begin(),grid.end());
  double delta = grid.at(1) - grid.at(0);
  if (delta == 0) {
    return 1;
  } else {
    return delta;
  }
}

int InterpolateLikelihood::getUpperBoundaryIndex(std::vector<double> grid, double v) {
  int zeroDeltaCount = 0;
  for (int i = 0; i < grid.size(); i++) {
    if (v <= grid.at(i)) {
      if (i == 0) {
        return 1;
      }
      return i;
    }
  }

  //i.e if the variable is equal to the entire grid
  //then there is no change and the any value can work.
  return -1;
}

void InterpolateLikelihood::SaveImageOfGraphs(std::vector<std::vector<std::shared_ptr<TGraphErrors>>> grs, std::string type = "mean") {
  for (int p = 1; p < 2; p++) {
    //std::cout << "SAVE" << std::endl;
    std::string pathName = "./Images/MeanInterpolation/";
    //Create Path Name and File names
    std::string particleName = measurements.at(0).particle.at(p);
    mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    std::string fileName = pathName + particleName + type + "Spline";
    std::string namePng = fileName + ".png";
    std::string nameRoot = fileName + ".root";
    std::string namePdf = fileName + ".pdf";

    std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>(particleName.c_str(), particleName.c_str(), 1200, 1200));
    c1->DivideSquare(grs.at(p).size());
    for (int i = 0; i < grs.at(p).size(); i++) {
      c1->cd(i + 1);
      grs[p][i]->SetName(std::to_string(i + 1).c_str());
      grs[p][i]->SetMaximum(1.5);
      grs[p][i]->SetMinimum(0.0);
      grs[p][i]->Draw("AC");
      grs[p][i]->Paint("AC");
      gPad->Modified();
      gPad->Update();
    }
    c1->cd(0);
    c1->SaveAs(namePng.c_str());
    c1->SaveAs(nameRoot.c_str());
    c1->SaveAs(namePdf.c_str());
  }
  return;
}

void InterpolateLikelihood::SaveImageOfEachGraph(std::vector<std::vector<std::shared_ptr<TGraphErrors>>> grs, std::string type = "mean") {
  std::string pathName = "./Images/MeanInterpolationPixels/";
  //Create Path Name and File names
  mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  //std::string fileName = pathName+ particleName + type+"Spline";
  //std::string namePng = fileName + ".png";
  //std::string nameRoot = fileName + ".root";
  //std::string namePdf = fileName + ".pdf";

  for (int i = 0; i < grs.at(0).size(); i++) {
    mkdir((pathName + std::to_string(i)).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    for (int p = 1; p < 2; p++) //2 particles
    {
      std::string particleName = measurements.at(0).particle.at(p);

      std::string fileName = pathName + std::to_string(i) + "/" + std::to_string(i) + particleName + type;
      std::string namePng = fileName + ".png";
      std::string nameRoot = fileName + ".root";
      std::string namePdf = fileName + ".pdf";
      std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>((std::to_string(i) + particleName).c_str(),
                                                            (std::to_string(i) + particleName).c_str(),
                                                            1200, 1200));
      grs[p][i]->SetName(std::to_string(i + 1).c_str());
      grs[p][i]->Draw("AC");
      grs[p][i]->GetXaxis()->SetTitle(type.c_str());
      grs[p][i]->GetYaxis()->SetTitle("Mean Detected Photons");

      gPad->Modified();
      gPad->Update();
      c1->cd(0);
      c1->SaveAs(namePng.c_str());
      c1->SaveAs(nameRoot.c_str());
      c1->SaveAs(namePdf.c_str());
    }
  }
}

void InterpolateLikelihood::SaveImageOf2DGraphs(std::vector<std::vector<std::shared_ptr<TGraph2D>>> grs, std::string type = "mean") {
  for (int p = 0; p < 1; p++) {
    std::string pathName = "./Images/MeanInterpolation/";

    //Create Path Name and File names
    std::string particleName = measurements.at(0).particle.at(p);
    mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    pathName = pathName + type + "/";
    mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    std::string fileName = pathName + particleName + type + "2DSpline";
    std::string namePng = fileName + ".png";
    std::string nameRoot = fileName + ".root";
    std::string namePdf = fileName + ".pdf";

    std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>(particleName.c_str(), particleName.c_str(), 1200, 1200));
    c1->DivideSquare(grs.at(p).size());
    for (int i = 0; i < (grs.at(p).size()); i++) {
      c1->cd(i + 1);
      grs[p][i]->SetName(std::to_string(i + 1).c_str());
      gStyle->SetPalette(1);
      grs[p][i]->Draw("COL4Z");
      grs[p][i]->SetMaximum(1.5);
      grs[p][i]->SetMinimum(0.0);
      grs[p][i]->Paint("AC");
      grs[p][i]->GetXaxis()->SetRangeUser(measurements.front().p, measurements.back().p);
      if (type.find("TEST") != std::string::npos) {
        grs[p][i]->SetMaximum(256);
        grs[p][i]->SetMinimum(0.0);
      }
      gPad->Modified();
      gPad->Update();
    }
    c1->cd(0);
    c1->SaveAs(namePng.c_str());
    c1->SaveAs(nameRoot.c_str());
    c1->SaveAs(namePdf.c_str());
  }
  return;
}

void InterpolateLikelihood::SaveImageOfEach2DGraphs(std::vector<std::vector<std::shared_ptr<TGraph2D>>> grs, std::string type = "mean") {
  std::string pathName = "./Images/MeanInterpolationPixels/";
  mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  for (int i = 0; i < grs.at(0).size(); i++) {
    mkdir((pathName + std::to_string(i)).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    for (int p = 0; p < grs.size(); p++) //2 particles
    {
      std::string particleName = measurements.at(0).particle.at(p);

      std::string fileName = pathName + std::to_string(i) + "/" + std::to_string(i) + particleName + type;
      std::string namePng = fileName + ".png";
      std::string nameRoot = fileName + ".root";
      std::string namePdf = fileName + ".pdf";
      std::shared_ptr<TCanvas> c1(std::make_shared<TCanvas>((std::to_string(i) + particleName).c_str(),
                                                            (std::to_string(i) + particleName).c_str(),
                                                            1200, 1200));
      grs[p][i]->SetName(std::to_string(i + 1).c_str());
      grs[p][i]->SetMaximum(1.5);
      grs[p][i]->SetMinimum(0.0);
      grs[p][i]->GetXaxis()->SetTitle(type.c_str());
      grs[p][i]->GetYaxis()->SetTitle("Mean Detected Photons");

      gStyle->SetPalette(1);
      grs[p][i]->Draw("colz");
      gPad->Modified();
      gPad->Update();
      c1->cd(0);
      c1->SaveAs(namePng.c_str());
      c1->SaveAs(nameRoot.c_str());
      c1->SaveAs(namePdf.c_str());
    }
  }
}

void InterpolateLikelihood::printVariablLists() {

  std::string pathName = "./Textfiles/";
  mkdir((pathName).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  std::string fileName = pathName;

  fileName = fileName + "variableList.txt";
  int w = 20;
  std::ofstream ofile(fileName.c_str());
  if (ofile.is_open()) {
    ofile << "#x" << std::setw(w);
    ofile << "y" << std::setw(w);
    ofile << "phi" << std::setw(w);
    ofile << "theta" << std::setw(w);
    ofile << "beta"
          << "\n";

    for (int i = 0; i < xList.size(); i++) {
      ofile << xList.at(i) << std::setw(w);
      ofile << yList.at(i) << std::setw(w);
      ofile << phiList.at(i) << std::setw(w);
      ofile << thetaList.at(i) << std::setw(w);
      ofile << betaList.at(i) << "\n";
      //ofile << i << std::setw(w) << i << std::setw(w) << i << "\n";
    }
    ofile.close();
  }
}

void InterpolateLikelihood::readVariableLists() {

  std::string pathName = "./Textfiles/";
  std::string fileName = pathName;
  fileName = fileName + "variableList.txt";
  std::ifstream infile(fileName.c_str());
  if (infile.is_open()) {
    while (infile.good()) {
      measuredBatch currentBatch;
      std::string line;
      getline(infile, line);

      std::istringstream iss(line);
      getline(infile, line);
      if (line.find("#") == std::string::npos &&
          line.find("*") == std::string::npos && !line.empty()) {
        double x = 0.0;
        double y = 0.0;
        double theta = 0.0;
        double phi = 0.0;
        double beta = 0.0;
        iss >> x >> y >> phi >> theta >> beta;
        xList.push_back(x);
        currentBatch.x = x;

        yList.push_back(y);
        currentBatch.y = y;

        thetaList.push_back(theta);
        currentBatch.theta = theta;

        phiList.push_back(phi);
        currentBatch.phi = phi;

        betaList.push_back(beta);
        currentBatch.beta = beta;
      }
      measurements.push_back(currentBatch);
    }
  }
}

void InterpolateLikelihood::findVariableGivenMean(int pixel, double mean,
                                                  double &x, double &y,
                                                  double &theta, double &phi,
                                                  double &beta,
                                                  std::vector<bool> variableSwitch) {
  std::vector<double> variables;
  variables.push_back(x);
  variables.push_back(y);
  variables.push_back(theta);
  variables.push_back(phi);
  variables.push_back(beta);

  double delta = 1E9;
  int savedIndex = 0;
  for (int j = 0; j < measurements.size(); j++) {
    if (!variableSwitch.at(0) &&
        (variables.at(1) == measurements.at(j).y) &&
        (variables.at(2) == measurements.at(j).theta) &&
        (variables.at(3) == measurements.at(j).phi) &&
        (variables.at(4) == measurements.at(j).beta)) {
      double tempDelta = std::abs(mean - measurements.at(j).mean[0][pixel]);
      if (tempDelta < delta) {
        savedIndex = j;
        delta = tempDelta;
      }
    }

    else if ((variables.at(0) == measurements.at(j).x) &&
             (variables.at(2) == measurements.at(j).theta) &&
             (variables.at(3) == measurements.at(j).phi) &&
             (variables.at(4) == measurements.at(j).beta)) {
      double tempDelta = std::abs(mean - measurements.at(j).mean[0][pixel]);
      if (tempDelta < delta) {
        savedIndex = j;
        delta = tempDelta;
      }
    }

    else if ((variables.at(1) == measurements.at(j).y) &&
             (variables.at(0) == measurements.at(j).x) &&
             (variables.at(3) == measurements.at(j).phi) &&
             (variables.at(4) == measurements.at(j).beta)) {
      double tempDelta = std::abs(mean - measurements.at(j).mean[0][pixel]);
      if (tempDelta < delta) {
        savedIndex = j;
        delta = tempDelta;
      }
    } else if ((variables.at(1) == measurements.at(j).y) &&
               (variables.at(2) == measurements.at(j).theta) &&
               (variables.at(0) == measurements.at(j).x) &&
               (variables.at(4) == measurements.at(j).beta)) {
      double tempDelta = std::abs(mean - measurements.at(j).mean[0][pixel]);
      if (tempDelta < delta) {
        savedIndex = j;
        delta = tempDelta;
      }
    } else if ((variables.at(1) == measurements.at(j).y) &&
               (variables.at(2) == measurements.at(j).theta) &&
               (variables.at(3) == measurements.at(j).phi) &&
               (variables.at(0) == measurements.at(j).x)) {
      double tempDelta = std::abs(mean - measurements.at(j).mean[0][pixel]);
      if (tempDelta < delta) {
        savedIndex = j;
        delta = tempDelta;
      }
    }
  }
  x = measurements.at(savedIndex).x;
  y = measurements.at(savedIndex).y;
  theta = measurements.at(savedIndex).theta;
  phi = measurements.at(savedIndex).phi;
  beta = measurements.at(savedIndex).beta;
  phi = measurements.at(savedIndex).phi;
  beta = measurements.at(savedIndex).beta;
}

void InterpolateLikelihood::testIndexFinder(std::vector<double> variables) {
  int upperIndex = 0;
  int lowerIndex = 0;

  //checkVector(xList);

  getVariableBoundariesList(variables, upperIndex, lowerIndex);

  std::vector<double> vUpp = getVariablesFromIndex(upperIndex);
  std::vector<double> vLow = getVariablesFromIndex(lowerIndex);
}

//std::vector<double> InterpolateLikelihood::interpolationTriCubicXYTheta(std::vector<double> variables) {
//std::vector<bool> vSwitch(5, false);
//vSwitch[0] = true; //x
//vSwitch[1] = true; //y
//vSwitch[2] = true; //theta

//std::cout << "Check Variables" << std::endl;
//std::cout << "Input Variables" << std::endl;
//checkVector(variables);

//int upperIndex = 0;
//int lowerIndex = 0;
//getVariableBoundariesList(variables, upperIndex, lowerIndex);
//std::vector<double> vUpp = getVariablesFromIndex(upperIndex);
//std::vector<double> vLow = getVariablesFromIndex(lowerIndex);
//std::cout << "Upper Variables" << std::endl;
//checkVector(vUpp);
//std::cout << "Lower Variables" << std::endl;
//checkVector(vLow);

//std::vector<double> vDiff = getDifferenceBetween2Vectors(vUpp, vLow);
//for (int i = 0; i < vDiff.size(); i++) {
//if (vDiff.at(i) == 0 || vSwitch[i]) {
//vDiff.at(i) = 1;
//}
//}
//std::vector<std::vector<double>> vAll = getAllPoints(vUpp, vLow);
//vAll = removeVariable(vAll, vSwitch);
//removeDuplicates(vAll);
//std::vector<std::vector<double>> interpolateMeans;
//for (int points = 0; points < vAll.size(); points++) {
//std::vector<int> pointsIndex;
//std::vector<double> variableList;
//for (int v = 0; v < vAll[points].size(); v++) {
////std::cout << vAll[points][v] << "    ";
//}
////std::cout << std::endl;
//getListOfVariableIndex(vSwitch, vAll[points], pointsIndex, variableList);
//interpolateMeans.push_back(getTriCubicInterpolationValue(pointsIndex, variables, vSwitch));
//}

//double deltaWeight = 1;
//for (int i = 0; i < vDiff.size(); i++) {
//if (!vSwitch[i]) {
//deltaWeight = vDiff.at(i) * deltaWeight;
//}
//}

//std::vector<double> means(measurements.at(0).mean[0].size(), 0.0);
//for (int i = 0; i < interpolateMeans.size(); i++) {
//int fI = findFurthest(i, vAll);
//std::vector<double> vDelta;
//vDelta.push_back(fabs(vAll[fI][0] - variables.at(0)));
//vDelta.push_back(fabs(vAll[fI][1] - variables.at(1)));
//vDelta.push_back(fabs(vAll[fI][2] - variables.at(2)));
//vDelta.push_back(fabs(vAll[fI][3] - variables.at(3)));
//vDelta.push_back(fabs(vAll[fI][4] - variables.at(4)));

//double vDeltaProd = 1;
//for (int d = 0; d < vDelta.size(); d++) {
//if (vSwitch[d]) {
//vDelta.at(d) = 1;
//}
//vDeltaProd = vDeltaProd * vDelta.at(d);
//}

//double weight = vDeltaProd / deltaWeight;
//for (int pixel = 0; pixel < interpolateMeans[i].size(); pixel++) {
//means.at(pixel) = means.at(pixel) + interpolateMeans[i][pixel] * weight;
//}
//}
//return means;
//}

//std::vector<double> InterpolateLikelihood::getTriCubicInterpolationValue(std::vector<int> pointsIndex, std::vector<double> variables, std::vector<bool> vSwitch) {

////Check to see if the points are consistently spaced and store the bound and point information

//std::vector<std::vector<double>> pointsVariablesAll;
//std::vector<double> storedDelta(3, 0.0); //must be 3 because it is triCubic Interpolation
//std::vector<int> bound(3, 0);
//for (int point = 0; point < pointsIndex.size(); point++) {
//std::vector<double> vSwitchDouble(vSwitch.begin(), vSwitch.end());
//std::vector<double> variableForPoint;
//variableForPoint.push_back(measurements[pointsIndex[point]].x);
//variableForPoint.push_back(measurements[pointsIndex[point]].y);
//variableForPoint.push_back(measurements[pointsIndex[point]].theta);
//variableForPoint.push_back(measurements[pointsIndex[point]].phi);
//variableForPoint.push_back(measurements[pointsIndex[point]].beta);
//variableForPoint.push_back(pointsIndex[point]);

//variableForPoint.erase(variableForPoint.begin() + 3);
//variableForPoint.erase(variableForPoint.begin() + 3);

//pointsVariablesAll.push_back(variableForPoint);
//if (point > 0) {
//std::vector<double> vDiff = getDifferenceBetween2Vectors(pointsVariablesAll[point - 1], pointsVariablesAll[point]);
//for (int i = 0; i < vDiff.size(); i++) {
//bound[i]++;
//if (storedDelta[i] == 0) {
//storedDelta[i] = std::fabs(vDiff[i]);
//}
//}
//}
//}
//std::sort(pointsVariablesAll.begin(), pointsVariablesAll.end(), NColumnOnlyCmp(0));
//std::sort(pointsVariablesAll.begin(), pointsVariablesAll.end(), NColumnOnlyCmp(1));
//std::sort(pointsVariablesAll.begin(), pointsVariablesAll.end(), NColumnOnlyCmp(2));
//std::vector<double> interpolatedMeanForPixels;
//for (int pixels = 0; pixels < measurements[0].mean[0].size(); pixels++) {
//std::vector<double> meansForPixel;
//for (int points = 0; points < pointsVariablesAll.size(); points++) {
//meansForPixel.push_back(measurements[pointsVariablesAll[points].back()].mean[0][pixels]);
//}
//gte::IntpTricubic3<double> tempTriCubic(5, 5, 15, pointsVariablesAll.front()[0],
//storedDelta[0], pointsVariablesAll.front()[1], storedDelta[1], pointsVariablesAll.front()[2], storedDelta[2],
//&meansForPixel[0], false);

//double eval = tempTriCubic(variables[0], variables[1], variables[2]);
//if (eval > 0) {
//interpolatedMeanForPixels.push_back(eval);
//} else if (eval <= 0) {
//interpolatedMeanForPixels.push_back(1E-20);
//}
//}
//return interpolatedMeanForPixels;
//}

std::vector<double> InterpolateLikelihood::interpolationSwitch(std::vector<double> variable,
                                                               int intSwitch, std::vector<bool> vSwitch,
                                                               double pixelSideCount, double totalDistance,
                                                               double detectorSizex, double detectorSizey,
                                                               double pixelSideCountx, double pixelSideCounty) {
  int interpolateMeanLinearly = 0;
  int interpolateMeansLinearlyWithSpline = 1;
  int interpolateMeansLinearlyWithGaussian = 2;
  std::cout << intSwitch << "   " << interpolateMeansLinearlyWithGaussian << std::endl;

  double xSteps = 0;
  double ySteps = 0;
  double phiSteps = 0;

  calcPixelPositionSteps(pixelSideCountx,
                         detectorSizex,
                         variable[0],
                         xSteps);

  calcPixelPositionSteps(pixelSideCounty,
                         detectorSizey,
                         variable[1],
                         ySteps);
  calcPixelPhiSteps(variable[3], phiSteps);

  if (xSteps != 0 || ySteps != 0) {
    shiftBeamPosition(detectorSizex, detectorSizey,
                      pixelSideCountx, pixelSideCounty,
                      xSteps, ySteps,
                      variable);
  }

  if (phiSteps != 0) {
    shiftBeamPhi(detectorSizex, detectorSizey,
                 pixelSideCountx, pixelSideCounty,
                 phiSteps, variable, xSteps, ySteps);
  }
  std::vector<double> pixelMap;
  if (intSwitch == interpolateMeanLinearly) {
    pixelMap = this->interpolateMeanLinearly(variable);
  } else if (intSwitch == interpolateMeansLinearlyWithSpline) {
    pixelMap = this->interpolateMeansLinearlyWithSpline(variable, vSwitch);
  } else if (intSwitch == interpolateMeansLinearlyWithGaussian) {
    std::cout << "GAUSSIAN" << std::endl;
    pixelMap = this->interpolateMeansLinearlyWithGaussianFit(variable, vSwitch);
  } else {
    return pixelMap;
  }
  pixelMap = shiftPixelMapPhiPosition(pixelSideCount, phiSteps, xSteps, ySteps, pixelMap);
  return pixelMap;
}

//POSITION
void InterpolateLikelihood::calcPixelPositionSteps(double pixelSideCount,
                                                   double detectorSize,
                                                   double x,
                                                   double &xSteps) {
  if (x >= detectorSize / 2 + pixelSideCount) {
    xSteps = TMath::Floor((x - detectorSize / 2) / pixelSideCount);
  } else if (x <= detectorSize / 2) {
    xSteps = TMath::Floor((x - detectorSize / 2) / pixelSideCount);
  }
  return;
}

void InterpolateLikelihood::shiftBeamPosition(double detectorSizex, double detectorSizey,
                                              double pixelSideCountx, double pixelSideCounty,
                                              double xSteps, double ySteps,
                                              std::vector<double> &variables) {

  if ((variables[0] <= detectorSizex / 2 + pixelSideCountx) & (variables[0] >= detectorSizex / 2) &
      (variables[1] <= detectorSizey / 2 + pixelSideCounty) & (variables[1] >= detectorSizey / 2)) {
    return;
  }
  variables[0] = variables[0] - (xSteps)*pixelSideCountx;
  double delx = fabs(detectorSizex / 2 - variables[0]); //if delta is negative correct for position
  variables[0] = detectorSizex / 2 + delx;

  variables[1] = variables[1] - (ySteps)*pixelSideCounty;
  double dely = fabs(detectorSizey / 2 - variables[1]); //if delta is negative correct for position
  variables[1] = detectorSizey / 2 + dely;
  return;
}

void InterpolateLikelihood::shiftPixelMapPosition(double xSteps, double ySteps, double pixelSideCount, std::vector<double> &pixelMap) {
  std::vector<double> temp;
  double tempXSteps = (xSteps);
  double tempYSteps = (ySteps);
  for (int pixel = 0; pixel < pixelMap.size(); pixel++) {
    int shiftedPixel = pixel + pixelSideCount * (tempYSteps) - (tempXSteps);
    int startPixel = (TMath::Floor(shiftedPixel / pixelSideCount)) * pixelSideCount; // + (tempXSteps);
    int endPixel = (pixelSideCount - 1) + startPixel - (tempXSteps);
    //std::cout << "O:" << pixel << "    ";
    //std::cout << "s:" << startPixel << "    ";
    //std::cout << "n:" << shiftedPixel << "   ";
    //std::cout << "e:" << endPixel << std::endl;
    if ((shiftedPixel >= startPixel & shiftedPixel <= endPixel) &
        (shiftedPixel >= 0 & shiftedPixel <= pixelSideCount * pixelSideCount - 1)) {
      temp.push_back(pixelMap[shiftedPixel]);
    } else {
      temp.push_back(1E-20); //no signal
    }
  }
  pixelMap = temp;
  return;
}

//THETA
void InterpolateLikelihood::calcPixelThetaSteps(double totaldistance,
                                                double pixelSideCount,
                                                double theta,
                                                double &thetaSteps) {
  if (theta < 0) {
    thetaSteps = TMath::Floor(totaldistance * tan(theta) / pixelSideCount);
  }
  return;
}

void InterpolateLikelihood::shiftBeamTheta(double &theta) {
  if (theta < 0) {
    theta = -theta;
    return;
  }
  return;
}

void InterpolateLikelihood::shiftPixelMapTheta(double pixelSideCount, double xSteps, double ySteps, double thetaSteps, std::vector<double> &pixelMap) {
  if (thetaSteps != 0) {
    std::vector<double> temp;
    for (int pixel = 0; pixel < pixelMap.size(); pixel++) {
      int shifted = pixelSideCount * pixelSideCount - (pixel + 1);
      int startPixel = (TMath::Floor(shifted / pixelSideCount)) * pixelSideCount;
      int endPixel = (pixelSideCount - 1) + startPixel;
      temp.push_back(pixelMap[shifted]);
    }
    pixelMap = temp;
  }
  return;
}

//PHI
void InterpolateLikelihood::calcPixelPhiSteps(double phi, double &phiSteps) {
  phiSteps = TMath::Floor(phi / (.25 * TMath::Pi()));
}

void InterpolateLikelihood::shiftBeamPhi(double detectorSizex, double detectorSizey,
                                         double pixelSideCountx, double pixelSideCounty,
                                         double phiSteps, std::vector<double> &variable,
                                         double &xSteps, double &ySteps) {
  if (phiSteps == 1 || phiSteps == 3 || phiSteps == 5 || phiSteps == 7) {
    variable[3] = (phiSteps + 1) * (TMath::Pi() / 4) - variable[3];
  } else {
    variable[3] = variable[3] - phiSteps * (TMath::Pi() / 4);
  }
  if (phiSteps == 1 || phiSteps == 6) {
    transposeXY(variable, xSteps, ySteps);
  } else if (phiSteps == 2 || phiSteps == 5) {
    transposeXY(variable, xSteps, ySteps);
    flipVSteps(variable[1], xSteps, pixelSideCounty, detectorSizey);
  } else if (phiSteps == 3 || phiSteps == 4) {
    flipVSteps(variable[0], xSteps, pixelSideCountx, detectorSizex);
  }
  if (phiSteps > 3 && phiSteps < 8) {
    if (phiSteps == 5 || phiSteps == 6) {
      flipVSteps(variable[0], ySteps, pixelSideCountx, detectorSizex);
    } else {
      flipVSteps(variable[1], ySteps, pixelSideCounty, detectorSizey);
    }
  }
}

void InterpolateLikelihood::transposeXY(std::vector<double> &variable, double &xSteps, double &ySteps) {
  double temp = 0;
  temp = variable[0];
  variable[0] = variable[1];
  variable[1] = temp;
}

void InterpolateLikelihood::flipVSteps(double &variable, double &vSteps, double pixelSideCountv, double detectorSizev) {
  double delta = pixelSideCountv + 0.5 * detectorSizev - variable;
  if (vSteps != 0) {
    vSteps = vSteps + 1;
  } else if (vSteps == 0) {
    return;
  }
  variable = delta + 0.5 * detectorSizev;
  return;
}

std::vector<double> InterpolateLikelihood::shiftPixelMapPhiPosition(double pixelSideCount, double phiSteps, double xSteps, double ySteps, std::vector<double> pixelMap) {

  if (phiSteps == 1 || phiSteps == 6) {
    transposePixelMap(pixelSideCount, xSteps, ySteps, pixelMap);
  } else if (phiSteps == 2 || phiSteps == 5) {
    transposePixelMap(pixelSideCount, xSteps, ySteps, pixelMap);
    flipPixelMapAlongAxis(pixelSideCount, xSteps, ySteps, pixelMap, 'y');
  } else if (phiSteps == 3 || phiSteps == 4) {
    flipPixelMapAlongAxis(pixelSideCount, xSteps, ySteps, pixelMap, 'y');
  }
  if (phiSteps > 3 && phiSteps < 8) {
    flipPixelMapAlongAxis(pixelSideCount, xSteps, ySteps, pixelMap, 'x');
  }
  shiftPixelMapPosition(xSteps, ySteps, pixelSideCount, pixelMap);
  return pixelMap;
}

void InterpolateLikelihood::transposePixelMap(double pixelSideCount, double &xSteps, double &ySteps, std::vector<double> &pixelMap) {
  std::vector<double> temp;
  for (int pixel = 0; pixel < pixelMap.size(); pixel++) {
    double level = TMath::Floor(pixel / pixelSideCount);
    double diagonalPixel = (pixelSideCount * (level + 1) - 1) - level;
    double distanceFromDiagonal = diagonalPixel - pixel;
    double shiftedPixel = distanceFromDiagonal * (pixelSideCount + 1) + pixel;
    temp.push_back(pixelMap[shiftedPixel]);
  }
  pixelMap = temp;
}

void InterpolateLikelihood::flipPixelMapAlongAxis(double pixelSideCount, double &xSteps, double &ySteps, std::vector<double> &pixelMap, char c) {
  std::vector<double> temp;
  if (c == 'y') {
    for (int pixel = 0; pixel < pixelMap.size(); pixel++) {
      double level = TMath::Floor(pixel / pixelSideCount);
      double center = ((level + 1) * pixelSideCount - 1) - pixelSideCount / 2;
      double shiftedPixel = 2 * center + 1 - pixel;
      temp.push_back(pixelMap[shiftedPixel]);
    }
    pixelMap = temp;
  } else if (c == 'x') {
    double allPixels = pixelSideCount * pixelSideCount;
    for (int pixel = 0; pixel < pixelMap.size(); pixel++) {
      double level = TMath::Floor(pixel / pixelSideCount);
      double shiftedPixel = (allPixels - pixelSideCount) + 2 * (pixel - level * pixelSideCount) - pixel;
      temp.push_back(pixelMap[shiftedPixel]);
    }
    pixelMap = temp;
  }
  return;
}
