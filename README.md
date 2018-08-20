# Ring Imaging Cherenkov Detector Simulation code breakdown

This code is meant to simulate the processes involved in the detection of light emitted by aerogel by a multi-anode **PMT**. The code works by:

 - Simulating a beam of particles (momentum, position, and the angle of the beam can be adjusted) 
 - A detector is also simulated which consists of:
	 1. A multianode PMT or a set of a multianode PMT
	 2. Two pieces of aerogel used to create the Cherenkov light
 - An event is considered to be when a single particle impinges the detector 
	 * This means that even if no light is emitted (i.e the particle is not greater than the Cherenkov threshold) it is still considered an event
 
 ## Classes List
 - **PositionClass**: Determines the position of all objects in the code including the position of the beam and the detector components (aerogel and mPMT)
 - **DirectionClass**: Determines the angle of all objects in the code 
 - **BeamClass**: Generates particles for each event, can set positions, angles, and momentum 
 - **AerogelClass**: This is where the particle creates the photons within each aerogel and where the event class is created.
 - **RandomEventClass**: The class that contains the number of photons produced from a single particle. The event contains:
	 - Number of photons (Poisson distributed)
	 - Cherenkov Angle <Vector>
	 - Interaction positions  <Vector>: the point at which the photon is created
	 - Theta (should be relabelled as Phi)  <Vector>: The random angle (uniformly distributed in 0-2<img src="https://latex.codecogs.com/gif.latex?\pi"/>)
	 - Photon Energy <Vector>: uniformly distributed energy for each photon
 - **DetectedEventClass**: A class that represents when a photon is detected by a PMT has the same form as the RandomEventClass
 - **PMTClass**: Sets the size, position, and the number of pixels of a mPMT. Also applies a crosstalk mask to the PMT to simulate crosstalk and also reads the quantum efficiency of the detector from the quantumEfficiency.csv file
 - **DetectorClass**: Class used to simulate the detector components. Used to arrange the geometry of the detector and simulate the analysis that goes on in the detector.
	 - Very complicated read separate readme for the DetectorClass

**OTHER:** The graphFunctionsNamespace.h contains all the graphing code used to produced plots of the histogram. 

## Analysis and Use Case

The purpose of the code is to use the simulated likelihood distribution of the different events to separate between <img src="https://latex.codecogs.com/gif.latex?\pi"/> and <<img src="https://latex.codecogs.com/gif.latex?K"/> events that will be truly measured. 
Steps of Simulation (for a single particle):
 1.  The particle will generate photons at a certain Cherenkov angle and a random interaction position (z position), phi angle, and energy. 
 2. These are then detected randomly (based off the quantum efficiency) by the PMT
 3. You will have a distribution of pixels that are either hit **OR** not hit 
 4. Based off of this distribution of hits (and the known momentum, position, and angle) the likelihood of the <img src="https://latex.codecogs.com/gif.latex?\pi"/> or <img src="https://latex.codecogs.com/gif.latex?K"/> hypothesis is tested

### Likelihood and Separation value

The likelihood is constructed with the hypothesis that a pixel of the mPMT can either have a probability of a photon detection or a probability that the photon is not detected, per pixel.  

<p align="center"><img align="middle" src="https://latex.codecogs.com/svg.latex?P_{Detected}(\mu_i)=ln(1-e^{-\mu_i})~~~N_{Detected}\geq1"/></p>
<p align="center"><img align="middle" src="https://latex.codecogs.com/svg.latex?P_{Not~Detected}(\mu_i)=ln(e^{-\mu_i})~~~N_{Detected}<1"/>
</p>
The equation above describes these probabilities where $\mu_i$ is the mean photoelectrons of that pixel (i) for a given particle. The log likelihood is then a sum of multiple detection probabilities for a particle, as shown in the equation below.

<p align="center"><img align="middle" src="https://latex.codecogs.com/svg.latex?\ln{L(\mu_{particle})}=-2\sum_{i=0}^{pixels}ln(P(\mu_i,N_{Detected~or~Not}))"/></p>

This can be done per beam event, which is a description of a fixed particle with a known momentum and position. This becomes a powerful tool in estimating beam parameters (such as momentum and position) as well as determining the probability of separation of the detector.

An example of determining the probability of separation can be performed for a beam of 5GeV pions as compared to a beam of 5GeV kaons perpendicular to the mPMT. This can be calculated by generating a histogram of the difference of log-likelihood of the pion events for the pion photoelectron means and kaon photoelectron means. The same can be done for the kaon events and a second histogram can be generated, the histograms are shown below. 
<p align="center"><img src="https://latex.codecogs.com/svg.latex?-2ln(\frac{L_{\pi}}{L_{K}})=-2ln(L(\mu_{\pi}))+2ln(L(\mu_{K}))"/></p>

<p align="center"><img align="middle" src="https://github.com/KhalidGameil/RICH-SIMULATION/blob/master/docImages/SepVal.png" width="50%"></p>

From these generated histograms the separation factor, which is defined as the overlap between the kaon and pion histograms at the $90\%$ area mark for the pion histogram, can be calculated. The separation value can give a metric to analyze the particle identification as a function of momentum (or other beam parameters). The figure below shows an analysis of the separation factor, for a simulated pion beam and a simulated kaon beam, as a function of momentum.

<p align="center"><img align="middle" src="https://github.com/KhalidGameil/RICH-SIMULATION/blob/master/docImages/separationMomentum.png" width="50%"></p>

## Important notes

 - List item

## Good Slide Shows
- w Parameter explanation and likelihood: Nuprisim_4new.pdf
- 
