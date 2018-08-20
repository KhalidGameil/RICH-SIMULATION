ing Imaging Cherenkov Detector Simulation code breakdown

This code is meant to simulate the processes involved in the detection of light emitted by aerogel by a multi-anode **PMT**. The code works by:

 - Simulating a beam of particles (momentum, position, and the angle of the beam can be adjusted) 
 - A detector is also simulated which consists of:
	 1. A multianode PMT or a set of a multianode PMT
	 2. Two pieces of aerogel used to create the cherenkov light
 - An event is considered to be when a single particle impinges the detector 
	 * This means that even if no light is emitted (i.e the particle is not greater than the cherenkov threshold) it is still considered an event
 
 ## Classes List
 - **PositionClass**: Determines position of all objects in the code including the position of the beam and the detector components (aerogel and mPMT)
 - **DirectionClass**: Determines the angle of all objects in the code 
 - **BeamClass**: Generates particles for each event, can set positions, angles, and momentum 
 - **AerogelClass**: This is where the particle creates the photons within each aerogel and where the event class is created.
 - **RandomEventClass**: The class that contains the number of photons produced from a single particle. The event contains:
	 - Number of photons (poisson distributed)
	 - Cherenkov Angle <Vector>
	 - Interaction positions  <Vector>: the point at which the photon is created
	 - Theta (should be relabelled as Phi)  <Vector>: The random angle (uniform distributed in 0-2$\pi$)
	 - Photon Energy <Vector>: uniformly distributed energy for each photon
 - **DetectedEventClass**: A class that represents when a photon is detected by a PMT has the same form as the RandomEventClass
 - **PMTClass**: Sets the size, position, and number of pixels of a mPMT. Also applies a crosstalk mask to the PMT to simulate cross talk and also reads the quatumEfficiency of the detector from the quantumEfficiency.csv file
 - **DetectorClass**: Class used to simulate the detector components. Used to arrange the geometry of the detector and simulate the analysis that goes on in the detector.
	 - Very complicated read separate readme for the DetectorClass

**OTHER:** The graphFunctionsNamespace.h contains all the graphing code used to produced plots of the histogram. 

## Analysis and Use Case

The purpose of the code is to use the simulated likelihood distribution of the different events to seperate between $\pi$ and $K$ events that will be truly measured. 
Steps of Simulation (for single particle):
 1.  The particle will generate photons at a certain cherenkov angle and a random interaction position (z position), phi angle, and energy. 
 2. These are then detected randomly (based off the quantum efficiency) by the PMT
 3. You will have a distribution of pixels that are either hit **OR** not hit 
 4. Based off of this distribution of hits (and the known momentum, position, and angle) the likelihood of the $\pi$ or $K$ hypothesis is tested

### Likelihood and Seperation value

The likelihood is constructed with the hypothesis that a pixel of the mPMT can either have a probability of a photon detection or a probability that the photon is not detected, per pixel.  
$$
P_{Detected} (\mu_i) = ln(1-e^{-\mu_i}) \qquad N_{Detected} \geq 1
$$
$$
P_{Not\ Detected} (\mu_i) = ln(e^{-\mu_i}) \qquad N_{Detected} < 1
$$
The equation above describes these probabilities where $\mu_i$ is the mean photoelectrons of that pixel (i) for a given particle. The log likelihood is then a sum of multiple detection probabilities for a particle, as shown in the equation below.

$$
ln(L(\mu_{particle})) = -2 \sum_{i=0}^{pixels} ln(P(\mu_i,N_{Detected~or~Not}))
$$

This can be done per beam event, which is a description of a fixed particle with known momentum and position. This becomes a powerful tool in estimating beam parameters (such as momentum and position) as well as determining the probability of separation of the detector.

An example of determining the probability of separation can be performed for a beam of 5GeV pions as compared to a beam of 5GeV kaons perpendicular to the mPMT. This can be calculated by generating a histogram of the difference of log-likelihood of the pion events for the pion photoelectron means and kaon photoelectron means. The same can be done for the kaon events and a second histogram can be generated, the histograms are shown below. 

$$
-2ln(\frac{L_{\pi}}{L_{K}}) = -2ln(L(\mu_{\pi})) + 2ln(L(\mu_{K})) 
$$



## Important notes

 - List item

## Good Slide Shows
- w Parameter explanation and likelihood: Nuprisim_4new.pdf
- 
