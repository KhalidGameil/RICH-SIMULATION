# Ring Imaging Cherenkov Detector Simulation code breakdown

This code is meant to simulate the processes involved in the detection of light emitted by aerogel by a multi-anode **PMT**. The code works by:

 - Simulating a beam of particles (momentum, position, and the angle of the beam can be adjusted) 
 - A detector is also simulated which consists of:
	 1. A multianode PMT or a set of a multi-anode PMT
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
The equation above describes these probabilities where <img src="https://latex.codecogs.com/gif.latex?\pi"/> and <<img src="https://latex.codecogs.com/gif.latex?\mu_i"/> is the mean photoelectrons of that pixel (i) for a given particle. The log likelihood is then a sum of multiple detection probabilities for a particle, as shown in the equation below.

<p align="center"><img align="middle" src="https://latex.codecogs.com/svg.latex?\ln{L(\mu_{particle})}=-2\sum_{i=0}^{pixels}ln(P(\mu_i,N_{Detected~or~Not}))"/></p>

This can be done per beam event, which is a description of a particle with a known momentum and position. This becomes a powerful tool in estimating beam parameters (such as momentum and position) as well as determining the probability of separation of the detector.

An example of determining the probability of separation can be performed for a beam of 5GeV pions as compared to a beam of 5GeV kaons perpendicular to the mPMT. This can be calculated by generating a histogram of the difference of log-likelihood of the pion events for the pion photoelectron means and kaon photoelectron means. The same can be done for the kaon events and a second histogram can be generated, the histograms are shown below. 
<p align="center"><img src="https://latex.codecogs.com/svg.latex?-2ln(\frac{L_{\pi}}{L_{K}})=-2ln(L(\mu_{\pi}))+2ln(L(\mu_{K}))"/></p>

<p align="center"><img align="middle" src="https://github.com/KhalidGameil/RICH-SIMULATION/blob/master/docImages/SepVal.png" width="50%"></p>

From these generated histograms the separation factor, which is defined as the overlap between the kaon and pion histograms at the $90\%$ area mark for the pion histogram, can be calculated. The separation value can give a metric to analyze the particle identification as a function of momentum (or other beam parameters). The figure below shows an analysis of the separation factor, for a simulated pion beam and a simulated kaon beam, as a function of momentum.

<p align="center"><img align="middle" src="https://github.com/KhalidGameil/RICH-SIMULATION/blob/master/docImages/separationMomentum.png" width="50%"></p>

## Mean Number of Photon Interpolator
To speed up the analysis it was proposed that the mean number of photons be simulated for a wide range of scenarios and then act as a look-up table. This is the purpose of the "MeanGeneratorCode". The logic of the code is very simple:
1. Simulate multiple events for a given <img src="https://latex.codecogs.com/svg.latex?(\beta,x,y,\theta,\phi)"/>
2. Generate a pixel map of the mean number of photons generated
<p align="center"><img align="middle" src="https://github.com/KhalidGameil/RICH-SIMULATION/blob/master/docImages/full-13.png" width="50%"></p>

  3. Save the <img src="https://latex.codecogs.com/svg.latex?(\beta,x,y,\theta,\phi)"/> and the mean # of photons for each pixel (256 pixels) 
  4. Repeat for all other iterations of <img src="https://latex.codecogs.com/svg.latex?(\beta,x,y,\theta,\phi)"/>

### The ranges used for the interpolation
* x = 0.097 to 0.1030625m (5 steps)
* y = 0.097 to 0.1030625m (5 steps)
* <img src="https://latex.codecogs.com/svg.latex?\theta"/> = 0.00 to 0.40 (30)
* <img src="https://latex.codecogs.com/svg.latex?\phi"/> = 0.00 to <img src="https://latex.codecogs.com/svg.latex?\pi/4"/> (30)
* <img src="https://latex.codecogs.com/svg.latex?\beta"/> = 0.988 to 0.999 (20) and 0.999 to 1.0 (15) 
   * These two ranges were chosen because at higher beta values were found to vary highly between in mean number of photons
Below shows the mean number of photons as a function of <img src="https://latex.codecogs.com/svg.latex?\beta"/>:
<p align="center"><img align="middle" src="https://github.com/KhalidGameil/RICH-SIMULATION/blob/master/docImages/MeanInterpolationDataBank_database_2018Jan26.png" width="50%"></p>

After this was done the file was created and placed in the directory: `/neut/datasrv2a/kgameil/MeanInterpolationDataBank_database_2018Jan26.txt`

## W parameter: 
The w parameter was added to account for two things:
1. The generated means do not accurately reflect the number of photons that there will be in the edges
of the detector
2. The interpolation process does not account for correlations between x,y, and theta
<p align="center"><img src="https://latex.codecogs.com/svg.latex?w=\mu-\mu_0"/></p>
Where <img src="https://latex.codecogs.com/svg.latex?\mu"/> is the true mean and <img src="https://latex.codecogs.com/svg.latex?\mu_0"/> is the interpolated mean.

This has a direct effect on the log likelihood for a single particle. 
<p align="center"><img align="middle" src="https://github.com/KhalidGameil/RICH-SIMULATION/blob/master/docImages/equations.png" width="70%"></p>

If it is assumed that the minimum likelihood is at 0 then their are two possible equations for w:
<p align="center"><img src="https://latex.codecogs.com/svg.latex?w_{Hit}=\frac{\sigma_w^2e^{-\mu-w}}{e^{-\mu-w}-1}"/></p>
<p align="center"><img src="https://latex.codecogs.com/svg.latex?w_{No~Hit}=\sigma_w^2"/></p>

### Testing the W parameter:
The effectiveness of the w parameter was tested by using a minimization of the variables <img src="https://latex.codecogs.com/svg.latex?(\beta,x,y,\theta,\phi)"/> values for a given configuration. By viewing the distribution of the fitted values compared to the known true value (implemented in the code) the determination of the w parameter's effectiveness can be proven. 

|   | W Parameter OFF  | W Parameter ON  |
|---|---|---|
| Linear Interpolate Mean number of Photons (default) | <img align="middle" src="https://github.com/KhalidGameil/RICH-SIMULATION/blob/master/docImages/PionwParameterOFF_LINEAR_INTERPOLATION_10000_Events.png" width="100%">  | <img align="middle" src="https://github.com/KhalidGameil/RICH-SIMULATION/blob/master/docImages/PionwParameterON_LINEAR_INTERPOLATION_10000_Events.png" width="100%">  |
| Spline Interpolate Mean number of Photons  | <img align="middle" src="https://github.com/KhalidGameil/RICH-SIMULATION/blob/master/docImages/PionwParameterOFF_SPLINE_INTERPOLATION_10000_Events.png" width="100%">  | <img align="middle" src="https://github.com/KhalidGameil/RICH-SIMULATION/blob/master/docImages/PionwParameterON_SPLINE_INTERPOLATION_10000_Events.png" width="100%">  |

The above plots show a sample of the w parameter test. Two things are being tested in the above plots: the w parameter and the interpolation method of the mean number of photons (linear vs spline). It can be seen that the combination of the linear interpolation results in the most Gaussian solution without any false minima or residual effects from the interpolation. The test was performed in a location which requires interpolation between all five parameters.

## Main Files in analysis_wparameter folder
1. main_betaScan.txt: the main file to produce SCAN graphs of the beta minimization from Minuit
2. main_betaSweep.txt: produces loglikelihood histogram for different beta values 
3. main_meanVsbeta.txt: produces beta value vs mean Number of photon list
4. main_parameter_histogram_sweep.txt: produces files of fitted histogram parameters as well as whether minuit returned any errors or warnings

# analysis_SeperationParameter
Tests the separation parameter as a function of <img src="https://latex.codecogs.com/svg.latex?(\beta,x,y,\theta,\phi)"/>
