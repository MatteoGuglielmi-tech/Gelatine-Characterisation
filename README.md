  
# Ultrasound techniques for medical applications assignment

<!--toc:start-->
- [Ultrasound techniques for medical applications assignment](#ultrasound-techniques-for-medical-applications-assignment)
  - [Problem Statement](#problem-statement)
  - [Results](#results)
<!--toc:end-->

## Problem Statement
In this specific repository, you can find a possible solution to characterize the acoustic properties of an unknown medium. The experiments were conducted inside the [ULTRa Laboratory](https://sites.google.com/view/drlibertariodemi/ultrasound-lab) in University of Trento, Italy. In particular, the code is meant to carry out the following requests :
- compute an estimation of the acoustic SOS in the medium
- estimate the attenuation coefficient of the gelatine as a function of depth and frequency

## Results
The estimated speed of sound, averaged across 5 files to smooth out some noise, turned out to be $1514 m/s$ with a variance of $~23m/s$.  
Concerning the coefficient attenuation, the most suitable model turned out to be the one with an exponential kernel. The fitted curve tells that with the increase in frequency, the attenuation coefficient tends to increase which is reasonable since the wavelengths of sound at higher frequency are shorter, i.e. the sound waves have a higher number of cycles in a given distance. This makes it easier for the sound waves to be absorbed or scattered by the particles in the medium, leading to a greater amount of attenuation.  

Further experiments show that the SOS evolves around the $1500 [m/s]$ as previously mentioned. For the attenuation coefficient, the results turned out to be very noisy across the different recorded samples.
