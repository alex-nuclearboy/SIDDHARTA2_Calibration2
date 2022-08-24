# The SIDDHARTA-2 Calibration

The SIDDHARTA-2 experiment at the DAΦNE collider aims to perform the first kaonic deuterium X-ray transitions to the fundamental, 1s, level measurement. 

A kaonic atom is formed when a negatively charged kaon enters a target, stops inside and replaces an electron in an atom, forming a highly excited (mass of the kaon much larger than the electron one) kaonic atom. 
The kaonic atom de-excites towards lower energy levels and emits radiation. 
When reaching the lowest level, 1s, the kaon interacts with the nucleus also by the strong interaction, which induces a shift and a broadening of the width of the level with respect to the purely electromagnetic calculated values. 
The X-ray transitions to the fundamental level alow to extract the shift and width and these quantities are then related to low-energy QCD.

To achieve this challenging goal the experimental apparatus is equipped with 384 Silicon Drift Detectors (SDDs) distributed around its cryogenic gaseous target. 
The SDDs developed by the SIDDHARTA-2 collaboration are suitable for high precision kaonic atoms spectroscopy, thanks to their high energy and time resolutions combined with their radiation hardness. 
The energy response of each detector must be calibrated and monitored to keep the systematic error, due to processes such as gain fluctuations, at the level of 2-3 eV. 

## The Silicon Drift Detectors for the SIDDHARTA-2 experiment

A monolithic SDDs array consists of eight square SDD cells, each with an active area of 8×8 mm2. 
The 450 µm thick silicon bulk allows to achieve an efficiency of almost 100% for Xrays in the energy between 5-12 keV, which corresponds to the region of interest for the kaonic deuterium measurement. 
The individual SDD cells are arranged in a 2×4 array with a 1 mm dead region along the device’s borders. 
The silicon wafer is glued on an alumina ceramic carrier, which provides the polarization voltages, and screwed on an aluminium holder. 
A very significant improvement is the change of the preamplifier system from the J-FET to a complementary metal-oxide semiconductor integrated charge sensing amplifier, named CUBE. 
For each SDD cell, the CUBE is placed on the ceramic carrier as close as possible to the collecting anode and connected to it with a bonding wire. 
Thanks to the new preamplifier CUBE, the SDDs’ performances are stable even when exposed to high and variable charged particle rates and a faster drift time is achievable thanks to the possibility to cool down to 100 K. 

The signals are processed by a dedicated front-end electronics and DAQ system. 
The output of the CUBE is connected to a common ASIC called SFERA (SDDs Front-End Readout ASIC), an Integrated Circuit performing analog shaping and peak detection of the signals. 

##  Calibration method

The SDDs are placed around the cryogenic target cell made of high purity aluminium structure and 75 µm thick Kapton walls. 
The energy calibration of the SDDs is performed using two X-ray tubes installed on two sides of the vacuum chamber, and a multi element target made of high purity titanium and copper strips placed on the target cell walls. 
The X-ray tubes induce the fluorescence emission of the target elements and the characteristic Kα and Kβ transitions are detected by the SDDs. 
The titanium and copper lines are clearly visible together with the Mn Kα, Fe Kα and Zn Kα lines produced by the accidental excitation of other components of the setup. 
The energy response function of the detector is predominantly a Gaussian curve for every fluorescence X-ray peak; however the response has a low energy component due to the incomplete charge collection and electron-hole recombination. 
Thus, the total peak fit function is formed by two contributions: 

    Gauss function (the main contribution to the peak shape);
    
    Tail function (an exponential function to reproduce the incomplete charge collection).

Only the Ti Kα and Cu Kα peaks are exploited to calibrate the detectors, since they have the highest signal-to-background ratio. 
The Fe Kα peak, not included in the determination of calibration parameters, can be used to evaluate the goodness of the calibration itself. 
After the sum of the calibrated spectra for all the SDDs, a fit was performed to obtain the energy position of the Fe Kα peak.
