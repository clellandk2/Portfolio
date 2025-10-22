## Overview

A code I wrote to subtract out multiphonon effects from inelastic neutron scattering data from the VISION beamline at ORNL (SNS). Intended to be used when the momentum transfer (Q) is not explicitly recorded by the instrument.
If that's gibberish to you, you probably have no use for this.

Given the wide range in coding experience of ORNL users, I wrote this code **primarily with readability in mind**. I hope it can be understood so long as the reader is familiar with python basics. That said, it could be made to run faster with better data structures/memory management. If you know what you're doing, feel free to make modifications, so long as the original author is credited (see license). For my use cases (~2000 data points), the runtime was 1-2 minutes.

To adapt to other indirect geometry instruments, modify the line "E1 = 3.6" in the getW function with the appropriate final collection energy (in meV). If Q dependence is explicitly recorded by the instrument, there are better ways to do this calculation.


## Usage
The code can be run as a normal python script or installed as a package. If installed, the command is 

        sample <temperature> <material type>

which is equivalent to running  

        python VISION_Multiphonon_Correction.py <temperature> <material type>


## Notes about Input Files

### File Formatting

This code assumes that VISION data has been exported from **MANTID** via the SaveAscii algorithm. Examples can be found in the Example_Data folder. 

It is currently setup to read only from **VISION's spectrum 2** (the Q-averaged spectrum). Rewriting the scrapexy function to read other formats will work, but you should use the getDOS script if you have complete wavevector data. 

### Arguments

The arguments given to the function are:
 - temperature (temperature array, integer values in Kelvin)
 - material type (experiment descriptor)

Temperature can consist of multiple values to analyze multiple datasets and generate overplots. Material type may only be one string. Once these arguments are given, the code will search for files in the Example_Data directory with the names \<temp1\>\<material type\>SQE.txt, \<temp2\>\<material type\>SQE.txt, etc.

Below is an example command to anlyze and overplot three data files:

        sample 100 200 300 Stoich

Type 

        sample --help

for additional information.

The code does **not** include error checking or graceful handling of invalid file names so specifying nonexistant files will simply fail.  

The temperature values given as arguments are also used in mathematical calculations and should ideally be precise to ~1 K. 






## Output 

The code produces a number of graphs which are useful for various scientific purposes. The graphs are not displayed, but saved as .png files in the Example_Data directory. 





