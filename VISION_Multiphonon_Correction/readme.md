## Overview

A code I wrote to subtract out multiphonon effects from inelastic neutron scattering data from the VISION beamline at ORNL (SNS). Intended to be used when the momentum transfer (Q) is not explicitly recorded by the instrument.
If that's gibberish to you, you probably have no use for this.

Given the wide range in coding experience of ORNL users, I wrote this code **primarily with readability in mind**. I hope it can be understood so long as the reader is familiar with python basics. That said, it could be made to run faster with better data structures/memory management. If you know what you're doing, feel free to make modifications, so long as the original author is credited (see license). For my use cases (~2000 data points), the runtime was 1-2 minutes and this was unnecessary.

To adapt to other indirect geometry instruments, modify the line "E1 = 3.6" in the getW function with the appropriate final collection energy (in meV). If Q dependence is explicitly recorded by the instrument, there are better ways to do this calculation.



\### Notes about Input Files



\- \*\*File Formatting:\*\*



&nbsp; This code assumes that VISION data has been exported from \*\*MANTID\*\* via the SaveAscii algorithm. It is currently hard-coded to read only from \*\*VISION's spectrum 2\*\* (the Q-averaged spectrum).



&nbsp; Rewriting the scrapexy function to read other formats should work, but you should probably just use getDOS if it's a direct geometry instrument where you have Q data.



\- \*\*Input Files:\*\*



&nbsp; Files are specified by:



&nbsp; - tempArray (temperature array)

&nbsp; - type (experiment descriptor)



&nbsp; which are given as command line arguments.



&nbsp; Temperature values are also used in calculations and should ideally be precise to ~1 K. Files must be named:



&nbsp; (temp)(type)SQE.txt



&nbsp; The code does \*\*not\*\* include error checking or graceful handling, because a temperature typo could lead to incorrect analysis and should be quickly noticed ("fail loudly").



\- \*\*Command-line Arguments:\*\*



&nbsp; tempArray and type are passed in via command-line arguments. Example:



&nbsp; ```bash

&nbsp; python VISION\_Multiphonon\_Correction.py 100 200 300 Stoich



