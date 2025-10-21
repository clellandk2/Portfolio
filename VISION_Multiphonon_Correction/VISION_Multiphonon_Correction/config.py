"""
These are the important material parameters, energy cutoffs, and other settings needed to compute the multiphonon contrubution and adjust the neutron data accordingly.  
"""
class VisionConfig:
    # Average particle mass (in kg)
    Mass = 4.40*(10**-26)                         # NASICON x = 2.1 sample
    # Mass = 8.46*(10**-26)                       # Vanadium is typically used as a callibration sample

    # Number of multiphonon corrections to be included (e.g. two phonon scattering, three phonon scattering...), 3-6 is recommended 
    MP = 3  

    MPConvBounds = [0.985, 1.01]  
    # The code determines whether to perform a subsequent iteration of recalculating the one phonon spectra based on the ratio of the area of the original SQE spectra and the area of the new calculated SQE spectra based on the current estimate of the one phonon spectra. These are the thresholds for that ratio which define convergence. 

    # energy cutoffs (in meV) 
    globalvmin = 1               # chosen to remove the elastic peak
    globalvmax = 200             # maximum plot range
    globalgammamax = 140         # Where the vibrational density of states (DOS) *should* go to 0, it probably will level out to some constant in VISION data though becuase of the additional background signal (see below).  

    # Energy binning - The code automatically interpolates the data into energy bins of this size: 
    standardBinning = 0.5    # make sure your cutoff values are nicely divisible by this, I mainly use 0.5 and 0.25

    # Additional background subtraction  (becuase VISION has a lot of time independent background that changes based on the sample and isn't necesarily neatly subtracted out with the vanadium background)

    backgroundShape = 2         # determines the shape of the additional background subtraction; 1-arctan, 2-tanh, 3-x^(1/3) (Some beam scientist said to try this one time), 0 or anything else - constant

    background = [0.1,0.1,0.1,0.1,0.1]    # Scale factor of the background for each temperature. Choose Parameters that make the background go to 0 in the high energy limit for each SQE/DOS spectra

    InheritScale = "Yes"   # determines whether the code preserves the arbitrary intensity scaling in the input files. Strongly recomended for ease of comparing files. 
