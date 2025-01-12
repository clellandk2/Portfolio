A code I wrote to subtract out multiphonon effects from inelastic neutron scattering data from the VISION beamline at ORNL (SNS). Intended to be used when Q transfer is not explicitely recorded by the instrument.
If that's gibberish to you, you probably have no use for this. 

Given the wide range in coding experience of ORNL users, I wrote this code primarily with readability in mind. I hope it can be understood so long as the reader is familiar with python basics. 
But that said, it could be made to run faster with beter data structures/memory management. If you know what you're doing, feel free to make modifications, so long as the original author is creditted (see liscence). For my use cases (~2000 data points), the run time was 1-2 minutes and this was unnecessary. 

To adapt to other indirect geometry instruments, modify the line "E1 = 3.6" in the getW function with the appropriate final collection energy (in meV). If Q dependence is explicitely recorded by the instrument, there are better ways to do this calculation. 
