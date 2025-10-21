#General Notes: This code is intended for use in analyzing neutron scattering data taken on the VISION beamline at the SNS at ORNL. With some modifications, it could be adapted to other indirect beamlines. 
#Lists are used as the primary data strucutre in an attempt to make things more easily understandable for users who may be less familiar with coding. For advanced users who require higher order MP corrections or are interested in a very large energy range, switching over to dictionaries should decrease the runtime, but would be a non trivial effort. You may also want to look into libraries for memory efficient calculation of discrete convolution products. 
#In my own experience with E<200meV, 0.5 meV binning, and MP = 3 this usually took around a minute. 

#The math done here closely follows the appendix of "Phonon density of states in vanadium" (Sears, 1995), which should be consulted for reference. 

#from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import sys
import argparse
from VISION_Multiphonon_Correction.config import VisionConfig
from matplotlib.ticker import AutoMinorLocator


#Code to change directory since MANTID defaults can be annoying (modify as appropriate):
change = 1 #0 is don't
startdir = os.getcwd() 
if (change ==1):
    #directory = startdir + "/data/Python Scripts/VISION MP/Vanadium Tests"
    directory = startdir + "/Example_Data"
    os.chdir(directory)



globalDOSmax = VisionConfig.globalvmax   #only for plotting

#Function to convert back and forth between meV and J
def convertE(value,mode):
    if (mode == 0):    #meV to J
        return (value * 1.602*10**-22)
    elif (mode == 1):  #J to meV
        return (value / (1.602*10**-22))
    else:
        print("No. Bad.")
        return 0

#Unit Conversions
kb = convertE(1.38*10**-23,1)   
h=convertE(6.626*10**-34,1)  
VisionConfig.Mass = VisionConfig.Mass * (62.42)           #The 62.42 number is a weird unit conversion
VisionConfig.globalvmin = VisionConfig.globalvmin / h     #convert to frequency 
VisionConfig.globalvmax = VisionConfig.globalvmax / h   
#globalDOSmax = globalDOSmax / h   
VisionConfig.globalgammamax = VisionConfig.globalgammamax / h   

#trapezoid rule integration
def Nintegrate(xdata,ydata):
    #gonna assume a nicely ordered array, but there is a sanity check 
    total = 0
    if (len(xdata) != len(ydata)):
        print ("Coordinate arrays are not the same size!\n")
    for i in range(len(xdata)):
        x = xdata[i]
        if (i ==0):
            prevx = x
            continue
        if (prevx >= x):
            print ("This doesn't sort x data, do it yourself (no duplicates either)\n") 
            #this might come up later but hopefully not
        area = (x - prevx)*(ydata[i]+(ydata[i-1]-ydata[i])/2)
        total += area
        if (total < 0):
            dumb = 1
        prevx = x
    return total

#reads data from a text file
def scrapexy(dataFile, constantBackground, kBT):
    #I just want averaged signal, so channel 2 (for VISION!)
    data = open(dataFile, 'r')
    xdata = []
    ydata = []           
    printAlready = 1
    channelcount=0
    for line in data:
        #print(line)
        if (len(line.split(",")) == 1):
            channelcount+=1
            continue
        if (channelcount < 3):
            continue
        if (channelcount > 3):
            print ("More than 3 channels, adjust code accordingly")
            continue
        xvalue = line.split(",")[0]
        yvalue = line.split(",")[1]
        #print (xvalue)
        #print (yvalue)
        try :
            float(xvalue)
        except: 
            print("Bad data (x) \n")
            continue    
        try :
            float(yvalue)
        except: 
            print("Bad data (y) \n")
            continue   
        xvalue = float(xvalue) / h   #conversion to frequency, important
        if (xvalue > VisionConfig.globalvmax):
            continue
        if (xvalue < VisionConfig.globalvmin):
            continue
        xdata.append(xvalue) 
        if (VisionConfig.backgroundShape == 1):  
            finalBackground = constantBackground * np.arctan(xvalue*h/kBT)   #arctan
            if (printAlready == 1):
                print("\n Using Arctan Background \n")
                printAlready = 0
        elif (VisionConfig.backgroundShape == 2):  
            finalBackground = constantBackground * np.tanh(xvalue*h/kBT)   #hyperbolic tangent
            if (printAlready == 1):
                print("\n Using Tanh Background \n")
                printAlready = 0
        elif (VisionConfig.backgroundShape == 3): 
            finalBackground = constantBackground * (xvalue * h/kBT)**(1/3) / 1000  #cubic, never got this to fit well, don't think this is the effect from VISION either
            if (printAlready == 1):
                print("\n Using Cubic Background \n")
                printAlready = 0
        else:
            finalBackground = constantBackground   #constant (not a great approximation at low energy)
            if (printAlready == 1):
                print("\n Using Constant Background \n")
                printAlready = 0
        ydata.append(float(yvalue) - finalBackground)    
    data.close()
    return xdata, ydata

#Computes Gamma0 as defined in the Sears (1995) paper this is all based on
def getgamma(vArray,gArray,v0):
    integrand = []
    xpositive = []
    vmin = VisionConfig.globalvmin 
    vmax = VisionConfig.globalgammamax   #the integral is technically to infinity, but we're assuming DOS is 0 after this point
    
    for i in range(len(vArray)): 
        v = vArray[i]
        if (v < vmin):
            continue #the integral is defined from 0 to Infinity, removing negative freq's here 
        if (v > vmax):
            continue #implements a cutoff frequency
        g = gArray[i]
        xpositive.append(v)
        number = (np.cosh(v/(2*v0))/np.sinh(v/(2*v0)))*g/v 
        integrand.append(number)
    return Nintegrate(xpositive,integrand)

#Self explanatory
def normalize (vArray, gArray):
    area = Nintegrate(vArray, gArray)
    #print ("area:", area)
    newArray = []
    for i in range(len(gArray)):
        newArray.append(gArray[i]/area)
    #print ("new area:", Nintegrate(vArray, newArray))
    return newArray  #just y data, x axis data has not changed

#Computes the A1 spectra defined in Sears (1995)
def getA1 (vArray,gArray,gamma0, v0):
    A1 = []
    for i in range(len(gArray)):
        if (round(vArray[i]*h,3) == 0):
            value = 0
        else:
            value = gArray[i]/vArray[i]
            value = value / (1-np.exp(-1*vArray[i]/v0))
            value = value/gamma0
        A1.append(value)
    return A1

#pads the x and y axis data arrays with energy values and 0s in case the arrays are missized.
#I wrote this when I only had the positive energy region and the convolution product shifted 
#the start point rightwards so it didn't line up with the previous data, which made plotting
#really annoying. I don't think it's still necessary, but I kept it in just in case. If it's not 
#doing anything, then it's not adding runtime, so I wouldn't remove it.
def add0s (x1data, xdata, A1y,A2y):    
    fixedA1 = []
    fixedA2 = []
    finalx = []
    
    x1Low = x1data[0]
    x1High = x1data[len(x1data)-1]
    x2Low = xdata[0]
    x2High = xdata[len(xdata)-1]
    #I'm assuming both x arrays are ordered from least to greatest
    if (x2Low < x1Low):
        xmin = round(x2Low,1)
    else:
        xmin = round(x1Low,1)
    #fortunately the data does keep the same binning window, but here's a check anyway
    deltax = round((xdata[1] - xdata[0]),1)
    if (deltax != round((x1data[1] - x1data[0]),1)):
        print ("Inconsistent binning, this should never happen.")
    
    if (x2High > x1High):
        xmax = round(x2High,1)
    else:
        xmax = round(x1High,1)
    x = xmin
    index1=0
    index2=0
    while (x <= xmax):
        if (x < x1Low):
            fixedA1.append(0)
        elif (x > x1High):
            fixedA1.append(0)
        else:
            fixedA1.append(A1y[index1])
            index1 +=1
        if (x < x2Low):
            fixedA2.append(0)
        elif (x > x2High):
            fixedA2.append(0)
        else:
            fixedA2.append(A2y[index2])
            index2 +=1 
        finalx.append(x)
        x+=deltax
       
    return finalx,fixedA1, fixedA2

#Computes the next A_n as defined in Sears (1995)
def getNextA(x1data, xdata, A1y, A2y):
    newx,A1y,A2y = add0s(x1data,xdata,A1y,A2y)
    if (len(A1y) != len(A2y)):
        print ("Lengths are wrong, this shouldn't happen")
    if (len(newx) != len(A2y)):
        print ("Lengths are wrong, this shouldn't happen")
    ANy = np.convolve(A1y,A2y)
    finalx = []
    deltax = round(newx[1] - newx[0],1)
    #check constant binning here
    for i in range(len(newx)-1):
        xrange = round((newx[i+1] - newx[i]),1)
        if (xrange != deltax):
            print ("Data does not have constant binning!")
            print(xrange)
    curx = round((newx[0]*2),1)
    finalx.append(curx)
    for k in range(len(ANy)-1):
        curx = round((curx + deltax),1)     #This round is necessary becuase computer hardware is imperfect, and sometimes on the SNS cluster it adds random small values to floating points that add up... 
        finalx.append(curx)             

    return finalx, ANy

#Linearly interpolates a function defined by x:(CoarseBins) and y:(MPY)
#into a new x axis binning given by FineBins, which is assumed to be denser
#BE SURE THE BINS HAVE THE SAME START AND END POINTS (the code should take care of this)
def Interp (FineBins,CoarseBins,MPY): 
    watcherror = 0
    newYbin = []
    #assumes they start in the same place
    coarseBinIndex=1
    for i in range(len(FineBins)):
        x = round(FineBins[i],1)
        cx = round(CoarseBins[coarseBinIndex],1)
        while (x > cx):
            coarseBinIndex += 1
            if (coarseBinIndex >= len(CoarseBins)):
                coarseBinIndex = len(CoarseBins) - 1
                watcherror = 1
            cx = round(CoarseBins[coarseBinIndex],1)
            #this will loop infinitely I believe if your fine bins and coarse bins don't have the same endpoint
            #so like, don't do that (and none of this will come up unless you go modify this code yourself)
        deltax = round(CoarseBins[coarseBinIndex] - CoarseBins[coarseBinIndex-1],1)    
        relx = (x-CoarseBins[coarseBinIndex-1])/deltax
        if (relx > 1.0000001):   #I use this number becuase machine precision can make this slightly bigger than 1 sometimes, without there being an actual problem 
            print ("This shouldn't be possible if data is ordered correctly.", relx)
        deltay = MPY[coarseBinIndex] - MPY[coarseBinIndex - 1]
        newYbin.append(MPY[coarseBinIndex - 1] + (relx * deltay))
    return newYbin


#computes the Debye-Waller Factor
def getW(energy, gamma0): 
    if (energy < 0):
        energy = energy * -1
    E1 = 3.6   #mev   #this is fixed for VISION, change if using another instrument
    NeutMass =  1.675 * (10**-27)    #in Kg
    NeutMass = NeutMass * (62.42)
    deg = 90
    AngleCorrection = 2*np.cos(deg * np.pi/180)
    Q = np.sqrt((8 * (np.pi**2) * NeutMass  / (h**2)  ) * (2*E1  + energy - AngleCorrection*np.sqrt(E1**2 + energy * E1)))
    vr = (Q*h/(2*np.pi))**2 / (2 * h * VisionConfig.Mass)
    W = vr * gamma0 / 2
    RMS = np.sqrt(vr*gamma0) / Q 
    #if (energy < 10):
     #   print("RMS:", RMS)   #RMS value if you want it
    return W,vr

#computes the S(Q) contribution to S(Q,E), as defined in Sears (1995)
def getSnQ(freq,n, gamma0):
    energy = freq * h
    W,notUsed = getW(energy, gamma0)
    SnQ = ((2*W)**n)/math.factorial(n)*np.exp(-2*W)
    return SnQ

#This was a function I wrote at some point to export the data as I went, but I never wound up using it
def writeIt(exportList, expfile):
    output= open(expfile, 'w')
    i=0
    while (i < len(exportList)):
        #output.write(str(int(i/2))+"\n")
        xdat = exportList[i]
        ydat = exportList[i+1]
        for j in range(len(xdat)):
            line = str(xdat[j]) + "," + str(ydat[j]) + ",0\n"
            output.write(line)
        i += 2

globalArray = []

#Computes the new DOS from SQE data and various parameters
def getNewDOS(SQEx, SQEy, gamma0,v0, DOSmax):
    DOSy = []
    mu = 0 #just in case this comes up for something... 
    for a in range(len(SQEx)):
        energy = SQEx[a]
        if (energy > DOSmax):
            DOSy.append(0)  #this is here to keep the x and y arrays the same size for plotting and such
            continue
        if (energy < DOSmax * -1):
            DOSy.append(0)
            continue
        W,vr = getW(energy, gamma0) 
        v=energy / h
        BOSE = 1-np.exp(-v/v0+mu)  #inverse of the bose factor, but I call it this anyway
        globalArray.append(BOSE)
        DOSy.append(v / vr * SQEy[a] / np.exp(-2*W) * BOSE)
    return DOSy

#Function for trimming the DOS 
def cutDOS(xdata, ydata, mode):
    xout = []
    yout = []
    threshold = VisionConfig.globalgammamax
    if (mode == 1):
        threshold = VisionConfig.globalgammamax * h
    for i in range(len(xdata)):
        x = xdata[i]
        y = ydata[i]
        if (x < -1*threshold):
            continue
        if (x < threshold):
            xout.append(x)
            yout.append(y)
    return xout,yout

#Takes a specified file, reads the data into arrays, and interpolates it into a constant, standardized binning which just makes everything so much easier
def standardizeBin(SQEFile, constantBackground, kBT):
    firstx, firsty = scrapexy(SQEFile+"SQE.txt", constantBackground, kBT)
    x = VisionConfig.globalvmin
    incriment = firstx[1] - firstx[0]
    incriment = VisionConfig.standardBinning / h
    masterxbin = []
    while (x <= firstx[len(firstx)-1]):
        masterxbin.append(x)
        x += incriment
    firsty = Interp(masterxbin, firstx, firsty) 
    
    #this is included in case you want to scale this up, but there's a different part of the code that does automatic rescaling so probably unnecessary
    for index in range(len(firsty)):
        firsty[index] = 1.0* firsty[index]
    
    return masterxbin, firsty
    
#Mirrors the positive energy data to presumptive negative values 
#Detailed balance can either be included or not (Should be included for SQE, but not for DOS)
def mirror(xdata, ydata,v0, DB):
    fullx = []
    fully = []
    i = len(xdata) - 1
    while (i >= 1):
        fullx.append(-1*xdata[i])
        if (DB == 1):
            corrected = np.exp(-1 * xdata[i] / v0) * ydata[i]  #detailed balance 
                #xdata better be in units of frequency
        else: 
            corrected = ydata[i]
        fully.append(corrected)
        i = i - 1
        
    xmin = -1 * xdata[0]
    xmax = xdata[0]
    if (xmax <= xmin):
        print ("Weird. Orignal binning and globalvmin was negative probably (It's assumed positive).") 
        return
    deltax = xdata[1] - xdata[0]
    if (deltax <= 0):
        print ("Data not ordered by increasing x.")
        return
    qensFillx = []
    qensFillx.append(xmin)
    while (round(xmin*h,1) < round(xmax*h, 1)):   
        xmin += deltax
        qensFillx.append(xmin)
    if (round(xmin*h,1) != round(xmax*h,1)):
        print ("Binning has to be constant, and you have to choose cutoffs so that 2*globalvmin is divisible by the x incriment")
        return 
        #this will not work for all choices of binning, the whole point of this is an identical x incriment everywhere
    if (DB == 0):
        zeropoint = len(qensFillx)//2
        firstx = qensFillx[:zeropoint]
        secondx = qensFillx[zeropoint:]
        qensFilly1 = Interp(firstx, [-1 * xdata[0], 0],[ydata[0],0])
        qensFilly2 = Interp(secondx,[0, xdata[0]],[0,ydata[0]])
        qensFilly = qensFilly1 + qensFilly2
        #this fancyness is so the DOS goes to 0 at 0
        #All the divide by 0 errors are handled by automatically setting the result to 0 
    else:
        qensFilly = Interp(qensFillx, [-1 * xdata[0],xdata[0]], [(np.exp(-1 * xdata[i] / v0) * ydata[i]), ydata[0]])
        
    qensFillx.pop(len(qensFillx) - 1)  #this is accounted for in xdata
    qensFilly.pop(len(qensFilly) - 1)
    
    fullx = fullx + qensFillx + xdata
    fully = fully + qensFilly + ydata
    return fullx, fully

#this is the big function that processes a SQE spectra and handles the MP corrections and everything
def processSQE(xdata, ySQE, v0, SQEFile, bigList, origy):         
    xdatac = []  #xdata converted (back to meV)
    for x in range(len(xdata)):
        xdatac.append(xdata[x] * h)
    savex = xdata
    savexc = xdatac
    
    gamma0 = 2.89*10**(-13)  #initial guess for gamma0, feel free to change but it should be about this order of magnitude
    
    xdatac, cutSQE = cutDOS(xdatac, ySQE, 1)
    xdata, unused = cutDOS(xdata, ySQE, 0)
    tempy = getNewDOS(xdatac, cutSQE, gamma0,v0, VisionConfig.globalgammamax)
    tempy = normalize(xdata, tempy)
   
    #gamma0 depends on DOS, DOS depends on gamma0, gotta do a cycle...
    convergencethreshold = 0.01
    giveup = 7  #max number of iterations  
    plt.clf()
    plottingDOSconv = "yes"
    print ("\nInitiating gamma0/DOS self consistency cycle...")
    print ("convergencethreshold =", convergencethreshold, ", giveup =", giveup, ", plottingDOSconv =", plottingDOSconv)
    for i in range(giveup):
        if (plottingDOSconv != "no"):
            tempcolors = ["Blue", "Green", "Gold" , "Orange" , "Red", "Black", "Brown", "Pink"]
            plt.plot(xdatac,tempy, color = tempcolors[i], label = i)
        gammaprev = gamma0
        gamma0=getgamma(xdata,tempy,v0)
        print("\nCycle", i, "; Gamma0 =", gamma0)
        if (abs((gamma0 - gammaprev)/gammaprev) < convergencethreshold):
            break
        elif (i == giveup - 1):
            print ("Gamma0 didn't converge, look into this further. Calculations may not be valid.")
        else:
            tempy = getNewDOS(xdatac, cutSQE, gamma0,v0, VisionConfig.globalgammamax)
            tempy = normalize(xdata, tempy)

    xdata = savex 
    xdatac = savexc 
    
    if (plottingDOSconv != "no"):
        print ("Plotting convergence...")
        plt.xlabel('Energy (meV)')
        plt.ylabel('DOS')
        plt.title("DOS/Gamma0 Convergence")
        plt.legend()
        plt.savefig('DOS Convergence')
        #note: this will overwrite itself if doing multiple files in the same run
        #alternatively, if you open the file up beforehand you can see it cycle though each of the convergence cycles
        
    print("Final Gamma0:", gamma0)
    ydata = getNewDOS(xdatac, ySQE, gamma0,v0, VisionConfig.globalgammamax)
    ydata = normalize(xdata, ydata)
    A1y = getA1(xdata,ydata, gamma0, v0)
    A1y = normalize(xdatac,A1y)
    
    #The normalization produced by the Sears (1995) equations does not match what comes out of the VISION algorithms
    #Likely, this is because the VISION data is probably normalized to 1 including all the strong noise at high energies (and NOT including the negative energy range). 
    #Olivier thought the normalization should match experiment and not just be the product of equations though (too idealized)
    #So to match the data to the initial scaling, the correction factor is computed here:

    #a = len(A1y) // 8 * 5
    a = len(A1y) // 2 + 40      #this is a completely arbitrary choice of point, the only thing that's important is it's non zero
    SnQ = getSnQ(xdata[a],1, gamma0)
    scalefactor = ySQE[a] / (A1y[a] * SnQ * 1000) #the 1000 is becuase I'm doing normalizations in meV not eV
    
    if (VisionConfig.InheritScale == "No"):
        scalefactor = 1
    print("\nRescalling factor:", scalefactor)
    print("\n")

    plt.clf()
    #plt.xticks(np.arange(-300,300,40))
    #plt.yticks(np.arange(0,0.05,0.01))
    #plt.xlim((-300,300))
    #plt.ylim((0,0.05))
    plt.xlabel('Energy (meV)')
    plt.ylabel('Intensity')
    plt.title("An's")

    colors = ["Blue", "Green", "Orange", "Red", "Purple", "Brown"]
    test = 0
    if (test == 1):
        testx,testy = getNextA([1,2,3,4,5],[1,2,3,4,5],[1,3,6,4,2])
        print(testx)
        print(testy)

    j=2
    newx = xdata
    newy = A1y
    total = {}  #sorry the names got confusing, "Total" is the data minus MP corrections
    sum = {}   #"sum" is the data plus MP corrections to assess convergence

    n=1
    SnQ = 1 
    SqList = []
    exportList = []    #this stores both x and y arrays
    for a in range(len(A1y)):
        #A1 is normalized
        SnQ = getSnQ(xdata[a],1, gamma0)
        SqList.append(A1y[a] * SnQ * 1000 * scalefactor) #the 1000 is becuase I'm doing normalizations in meV not eV
        if round(xdatac[a],1) in total:
            total[round(xdatac[a],1)] = total[round(xdatac[a],1)] + origy[a]
            sum[round(xdatac[a],1)] = sum[round(xdatac[a],1)] + SqList[a]
            print ("Shouldn't ever happen?") 
        else:
            total[round(xdatac[a],1)] = origy[a]
            sum[round(xdatac[a],1)] = SqList[a]

    exportList.append(xdatac)
    exportList.append(SqList)

    plt.plot(xdatac,A1y, color = "Indigo", label="A1")

    while (j <= VisionConfig.MP):
        SqList = []
        newx,newy = getNextA(xdata,newx,A1y,newy)
        newxc = []
        for b in range(len(newx)):
            newxc.append(newx[b] * h)
        newyc = normalize(newxc,newy) 
        for a in range(len(newyc)):
            SnQ = getSnQ(newx[a],j, gamma0)
            SqList.append(newyc[a] * SnQ * 1000 * scalefactor)  #1000 is becuase I'm normalizing in meV not eV
        for i in range(len(newxc)):
            x = newxc[i]
            if round(x,1) in total:
                total[round(x,1)] = total[round(x,1)] - SqList[i]
                sum[round(x,1)] = sum[round(x,1)] + SqList[i]
            else:
                total[round(x,1)] = 0 #zeroing MP contribution outside the graphing range 
                sum[round(x,1)] = 0
        dataTitle = "A" + str(j)
        plt.plot(newxc,newyc, color = colors[j-2], label=dataTitle)
        exportList.append(newxc)
        exportList.append(SqList)
        j += 1
    newx = newxc
    totalx = []
    totaly = []
    sumy = []
    z=0
    while (z < len(xdatac)):
        xtrack = round(xdatac[z],1)
        z += 1             
        totalx.append(xtrack)
        totaly.append(total[xtrack])
        sumy.append(sum[xtrack])
        
    origx = xdatac  
    plt.legend()
    plt.savefig(SQEFile+'A_nGraph')

    plt.clf()

    #plt.xticks(np.arange(0,300,40))
    #plt.yticks(np.arange(0,0.02,0.005))
    plt.xlim((-VisionConfig.globalvmax*h,VisionConfig.globalvmax*h))
    #plt.ylim((0,0.26))
    plt.xlabel('Energy (meV)')
    plt.ylabel('SQE correction')
    plt.title("SQE's")
   
    #this is all just plotting here
    j = 0
    plt.plot(origx,ySQE, color = "Indigo", label="Input Data") 
        
    while (j < len(exportList)):
        plotx = exportList[j]
        ploty = exportList[j+1]
        theN = "SQE_" + str(int(j/2 + 1))
        plt.plot(plotx,ploty, color = colors[int(j/2)], label=theN)
        j += 2 
        
    plt.legend()
    plt.savefig(SQEFile + 'S_nComponents') 

    plt.clf()

    #plt.xticks(np.arange(-300,300,40))
    #plt.yticks(np.arange(0,0.02,0.005))
    #plt.xlim((0,300))
    #plt.ylim((0,0.02))
    plt.xlabel('Energy (meV)')
    plt.ylabel('SQE')
    plt.title("SQE (Corrected)")
        
    plt.plot(origx,origy, color = "Black", label="Original SQE")
    plt.plot(totalx,totaly, color = "Green", label = "Estimated One Phonon")
    plt.plot(totalx,sumy, color = "Blue", label = "Sum of MP Contributions")
    
    plt.legend()
    plt.savefig(SQEFile + 'S_nSubtracted') 
    bigList.append(totalx)
    bigList.append(totaly)

    #Convergence criteria for the multiphonon iterations 
    showareas = 0    #you can print the three areas from the "Subtracted" graphs if you set this equal to 0
    area1 = Nintegrate(origx,origy)
    area3 = Nintegrate(totalx,sumy)
    convratio = area3 / area1
    if (showareas == 0):
        print ("CONVERGENCE CHECK!!")
        print ("Areas:")
        print (area1)
        #print (Nintegrate(totalx,totaly))
        print (area3)
        print("ratio:", convratio)
        print ("\n")
    stopiterations = 0
    if (convratio > VisionConfig.MPConvBounds[0]):
        if (convratio < VisionConfig.MPConvBounds[1]):
            stopiterations = 1
    
    DOSy = getNewDOS(origx,totaly, gamma0,v0, VisionConfig.globalvmax) 
    DOSy = normalize(origx,DOSy)
    plt.clf()
    #plt.xticks(np.arange(0,200,40))
    #plt.yticks(np.arange(0,0.015,0.005))
    plt.xlim((-1*globalDOSmax,globalDOSmax))
    #plt.ylim((0,0.015))
    plt.xlabel('Energy (meV)')
    plt.ylabel('DOS (normalized)')
    plt.title("Density of States")
    
    ydata = normalize(xdatac, ydata)
    
    plt.plot(xdatac,ydata, color = "Green", label = "Original")
    plt.plot(origx,DOSy, color = "Blue", label="Recalculated Phonon DOS")
    bigList.append(DOSy)   #(the x coordinates are the same)

    plt.legend()
    plt.savefig(SQEFile + 'DOSGraph') 

    print ("Done")
    return stopiterations 

#This makes some nice SQE and DOS overplots, as well as some less nice ones 
#For the nice SQE Graph I interpolated through the elastic peak to remove the vertical line at E = 0
def overplotSQEDOS(bigList, tempArray, type, Iteration):
    j = 0 

    colors = ["Blue", "Green", "Gold", "Orange", "Red", "Purple", "Brown"]
    #colors = ["Blue", "Gold", "Orange", "Red"]
    markers = ["o", "^", "s", "X", "d"]  #I guess the plotting stuff will crash if you try more than 5 temperatures, just add more to the arrays or do a fancier form of plotting
    plt.clf()
    #plt.xticks(np.arange(0,200,40))
    #plt.yticks(np.arange(0,0.015,0.005))
    #plt.xlim((0,200))
    #plt.ylim((0,0.015))
    plt.xlabel('Energy (meV)')
    plt.ylabel('SQEs')
    plt.title("SQE")
    while (j < len(bigList)):
        plt.plot(bigList[j],bigList[j+1], color = colors[int(j/3)], label= tempArray[int(j/3)] +"K SQE")
        j += 3
    plt.legend()
    plt.savefig(type + " SQE Overplot" + str(Iteration)) 
    #old lazy SQE plot
    
    j = 0 
    plt.clf()
    #plt.xticks(np.arange(0,200,40))
    #plt.yticks(np.arange(0,0.015,0.005))
    plt.xlim(0,globalDOSmax)
    #plt.ylim((0,0.015))
    plt.xlabel('Energy (meV)')
    plt.ylabel('DOS')
    plt.title("Density of States")
    while (j < len(bigList)):
        plt.plot(bigList[j],bigList[j+2], color = colors[int(j/3)], label= tempArray[int(j/3)] +"K DOS")
        #ax.errorbar(bigList[j],bigList[j+2], fmt = 'o', color = colors[int(j/3)], label= tempArray[int(j/3)] +"K DOS")
        
        j += 3
    
    plt.legend()
    plt.savefig(type + " DOS Overplot" + str(Iteration)) 
    #old lazy DOS plot
    
    plt.clf()
    plt.rcParams['font.size'] = 22
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111)

    j = 0
    while (j < len(bigList)):
        #plt.plot(bigList[j],bigList[j+2], color = colors[int(j/3)], label= tempArray[int(j/3)] +"K DOS")
        ax.plot(bigList[j],bigList[j+2], marker = markers[int(j/3)], markersize = 5, color = colors[int(j/3)], label= tempArray[int(j/3)] +"K DOS")
        
        j += 3
    ax.set_xlim(0, VisionConfig.globalvmax * h)
    #ax.set_xlim(0,25)
    #ax.set_ylim(-0.0005,0.004)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.text(200, 0.05, 'Density of States', weight='bold')
    ax.set_xlabel('Energy (meV)')
    ax.set_ylabel(r'Density of States (meV$^{-1}$)')
    plt.axhline(y=0, color = 'black', linestyle = '-')
    ax.legend(frameon=False)
    ax.tick_params(direction='in', which='both')
    plt.tight_layout()
    savestr = 'NiceDOSOverplot' + str(Iteration)
    plt.savefig(savestr)
    #Nice DOS Plot
    
   
    
    plt.clf()
    plt.rcParams['font.size'] = 22
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111)

    j = 0
    while (j < len(bigList)):
        #SQE Elastic Peak interpolation 
        xplot = bigList[j]
        yplot = bigList[j+1]
        yFinalPlot = []
        elasticX = []
        elasticY = []
        a = 0
        while (a < len(xplot)):
            energy = xplot[a] 
            if (energy < -1 * VisionConfig.globalvmin * h):
                yFinalPlot.append(yplot[a])
            elif (energy <= VisionConfig.globalvmin * h): 
                elasticX.append(energy)
                elasticY.append(yplot[a])
            else:
                break
            a += 1
        elasticY = Interp (elasticX,[elasticX[0],elasticX[-1]],[elasticY[0],elasticY[-1]]) 
        for b in elasticY:
            yFinalPlot.append(b)
        while (a < len(xplot)):
            energy = xplot[a] 
            if (energy <= VisionConfig.globalvmin * h):
                print("Data is not well ordered!")
            else: 
                yFinalPlot.append(yplot[a])
            a += 1
        
        #plt.plot(bigList[j],bigList[j+2], color = colors[int(j/3)], label= tempArray[int(j/3)] +"K DOS")
        ax.plot(xplot,yFinalPlot, marker = markers[int(j/3)], markersize = 5, color = colors[int(j/3)], label= tempArray[int(j/3)] +"K DOS")
        
        j += 3
    ax.set_xlim(-1 * VisionConfig.globalvmax * h, VisionConfig.globalvmax * h)
    #ax.set_xlim(0,25)
    #ax.set_ylim(-0.0005,0.004)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.text(200, 0.05, 'S(Q,E)', weight='bold')
    ax.set_xlabel('Energy (meV)')
    ax.set_ylabel('SQE')
    plt.axhline(y=0, color = 'black', linestyle = '-')
    ax.legend(frameon=False)
    ax.tick_params(direction='in', which='both')
    plt.tight_layout()
    savestr = 'NiceSQEOverplot' + str(Iteration)
    plt.savefig(savestr)
    #Nice SQE PLot
    
    return (Iteration + 1)

#If you ever want the final result as data instead of a graph
def writeIteration(thisData, iteration): 
    return 0   #I'm not exporting the data as text, but I included the skeleton for that here if you wanna do it 
    print ("Writing files...")
    
    print("Code this!") 
    for i in range(len(thisData)):
        if (i%3 ==0):
            energyData = thisData[i]
            print ("This is the energy (x axis) data at some temperature")
        if (i%3 ==1):
            SQEData = thisData[i]
            print ("This is the SQE y axis data at the same temperature")
        if (i%3 ==2):
            DOSData = thisData[i]
            print ("This is the DOS y axis data at the same temperature")
    
    print ("Done") 

#and this is the main function, which maintains the array of files and calls all the appropriate processing functions
def iterateSQEs(tempArray, type):
    print(os.getcwd())  #useful sanity check sometimes

    #initialize lists to save the SQEs and DOS
    bigList = []    #list of lists containing x data, ySQE data and yDOS data for each temperature
    origList = []   #Original Data, as opposed to bigList, which updates
    IterCount = 1   #counts iterations
    
    #keep track of which temperatures have converged, and which need subsequent iterations
    skipArray = []
    for i in tempArray:
        skipArray.append(0)

    giveup = 4 #max number of iterations  
    for k in range(giveup):
        for i in range(len(tempArray)):
            skipIt = skipArray[i]
            temp = tempArray[i]
            SQEFile = temp+type  
            print("File:", SQEFile, "SQE.txt")
            temp = int(temp) 
            v0 = kb*temp/h
            #print ("v0:", v0)
            if (k == 0):
                constantBackground = VisionConfig.background[i]
                #constantBackground = 0
                xdata, ySQE = standardizeBin(SQEFile, constantBackground, kb*temp)
                xdata, ySQE = mirror(xdata, ySQE, v0,1)
                origy = ySQE
                origList.append(ySQE)
            elif(skipIt == 1):
                print("Skipped! Already Converged!")
                xdata = bigList.pop(0)  
                ySQE = bigList.pop(0)
                keeptheDOs = bigList.pop(0)
                bigList.append(xdata)
                bigList.append(ySQE)
                bigList.append(keeptheDOs)
                continue
            else: 
                xdata = bigList.pop(0)  
                ySQE = bigList.pop(0)
                bigList.pop(0)
                
                origy = origList[i]
                
                SQEFile = SQEFile + str(k+1) 
                for e in range(len(xdata)):
                    xdata[e] = xdata[e] / h
            alldone = processSQE(xdata, ySQE, v0, SQEFile, bigList, origy)   #it just takes the file name to name more files
            if (alldone == 1):
                skipArray[i] = 1
        
        IterCount = overplotSQEDOS(bigList, tempArray, type, IterCount)
        writeIteration(bigList, k+1)        
        print ("Iteration #", (k+1))
        
    print("All Done")


def setupSQEs():
    parser = argparse.ArgumentParser(
        description="Performs a multiphonon correction on neutron scattering SQE data.", 
        epilog="Example usage:\npython VISION_Multiphonon_Correction.py 100 200 300 Stoich"
    )
    parser.add_argument(
        "temps",
        nargs="+", 
        help="List of temperatures, e.g., 100 200 300"
    )
    parser.add_argument(
        "type",
        help="Sample type descriptor, e.g., 'Stoich' or 'Na_doped'"
    )

    args = parser.parse_args()

    tempArray = args.temps
    sample_type = args.type

    iterateSQEs(tempArray, sample_type)

if __name__ == "__main__":
    setupSQEs()
