# mLS

## Introduction

mLS is a function that estimates the landslide magnitude (mLS) from an array of areas derived from a landslide inventory
following the methods of Tanyas and others (in prep). The function plots the best power-law fit for the medium
and large landslides and the frequency-area distribution of the analyzed inventory. It also returns the corresponding
cutoff (smallest area that follows power law), beta (power-law exponent), and mLS (landslide magnitude) values.


*Disclaimer:* This software is preliminary or provisional and is subject to 
revision. It is being provided to meet the need for timely best science. The 
software has not received final approval by the U.S. Geological Survey (USGS).
No warranty, expressed or implied, is made by the USGS or the U.S. Government as
to the functionality of the software and related material nor shall the fact of
release constitute any such warranty. The software is provided on the condition
that neither the USGS nor the U.S. Government shall be held liable for any
damages resulting from the authorized or unauthorized use of the software. 


## Installation and Dependencies

### To install this function:

Download mLS.m and add it to your matlab path


### This package also requires the following packages: 

#### plfit

mLS uses the plfit function of Clauset et al.(2009) to determine both the power-law exponent (beta)
and the cutoff value of the power-law. The plfit function can be downloaded from the following link: 
[http://www.santafe.edu/~aaronc/powerlaws/](http://www.santafe.edu/~aaronc/powerlaws/)

Add the plfit function to your path

#### EzyFit

mLS also uses a curve fitting toolbox, EzyFit. The 
referred toolbox needs to be downloaded from the following link: 
[http://www.fast.u-psud.fr/ezyfit/](http://www.fast.u-psud.fr/ezyfit/) and follow the given steps to intall
the toolbox (the explanation given below was taken from the cited website):
   1. Download and unzip the EzyFit Toolbox in a directory somewhere in your system. For instance, in a Windows installation, the directory Documents/MATLAB/ezyfit may be a good location. Do NOT install the toolbox in the Matlab directory itself (Program Files/Matlab directory in Windows). 
   2. Select 'Set Path' (available in the menu File in Matlab 7, or in the tab Home in Matlab 8). In the dialog box, click on 'Add Folder' (NOT 'with subfolders') and select the ezyfit directory. Click on 'Save' and 'Close'.


## Usage example

```
Insert usage example code here using the sample data
```