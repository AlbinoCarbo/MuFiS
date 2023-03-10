                                          NOTES FOR THE PRACTICAL USE OF MuFiS
                                                Albino Carbognani, Ph.D.
                                           Version 3.00.00 dated May 12, 2020
INTRODUCTION

MuFiS stands for "Multipurpose Fireball Software".

a) This software has been developed during the years 2017-2020 for the physical analysis of the bolides observed in the Earth's atmosphere within the PRISMA Project (http://www.prisma.inaf.it/).
b) It is a software that runs under 64-bit Windows inside a DOS or PowerShell window, developed and compiled with Matlab R2015b.
c) There are two modes of use:
   
1) BASIC: it is used for a first analysis of the fireball, when there are visual observations or when the speed data is noisy
2) ADVANCED: it is used when the observations are digital and with low uncertainty or when a physical model of the meteoroid is required (mass, dimensions, ablation)
and the delimitation on the ground of the strewn field.

BASIC mode is the default. It is possible to obtain the complete triangulation of the trajectory (altitude and velocity as a function of time, absolute and apparent magnitude as a function of time, coordinates of the apparent radiant, etc.), a first reconstruction of the heliocentric orbit with uncertainty estimate made using Monte Carlo and, if the final kinematic velocity is below of 6 km/s, a kinematic model of the dark flight with relative estimate of the point of fall to the ground of any residual meteoroid. As a kinematic model for velocity as a function of altitude Ceplecha model is used. If the latter model does not converge, we switch to a simpler linear model, suitable for small and fast racing cars that disintegrate at high altitudes. The same preliminary kinematic model, even if with different functionalities, is used to initialize the calculations of the dynamic model of the fireball. The isotherm law of the atmospheres, an exponential with a scale height of 7.64 km, is used as an atmospheric model for the dark flight in basic mode. With only the kinematic data the position of the strewn field is very uncertain, even by a factor of 2, therefore it is not recommended to go looking for any meteorites on the basis of kinematic analysis alone.

In the ADVANCED mode MuFiS always performs the triangulation as in the BASIC mode, but this is followed by the dynamic analysis of the motion of the bolide in the atmosphere to obtain speed to infinity, m/S ratio (from which, assuming an average density of 3500 kg/cm^3, mass and dimensions are found separately), ablation coefficient and velocity at the initial point.
This last parameter is very important because while the initial altitude is always well determined, the speed is not, therefore it must be left as a parameter to be determined within certain limits. COESA (1976) is used as atmospheric model for the dynamic model. The drag coefficient is set by the user (for mathematical reasons it cannot be separated from m/S), by default it is 0.58.
To estimate the uncertainties of the best fit parameters of the dynamic model, a Monte Carlo is used which creates 100 altitude and speed scenarios compatible with the observations.
The user can choose to raise this number. Up to 1000 Monte Carlo cycles the computation time is about 15 minutes on an i7 processor. If the final speed is below 6 km/s, MuFiS does a basic analysis of the dark flight by estimating the fall zone of the probable residual meteoroid. In the dark flight model the input data comes from the dynamic model and no longer from the kinematic model of the fireball, therefore the strewn field is more accurate. However, the state is not taken into account in dark flight actual atmosphere, such as density, pressure, temperature, and winds even if the COESA model is used. The resulting strewn field is only indicative. For more accurate calculations PRISMA uses another software that can better model the dark flight phase. Finally the heliocentric osculating orbit is calculated, but using the velocity at infinity and the uncertainties of the dynamical model.

c) The "IT20000000" folder is supplied with MuFiS which contains the files with the input data for a simulated fireball observed from 4 fictitious stations.
d) Files are text in ANSI format. Other formats could generate reading errors in the data input phase. For the structure of the input files, Standard or PRISMA, see further on.
e) This simulated fireball, whose parameters are known a priori, was used to test MuFiS.


PRACTICAL USE OF MuFiS with the simulated fireball

1-Copy the "IT20000000" folder to any location on the HDD. Inside this folder is the "Data" folder which contains:
 
   a) The files "IT20000000_n.txt" (n = 1, 2, 3, 4), contain 5 columns with the observations of time and position of the fireball seen from the n-th station with the following order:
     Time in JD; RA (??) to J2000; Dec (??) to J2000; Azimuth (??) to date; Height (??) at the date.

   b) The file "Stazioni_IT20000000.txt" contains the data on the observing stations, namely:
      Lat. (??) Long. (??) Height (m) Progressive integer % Any comments

   c) The IT20000000 folder also contains the "Results" folder with the results you should find if MuFis works correctly.
  
   NOTE 1: The order of the data columns in these files cannot be reversed, otherwise it will result in nonsense!
   NOTE 2: The first two stations of the "Stations" file are used for a preliminary triangulation using only two observers with the Ceplecha method.
           For this reason it is good that the first two stations are those that have observed the fireball in the best possible way (long and high trajectory on the horizon).
The subsequent triangulation with N > 2 stations occurs with the Borovicka method using the preliminary triangulation as a starting point.
  
   c) The IT20000000 folder also contains the "Results" folder with the results you should find if TrFis works correctly.
      Any other fireball for which input data is available can be analysed, provided that the format is Standard or PRISMA.
      In the files with the input data all the comments can be inserted after the % symbol at the beginning of the line: in this way Matlab does not read them.
  
2-Install the MATLAB Runtime R2015b libraries on the PC, downloadable from the MathWorks website:
   http://www.mathworks.com/products/compiler/mcr/index.html

3-Copy MuFiS.exe to any location on the HDD.

4-Start MuFiS running inside a Windows PowerShell session with the command: MuFiS.exe <Ret>.
The header will appear:

> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
> %                                                                     %
> %      MuFiS - Multipurpose Fireball Software  - PRISMA PROJECT       %
> %                                                                     %
> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
> %                                                                     %
> %                      by Albino Carbognani (INAF-OAS)                %
> %                                Ver 3.00.00                          %
> %                                                                     %
> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

5-Answer the questions that are asked progressively:

   a1-What is the name of the fireball folder (e.g. IT20180822)? IT20000000 <ret>
  
   a2-Input data format (PRISMA or standard, P/S)? S (on data input format see below)
  
The user is warned that there are two calculation possibilities, basic or advanced:

> BASIC computations: fireball triangulation, heliocentric orbit and dark flight (if applicable).
> Orbit and dark flight from kinematics.
> ADVANCED computations: fireball triangulation, dynamical model, simple dark flight (if applicable) and heliocentric orbit.
> Orbit and dark flight from dynamical model.

b1-Basic computations? Y/N (default: Y)
      If you answer 'Y' you have the basic analysis, ie the triangulation of the trajectory of the fireball and the orbit.
      This mode is the default one.
  
   b2-Basic computations? Y/N (default: Y)
      If you answer 'N' you get the complete analysis of triangulation, atmospheric dynamics, simplified dark flight/strewn field and orbit.
  
   c-Data processing mode in the meteoroid_dynamic function: delete or smooth (default: delete)? delete <ret>
  
   "mode = delete" (default) if you want to delete points scattered at the beginning or end of the observed trajectory.
         In this case it asks:
         1-Starting delay to avoid noisy data (fraction of the total flight time, default 0.05)? 0 <ret>
         2-Ending delay to avoid noisy data (fraction of the total flight time, default 0.05)? 0 <ret>
This is the default mode, the one that gives the best results. If you put 0 and 0 you get all the data.
  
   "mode = smooth", if you want to obtain an average trajectory by binning the data at 0.1 s intervals.
        This mode is good for "noisy" position curves and speeds. In this second mode, the software also calculates the acceleration observed in operation
        of the time and show/save the relative figures.
  
   d-Fixed drag coefficient for the fireball phase (default: 0.58)? 0.55 <ret>
  
   e-Full path of the fireball folder (e.g. C:\Users\Name\Documents\MATLAB\Fireball_data\)? ....... <ret>

   f-INFO: Creation of the folders to host results about trajectory, dynamics, dark flight and orbit. ')
     disp(' They are created inside fireball folder and if they don''t exist only.
  
6-The software automatically creates the following folders within IT20000000 (if they already exist, skip this step):

   a-"Trajectory", will contain the results of the geometric triangulation from N stations with the Borovicka method.
   b-"Dynamic", will contain the results of the atmospheric dynamics model solving the equations of motion for atmospheric drag and single body ablation.
   c-"Dark flight", will contain data on the trajectory and point of impact on the ground of any meteorite but the calculations are made with a simplified model that does not take into account the real state of the atmosphere.
   d-"Orbit", will contain data on the heliocentric orbit of the bolide with the uncertainties given by a Monte Carlo simulation (default N = 100 orbits).

7-Inside the "Trajectory", "Dynamic", "Dark Flight" and "Orbit" folders both the results files in .txt and .kml format and the figures of the graphics in .bmp format are saved.

8-The following are saved in the "Trajectory" folder:
   a) The .bmp figure with the plot of the geocentric trajectory in 3D
   b) The .bmp figures with the height and speed of the fireball as a function of time.
   c) The .bmp figure with the speed of the bolide as a function of the altitude.
   d) The .bmp figure with the absolute magnitude of the bolide as a function of time and the apparent magnitude of the bolide seen from station 1 as a function of time (only for the PRISMA type input format).
   e) A figure with the fit of the speed of the bolide as a function of the altitude with the empirical model of Ceplecha (1961). The value of the speed at infinity given by this
      simplified model allows to calculate the orbit in the "basic" mode. In the case of the "advanced" however, the speed at infinity is given by the dynamic model
of the meteoroid.
   f) Three .kml files to visualize the trajectory at the Ceplecha and the one at N stations, both on the ground and in the air, with Google Earth.
   g) Three .txt files reporting the results of the triangulation at the Ceplecha at two stations, the triangulation at N stations and the triangulation at N stations ordered in increasing time order.

9-Before starting the calculation of the dynamic model, MuFiS asks for the number of Monte Carlo scenarios to generate for the estimation of the uncertainties.
   By default it is 100 but it is advisable to put at least 1000.
  
10-The following are saved in the "Dynamic" folder:
   a) All the .bmp figures characterizing the single-body atmospheric model.
   b) The "Summary Dynamics" .txt file which contains all the numerical parameters of the dynamic model of the meteoroid.
   c) The "Dynamics" .txt file, which contains the input parameters for calculating the Dark Flight phase.
  
11-In the Dark Flight folder are saved:
    a) A "Simple_Dark_flight" .txt file with the results of the simplified dark flight calculation.
    b) A .bmp figure showing the trajectory of the dark flight seen in profile starting from the shutdown point of the bolide.
  
    NOTE: for the calculation of the complete dark flight, i.e. that also takes into account the winds and predicts the uncertainty of the strewn field,
    a more complete software not currently available in MuFiS is used. The difference between the two versions of dark flight (simplified and full)
    can reach several hundred meters.
  
12-In the "Orbit" folder are saved:
    a) A .css file with the heliocentric orbits of all the Monte Carlo clones that can be viewed with the Celestia software, a 3D simulation planetarium.
    b) A .css file with the nominal heliocentric orbit of the progenitor meteoroid of the bolide viewable with Celestia.
    c) A .txt file with the nominal orbital elements associated with the uncertainties that come out of the Monte Carlo computation.
    d) An image with the plot of the nominal orbit of the progenitor meteoroid of the bolide projected on the plane of the Ecliptic, the orbits of the planets and their position at the moment of the bolide.
    e) An image with the plot of the nominal orbit plus those of the Monte Carlo clones of the progenitor meteoroid of the bolide projected on the plane of the Ecliptic, the orbits of the planets and their position at the time of the bolide.
   
13-In the case of basic analysis, the processing stops at the triangulation, at the heliocentric orbit and - if necessary - at the dark flight, the "Dynamic" folder is not created.

INPUT FORMATS

MuFiS accepts two types of input formats. The standard one (which is activated with the letter S as input) is of the type:

 % Simulated Fireball:IT20000000
 %    
 % TS (degree): 217.6840167 
 %    
 % Time (s)                  RA (??)              Dec (??)             Azimut (??)          Altezza (??) 
 2457904.381446667500        220.307336897       0.176752764         176.286510199       45.035025925 
 2457904.381447509400        220.402238239       0.316915535         176.142948636       45.170612581 
 2457904.381448286100        220.498727393       0.459392132         175.996288993       45.308252280 
 2457904.381448975800        220.597141538       0.604677262         175.845976439       45.448411305 
 2457904.381449813000        220.697386299       0.752628064         175.692108906       45.590939125
 
 That of PRISMA (which is activated with the letter P in input) is instead of the type (the actual data begin after line n. 10):
 
 file_solution = ITTO01_201804_astro_solution.txt
file_error    = ITTO01_201804_astro_error.txt
file_covar    = ITTO01_201804_astro_covar.txt
C             = 8.408
k             = -0.317
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                          date          julian_date           az         s_az          alt        s_alt           ra         s_ra          dec        s_dec          mag        s_mag
                           [/]                  [/]        [deg]        [deg]        [deg]        [deg]        [deg]        [deg]        [deg]        [deg]          [/]          [/]
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
    2019-05-07T19:42:02.649880     2458611.32086400     189.3510       0.0041      11.4626       0.0180     160.3859       0.0241     -33.7917       0.0204       -4.124        0.129
    2019-05-07T19:42:02.716650     2458611.32086478     189.6228       0.0036      11.3339       0.0206     160.0465       0.0277     -33.8819       0.0234       -4.183        0.124
    2019-05-07T19:42:02.750050     2458611.32086516     189.7445       0.0035      11.3832       0.0194     159.9131       0.0260     -33.8156       0.0220       -3.823        0.122
    2019-05-07T19:42:02.816800     2458611.32086594     190.0219       0.0037      11.2405       0.0197     159.5638       0.0266     -33.9173       0.0225       -3.981        0.103
    2019-05-07T19:42:02.850200     2458611.32086632     190.1552       0.0038      11.2071       0.0195     159.4019       0.0264     -33.9308       0.0223       -4.040        0.119
    2019-05-07T19:42:02.916950     2458611.32086709     190.4706       0.0035      11.0578       0.0195     159.0055       0.0265     -34.0313       0.0223       -4.255        0.097
    2019-05-07T19:42:02.950350     2458611.32086748     190.5754       0.0034      11.1305       0.0188     158.8971       0.0256     -33.9428       0.0217       -4.348        0.086
    2019-05-07T19:42:03.017130     2458611.32086825     190.8710       0.0032      10.9410       0.0190     158.5148       0.0261     -34.0844       0.0220       -4.940        0.055

These data refer to a single observing station and must be contained in a text file like "ITYYYYYYMMDD_n.txt" where YYYY=year, MM=month, DD=day and n=1, 2, 3... is the number progressive that identifies the station. For both input formats it is necessary to prepare a separate file with the geographic coordinates of the stations, "Stazioni_ITYYYYMMDD.txt", with a structure like:

 % Simulated Fireball: IT20000000
 %    
 % Lat. (??) 		 Long. (??) 			 Height (m)				 Progressive number of the station
 45.081700 		  11.795000 		      15.00 		          1 
 44.817600 		  12.010600 		       9.00 		          2 
 44.917600 		  11.910600 		      29.00 		          3 
 45.176000 		  11.810600 		      50.00 		          4 
 
*******************
* MATLAB Compiler *
*******************

1. Prerequisites for Deployment 

. Verify the MATLAB Runtime is installed and ensure you    
  have installed version 9.0 (R2015b).   

. If the MATLAB Runtime is not installed, do the following:
  (1) enter
  
      >>mcrinstaller
      
      at MATLAB prompt. The MCRINSTALLER command displays the 
      location of the MATLAB Runtime installer.

  (2) run the MATLAB Runtime installer.

Or download the Windows 64-bit version of the MATLAB Runtime for R2015b 
from the MathWorks Web site by navigating to

   http://www.mathworks.com/products/compiler/mcr/index.html
   
   
For more information about the MATLAB Runtime and the MATLAB Runtime installer, see 
Package and Distribute in the MATLAB Compiler documentation  
in the MathWorks Documentation Center.    


NOTE: You will need administrator rights to run MCRInstaller. 


2. Files to Deploy and Package

Files to package for Standalone 
================================
-Pipeline_fireball_comp.exe
-MCRInstaller.exe 
   -if end users are unable to download the MATLAB Runtime using the above  
    link, include it when building your component by clicking 
    the "Runtime downloaded from web" link in the Deployment Tool
-This readme file 

3. Definitions

For information on deployment terminology, go to 
http://www.mathworks.com/help. Select MATLAB Compiler >   
Getting Started > About Application Deployment > 
Deployment Product Terms in the MathWorks Documentation 
Center.
