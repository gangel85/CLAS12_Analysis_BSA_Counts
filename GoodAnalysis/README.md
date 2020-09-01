**TO RUN THE CODE: **

Copy it on the ifarm.
Go on Ifarm  where the MainFunction is and to compile it run javac:

 javac -classpath ".:/w/hallb-scifs17exp/clas12/gangel/jar/jaw-2.0.jar:/w/hallb-scifs17exp/clas12/gangel/jar/clas-jcsg-6.2.0-SNAPSHOT.jar:/w/hallb-scifs17exp/clas12/gangel/jar/vecmath-1.3.1-2.jar" *.java

To run the compiled code:

nohup java -classpath "/lustre19/expphy/volatile/clas12/gangel/Analysis/:/w/hallb-scifs17exp/clas12/gangel/jar/jaw-2.0.jar:/w/hallb-scifs17exp/clas12/gangel/jar/clas-jcsg-6.2.0-SNAPSHOT.jar:/w/hallb-scifs17exp/clas12/gangel/jar/vecmath-1.3.1-2.jar" GoodAnalysis.MainFunction 0 0 /lustre19/expphy/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/ /lustre19/expphy/volatile/clas12/gangel/BSA_Analysis/ID_5/ 5


The first adress after classpath is the adress where the main folder minus a subfolder. Basically $PWD with a cd.. 
In this case the main is in : 
/lustre19/expphy/volatile/clas12/gangel/Analysis/GoodAnalyis/
So we need to provide the path up to a folder earlier: /lustre19/expphy/volatile/clas12/gangel/Analysis/
after the column (and separated by a new column) we have a path to all the external jar libraries used in my analysis code. 

After that we need to provide the name of the folder of the Mainfunction: GoodAnalysis followed by a dot and the file containing the main:
GoodAnalysis.MainFunction 

Be sure that in MainFunction the package is set to match the name of the folder "GoodAnalysis" otherwise java will complain and saying that Main not found.

The main functions argument are 5 and the following :

args[0] is integer with values  0 or 1 . 0 if you run the code on data, 1 if you run it on MC

args[1] is integer with values 0 or 1 . If 0 it runs on all the statistics, if 1 it run on a sub sample of 2M events used to debug.

args[2] ia a string, containing the path to the file to be used 

args[3] is a string containing the ouput directory of the analysis (be sure to to make it before running the code otherwise it will not save anything)

args[4] is an integer 0,1,2,3,4,5 . With the following menaings : 0 run the code with 1D binning in Q2. 1 run the code with 1D binning in xB. 2 run the code with 1D binning in z, 3 run the code with 1D binning in Pt, 4 run the code with 9 binning in XB and Q2 integrating z, 5 run the cude multidimensionally (9 bins in xB and Q2 and 7 bins in z and 7 bins in PT)


**OUTPUT OF THE CODE:**

The code will return a folder called Multibins divided in as many folders as the Q2  bins with subfloder equal to the xB bins.
Inside each folder a summary of the number of pions find for Z and PT bins (txt) as well as home hipo files.
These hipo files contains histograms of counts for the number of electrons, the helicity counts vs phi, the total counts of pions for each Z and PT bin, as well as the kinematic distributions inside the Z and PT bin. 
They are the following.

If you want to read them us the code:
BSA_DoubleCheck_2020.java

modify the file such to have the right String with the path of your folder , and the index corresponding to the index you used to run the code : 0,1,2,3,4,5.





**===== Analysis Code explanation ====**

The code define MultiBins object, these objects are double arrays of H1F, H2F or H3F, allowing to have a 5D binning. 
The first two variables are stored in array sand their allocation can be done using constand bin size, as well as defiying custommade polynomial bins shapes using the subrutine: PolyogonalBins in the MultiBins.java file . They have been set to match the Stefan's Diehl binning for BSA analysis but can be modified .
When you define a bin you need to privde the points in counter clock wise manner starting from the top right corner of the polygon.  
The class map will map that bin into an index i, j for the array. be sure that you modify that too.

The Class AnalysisSIDIS.java deal with filling those bins for all the relevant infromation and later save on an hipo files.

In order to enhance the Particle class,a ParticleREC.java class has been created. This class will add to the particle the Traj and Calorimeter banks so to be passed to the Fiducial cuts routine.

A class Photons is used to define good photons based on his internal cuts. 

THe charged hadrons are built using two classes: ChargedParticle, that read a ParticleREC and define the properites of that hadron by computing Phi Trento ,Z , PT, pseudorapidity (not working yet) 

Pi0 are created using: Pi0_Particle that take the Photons class and performe the combination of photons to provide the pi0. Several fitting procedures are available but only one has been implemented in the current status.



**Fiducial Cuts:**
The fiducial cuts are created trough the usage of a MAP. The MAP is created with a key tha tis the PID and a list of cuts that are going to be activated for that PID. For exampole now for PID 11 there are 3 cuts activated : Sampling Fraction, PCAL, DC. While for pid 211 only DC.  This is the standard implementation of the class but an user ca modify this map of cuts as wish , adding and removing particles or detectors cuts or adding a new .java file containing a new cut. The new cut should be written as " impements DetectorCut" so to inherit all the right properties.

The function Status loop over the  MAP and return the status of each detector cut  by invoking each detector status class.  If all of them are true it will return true for that PID, meaning that particle has passed all associated cuts. If at least one is false will return false and therefore the particle do not pass fiducial cuts. In this way every detector (DC.java, EC.java, SamplingFraction.java )  will apply their routiine to the PID and return the status : true false that is being read by the main Fiducial Class.






