# CLAS12_Analysis_BSA_Counts


TO RUN THE CODE: 
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

OUTPUT OF THE CODE:

The code will return a folder called Multibins divided in as many folders as the Q2  bins with subfloder equal to the xB bins.
Inside each folder a summary of the number of pions find for Z and PT bins (txt) as well as home hipo files.
These hipo files contains histograms of counts for the number of electrons, the helicity counts vs phi, the total counts of pions for each Z and PT bin, as well as the kinematic distributions inside the Z and PT bin. 
They are the following.

If you want to read them do the following: 

An example for the extraction of the BSA code appleid on this output is in the subfloder: BSA Analysis Code. 
