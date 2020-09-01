/**
 * 
 */
package GoodAnalysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.H3F;
//import org.jlab.io.base.DataEvent;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.physics.EventFilter;
import org.jlab.jnp.physics.LorentzVector;
import org.jlab.jnp.physics.Particle;
import org.jlab.jnp.physics.PhysicsEvent;
import org.jlab.jnp.reader.DataManager;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo.data.HipoEvent;
import org.jlab.jnp.hipo.data.HipoGroup;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.data.Schema;
/**
 * @author gangelini
 * Function to compute Pion BSA and Multiplicity raw counts (Counts in Phi, counts integrated in phi)
 */
public class MainFunction {

	// MultiBins is an object that allow to define custom binning up to 6Dimensions
	private static MultiBins InvariantMassBins,PionCounts,Counts_Phi,Helicity0,Helicity1,MissingM,MCParticles,Electron_counts ;

	// Analysis
	//0/1 (MC no, MC yes),0/1 (debug no, debug yes), 0/1  (Multi dimensional no=0 , yes =1 ), reading folder , writing folder 
	public static void main(String[] args) throws FileNotFoundException {



		// Select the particle you want to analyze
		int Contaunelettrone =0 ;
		int Contatuttielettroni=0;
		//	int ContaTuttiIpioni=0;
		int ContaiPioniBuoni=0;
		//	int contapionebanca=0;



		boolean MC=false;
		boolean debug=true;
		int goodeventsStop=2000000;
		String workdir = "/Volumes/My Passport/Stefan/";

		//String workdir = "/Volumes/My Passport/Problem/";

		String workdirout ="/Users/gangelini/work/Analysis_Results/Test_0/" ;	
		int valueBinStatus=2;

		/*

	 System.out.println(" Running the Analysis code by Giovanni Angelini  ");
	 System.out.println("The arguments need to be 0/1 (MC no, MC yes),0/1 (debug no, debug yes), 0/1  (Multi dimensional no=0 , yes =1 ), reading folder , writing folder ");
     boolean MC=false;
	 if (Integer.parseInt(args[0])==1) MC=true;
	 boolean debug=false;
     if (Integer.parseInt(args[1])==1) debug=true;
	 int goodeventsStop=50000000;
	 String workdir = args[2];
	 String workdirout = args[3];
	 int valueBinStatus=Integer.parseInt(args[4]);
     System.out.println(" ===================================");

     System.out.println("VALUE 4th variable input " + valueBinStatus);

		 */
		boolean PolygonalBinClas=true; // Define if for x and Q2 we arew using the polygon class 
		boolean SingleHadron= false;

		// CHOSE YOUR HADRON FOR ANALYIS 
		boolean pi0 = false;
		boolean piP = true;
		boolean piM = false;
		//int PionID = +321; // Analyzie Kaons
		int PionID = +211;	
		if (piM==true) PionID=-211;	

		//CUTS LIMITS
		double W2Cut = 2 ;
		double Q2Cut = 1;
		double YCut = 0.8;
		double E_Mom_Min = 2.1;
		double El_min_Angle =5;
		double Vertex_Max = 12;
		double Vertex_Min = -13;
		//Hadron:
		double Min_H_momentum = 1.25;
		double Max_H_momentum = 5;
		double Min_H_theta=5;
		double Max_H_theta=35;
		double Cut_MissingMass_H=1.5;



		boolean Q21D=false; boolean xB1D=false; boolean z1D=false; boolean PT1D=false;
		boolean Qx2D=false; boolean QxzPt=false;


		if (valueBinStatus==0) {Q21D=true;PolygonalBinClas=false;}
		if (valueBinStatus==1) {xB1D=true;PolygonalBinClas=false;}
		if (valueBinStatus==2) {z1D=true;PolygonalBinClas=false;}
		if (valueBinStatus==3) {PT1D=true;PolygonalBinClas=false;}
		if (valueBinStatus==4) {Qx2D=true;PolygonalBinClas=true;}
		if (valueBinStatus==5) {QxzPt=true;PolygonalBinClas=true;}


		int Q2_Bins=1; double Q_min= 1; double Q_max= 12; // 1 - 12 		  
		int Xb_Bins=1; double Xb_min= 0.01 ; double Xb_max= 0.9 ;  // 0 -0.9
		int Z_Bins= 7;  double Z_min= 0.3;   double Z_max= 1;  // 0.2 - 0.9
		int PT_Bins=1; double PT_min= 0;  double PT_max= 2;// 0.1 - 4
		int Mass_Bins= 100; double min_pion= 0.05;  double max_pion=0.25;	 
		int Phi_Bins= 12; double Phi_min=0.0; double Phi_max=360.0;

		//Common Bins Phi and Missing mass
		int bincount_phi = 12;
		double bin_phi[] = {0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 360.0};
		// To be incremented to 100 by using Axis class and self division:
		int binmass=20;
		double bin_mass[]= {0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25};

		if (valueBinStatus==0) {
			int bincount_q2 = 13;
			double bin_q2[] = {1.3, 1.65, 1.85, 2.05, 2.3, 2.6, 3.0, 3.5, 4.1, 4.9, 5.8, 6.8, 8.0, 11.0};
			int bincount_x = 1;
			double bin_x[] = {Xb_min,Xb_max};
			int bincount_z = 1;
			double bin_z[] = {Z_min,Z_max};
			int bincount_pt = 1;
			double bin_pt[] = {PT_min,PT_max};

			InitializeMultiBins( PolygonalBinClas,  bincount_q2, bincount_x,  bincount_z,  bincount_pt,  binmass,  bincount_phi,
					bin_q2,bin_x,  bin_z,  bin_pt, bin_mass,  bin_phi);
		}

		if (valueBinStatus==1) {
			int bincount_q2 = 1;
			double bin_q2[] = {Q_min,Q_max};
			int bincount_x = 11;
			double bin_x[] = {0.09, 0.145, 0.18, 0.215, 0.25, 0.29, 0.33, 0.37, 0.42, 0.48, 0.55, 0.70};
			int bincount_z = 1;
			double bin_z[] = {Z_min,Z_max};
			int bincount_pt = 1;
			double bin_pt[] = {PT_min,PT_max};

			InitializeMultiBins( PolygonalBinClas,  bincount_q2, bincount_x,  bincount_z,  bincount_pt,  binmass,  bincount_phi,
					bin_q2,bin_x,  bin_z,  bin_pt, bin_mass,  bin_phi);
		}

		if (valueBinStatus==2) {
			int bincount_q2 = 1;
			double bin_q2[] = {Q_min,Q_max};
			int bincount_x = 1;
			double bin_x[] = {Xb_min,Xb_max};
			int bincount_z = 12;
			double bin_z[] = {0.15, 0.21, 0.25, 0.29, 0.34, 0.40, 0.46, 0.52, 0.57, 0.63, 0.7, 0.79, 0.95};
			int bincount_pt = 1;
			double bin_pt[] = {PT_min,PT_max};
			InitializeMultiBins( PolygonalBinClas,  bincount_q2, bincount_x,  bincount_z,  bincount_pt,  binmass,  bincount_phi,
					bin_q2,bin_x,  bin_z,  bin_pt, bin_mass,  bin_phi);
		}
		if (valueBinStatus==3) {
			int bincount_q2 = 1;
			double bin_q2[] = {Q_min,Q_max};
			int bincount_x = 1;
			double bin_x[] = {Xb_min,Xb_max};
			int bincount_z = 1;
			double bin_z[] = {Z_min,Z_max};
			int bincount_pt = 14;
			double bin_pt[] = {0.0, 0.12, 0.22, 0.295, 0.36, 0.425, 0.5, 0.595, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.7};
			
			InitializeMultiBins( PolygonalBinClas,  bincount_q2, bincount_x,  bincount_z,  bincount_pt,  binmass,  bincount_phi,
					bin_q2,bin_x,  bin_z,  bin_pt, bin_mass,  bin_phi);
		}			


		if(valueBinStatus==4) {
			int bincount_q2 = 4;
			double bin_q2[] = {Q_min,Q_max};
			int bincount_x = 7;
			double bin_x[] = {Xb_min,Xb_max};
			int bincount_z = 1; 
			double bin_z[] = {Z_min,Z_max};
			int bincount_pt = 1; 
			double bin_pt[] = {PT_min,PT_max};

			InitializeMultiBins( PolygonalBinClas,  bincount_q2, bincount_x,  bincount_z,  bincount_pt,  binmass,  bincount_phi,
					bin_q2,bin_x,  bin_z,  bin_pt, bin_mass,  bin_phi);

		}

		if(valueBinStatus==5) {
			int bincount_q2 = 4;
			double bin_q2[] = {Q_min,Q_max};
			int bincount_x = 7;
			double bin_x[] = {Xb_min,Xb_max};
			int bincount_z = 7; 
			double bin_z[] = {0.15, 0.2, 0.24, 0.29, 0.36, 0.445, 0.55, 0.7};
			int bincount_pt =7; 
			double bin_pt[] = {0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0};

			InitializeMultiBins( PolygonalBinClas,  bincount_q2, bincount_x,  bincount_z,  bincount_pt,  binmass,  bincount_phi,
					bin_q2,bin_x,  bin_z,  bin_pt, bin_mass,  bin_phi);

		}

		//This is PT2 
		if(PolygonalBinClas==true) {

			// We will create extra bins , THe overflow and udnerfloow will be in Q2 =4 and xb 7.
			//The only utalize are 2 bins of Q2 and 4 of xB 
			
			Q2_Bins=4; Q_min= 1;     Q_max= 12; // 1 - 12 		  
			Xb_Bins=7; Xb_min= 0.01; Xb_max= 0.9 ;  // 0 -0.9
			Z_Bins= 1; Z_min= 0.3;   Z_max= 0.7;  // 0.2 - 0.9
			PT_Bins=1;  PT_min= 0;   PT_max= 2;// 0.1 - 4

		}


		// Used to look into invariant mass distribution
		H1F XMass = new H1F("XMass",200,0,2);

 
		Histos Histo = new Histos(workdirout); // To save on file electron and photons distributions
		File folder = new File(workdir); File[] listOfFiles = folder.listFiles();
		long count =0 ; // Total Event counters 
		long good_e_count=0; // Selected electrons
		long all_electrons=0; // Total Electrons reconstructed in events 
		long twophotons=0; //Events with 2 photons candidates to pi0 

		FiducialCuts SIDIS_Cuts = new FiducialCuts();  // This creates a fiducial cut with standardized cuts for my analysis
		SIDIS_Cuts.GetMap().show(); // print the map of cutsssociate with each PID. It removes Data in region of detector that cannot be trusted. 
		
		int helicity=-9; //  Helicity inizialization  
		double beamEnergy = 10.603; // Set the energy of the beam
		LorentzVector beam   = new LorentzVector(0.0,0.0,beamEnergy,beamEnergy); // Getting the Lorentz Vector associated to the beam
		LorentzVector target = new LorentzVector(0.0,0.0,0.0,0.93827); // Defining Lorentz Vector for the proton target
		
		//Filtering the events properly
		EventFilter filter   = new EventFilter();
		if (pi0==true ) filter.setFilter("11:22:22:X+:X-:Xn"); // Search in the events for electron , 2 photons , and any other particle 
		else if (piM==true ) filter.setFilter("11:-211:X+:X-:Xn");
		else if (piP==true ) filter.setFilter("11:211:X+:X-:Xn");

		System.out.println(" I am filtering the events using the scheme: ");System.out.println(filter);

		// Inizializing the Particles and Analysis tools 
		Pi0_Particle New_Pi0s = new Pi0_Particle(); // Pi0Particle 
		// Routine for the SIDIS analysis:
		AnalysisSIDIS Analysis = new AnalysisSIDIS(PolygonalBinClas,InvariantMassBins,Z_Bins, Z_min, Z_max, PT_Bins, PT_min, PT_max, SingleHadron); // Create an Analysis of Pi0
		
		// Analysis on the MC generated particles: 
		MCAnalysis LUND_pi = new MCAnalysis(InvariantMassBins, Z_Bins, Z_min, Z_max, PT_Bins, PT_min, PT_max, Phi_Bins); // Create an Analysis of MonteCarlo files
		LUND_pi.setBeam(beam);LUND_pi.setTarget(target); //Used to compute the PT, Phi, etc. from lund momenta


		int evento=0;
		// Reading the files 

	//	File File_txt = new File("/Users/gangelini/work/Analysis_Results/results.txt");
	//	File File_txt2 = new File("/Users/gangelini/work/Analysis_Results/el_ID_gio.txt");
	//	File File_ele = new File("/Users/gangelini/work/Analysis_Results/el_values_gio.txt");
	//	File File_pip = new File("/Users/gangelini/work/Analysis_Results/pip_values_gio.txt");

			for (File file : listOfFiles) {
				if (file.isFile()&& file.getName()!="skim4_005234.hipo") { // FILE 005234 make the code crashing one the banks are read
					System.out.println(" Reading  " + file.getName());
					HipoReader reader = new HipoReader(); reader.open(workdir+file.getName());
					Event event=new Event(); 
					SchemaFactory schema = reader.getSchemaFactory();
					Bank RRun = new Bank(schema.getSchema("RUN::config"));
					Bank REvent = new Bank(schema.getSchema("REC::Event"));
					Bank RCal = new Bank(schema.getSchema("REC::Calorimeter"));
					Bank RPart = new Bank(schema.getSchema("REC::Particle"));
					Bank RTraj = new Bank(schema.getSchema("REC::Traj"));
					Bank RCherenkov = new Bank(schema.getSchema("REC::Cherenkov"));
					
					double LorentzProduct=0, xB=0, Pq=0, Pl=0, y_var=0,epsilon=0,epsilon2= 0,gamfactor=0,epsilonnum=0,epsilonden=0;
					double e_mom=0, e_Theta =0,startime=0;
					 
					while(reader.hasNext()==true){
						count++;
						boolean checkEle = false;
						reader.nextEvent(event);
						event.read(REvent); event.read(RCal); event.read(RPart);event.read(RTraj); event.read(RRun);event.read(RCherenkov);
						evento = RRun.getInt("event", 0);
						if(count%100000 == 0) System.out.println(count + " events");
						if(good_e_count%50000 == 0) System.out.println(good_e_count + " good electrons");	
						
						helicity= REvent.getByte("helicity",0);
						startime=(double)REvent.getFloat("startTime", 0);
						double p_pi =0;
						
						if(RPart.getRows()>0 && RCal.getRows()>0 && RTraj.getRows()>0) // Be sure I have all the banks I need
						{
							PhysicsEvent phsEvent = DataManager.getPhysicsEvent(beamEnergy,RPart) ;
							int TriggerParticleID = RPart.getInt("pid",0); //Get the ID of the trigger
						
							//Check on the HTTC signal: 
							boolean HTTCSignal=false;
							double nphe =0;
							for(int i=0; i<RCherenkov.getRows(); i++) {
								if( RCherenkov.getInt("pindex", i)==0 && RCherenkov.getInt("detector", i)==15 && RCherenkov.getFloat("nphe", i) >2 ) {
									nphe=RCherenkov.getFloat("nphe", i);
									HTTCSignal=true;
								}
							}
							double min_mom=0, Q2=0, multifact=0;
							int indicepione=-9;
							//==> First selection : Filter, 1 electron only; Trigger is electron , HTTC NPE>2 
							if(  filter.isValid(phsEvent)==true && phsEvent.countByPid(11)==1&&TriggerParticleID==11 && HTTCSignal==true) {
								
								
								//Read the electron:
								Particle electron = phsEvent.getParticleByPid(11, 0);
								// Extend the electron by defining a Reconstructed particle, menaing attaching to it the Part,Trajectory and Calorimeter hits.
								ParticleREC electronRec = new ParticleREC(electron,RPart,RTraj,RCal);
								all_electrons++;// I found an electron 		
								LorentzVector vecW2 = new LorentzVector(0,0,0,0);
								LorentzVector vecQ2 = new LorentzVector(0,0,0,0);
								LorentzVector vecE = LorentzVector.from(electron.vector()); // electron.vector() returns the LorentzVector of the particle electron
								vecW2.add(target);vecW2.add(beam);vecW2.sub(vecE);vecQ2.add(beam); vecQ2.sub(vecE);
								Pq = -target.px()*vecQ2.px()-target.py()*vecQ2.py()-target.pz()*vecQ2.pz()+target.e()*vecQ2.e();
								Pl = -target.px()*beam.px()-target.py()*beam.py()-target.pz()*beam.pz()+target.e()*beam.e();
								xB = -vecQ2.mass2()/(2*Pq); 
								y_var = Pq/Pl;
								Q2= -vecQ2.mass2();
								epsilon2 = ( 1 - y_var - Q2 / (4*10.603*10.603) ) / ( 1 - y_var + 0.5*y_var*y_var + Q2 / (4*10.603*10.603) );
								gamfactor=(2*0.93827*xB)/Math.sqrt(-vecQ2.mass2());
								multifact= (0.25)*Math.pow(y_var, 2)* Math.pow(gamfactor, 2);
								epsilonnum=1-y_var- multifact ;
								epsilonden=1-y_var+(0.5)*Math.pow(y_var, 2)+multifact;
								epsilon=epsilonnum/epsilonden;
								
								
					
								e_mom=electron.p();
								e_Theta = Math.toDegrees(electron.theta());
								
								//===> 2 Cut: Kinemaitcs and electron angle and momentum
								if(vecW2.mass()>= W2Cut && -vecQ2.mass2()>= Q2Cut && y_var<= YCut && e_mom>= E_Mom_Min && electron.vz()> Vertex_Min && electron.vz()< Vertex_Max && Math.toDegrees(electron.theta())>El_min_Angle)
								{    	    
									Contatuttielettroni++;
									if(phsEvent.countByPid(11)>1) Contaunelettrone++;	
									
									//Checking that the electron passes the fiducial cuts:
									if(SIDIS_Cuts.Status(electronRec)==true) {	

										Electron_counts.FillH2(-vecQ2.mass2(),xB);// Filling the electron counts 
										good_e_count++; //Good electron SISDIS within all cuts	

							
										//Check the number of pions in the event
										int npioni = phsEvent.countByPid(PionID);
										boolean GoodPi=false;
										for (int kk=0; kk<npioni ; kk++) {
											Particle pion_positive2 = phsEvent.getParticleByPid(PionID, kk);
											//Defined the Particle REC with the info from the other banks attached
											ParticleREC Pione = new ParticleREC(pion_positive2,RPart,RTraj,RCal);
											Pione.changePid(PionID); //Be sure that the PID is set to 211
											
											// Looping over the REC::PARTICLE bank to identify the pion I am looking at:
											// For the purpose of Chi2PID 
											int indice=0;
											for(int i =0 ; i< RPart.getRows(); i++ ) {
												if((PionID)==RPart.getInt("pid", i))  {
													p_pi = Math.sqrt(RPart.getFloat("px", i)*RPart.getFloat("px", i)+RPart.getFloat("py", i)*RPart.getFloat("py", i)+RPart.getFloat("pz", i)*RPart.getFloat("pz", i));
													if ( pion_positive2.p() <= p_pi+0.0001 && pion_positive2.p() >= p_pi-0.0001 ) {
														indice=i; //Index of my pions in the REC::Particle 
													}											 
												}
											}

											double Chi2Pid2= RPart.getFloat("chi2pid",indice);
											double chi2pid_max=0;
											double chi2_pid_limit_lower= -3.0*0.88;
											if(Pione.p() < 2.44){ chi2pid_max = 3;}
											else{ chi2pid_max = (0.00869 + 14.98587 * Math.exp(-Pione.p() / 1.18236) + 1.81751 * Math.exp(-Pione.p() / 4.86394));}
											chi2pid_max=0.88*chi2pid_max;
											
											GoodPi=SIDIS_Cuts.Status(Pione); //Boolean for Fiducial cuts 
											
											// ====> Selection of pion: Vertex within 20 cm from electron , Chi2PID in the limits , FIDUCIAL CUTS Passed
											if(Math.abs(Pione.vz()-electron.vz())<20 &&  Chi2Pid2 > chi2_pid_limit_lower && Chi2Pid2 < chi2pid_max && GoodPi ) {	
												
												// For this good Pion generated trhe class ChargedParitcle, meaning get Z,PT,PHI, Pseudorapidity etc. 
												ChargedPi_Particle New_Pi= new ChargedPi_Particle(); // Charged Pi class
												New_Pi.getChargedPi(pion_positive2, beam,target,vecE,vecQ2,vecW2,xB,-vecQ2.mass2(),0.0);

												//Be sure the Z cut applies only at the integrated Z cases 
												boolean ZCUT=false;		    		 
												if(valueBinStatus==2 || valueBinStatus==5)ZCUT=true;
												else if(valueBinStatus!=2 && valueBinStatus !=5 && New_Pi.getZ()>Z_min)ZCUT=true;
												
											//==>Selection: XCUT, Pi Momentum, angle and missing mass 
												if(ZCUT&&  New_Pi.getMomentum()>Min_H_momentum && New_Pi.getMomentum()<Max_H_momentum && New_Pi.getTheta()>Min_H_theta &&New_Pi.getTheta()<Max_H_theta &&New_Pi.getMissMass()>Cut_MissingMass_H) {			
													ContaiPioniBuoni++; // Count the good pions
													//Run the analysis on this Pion (filling counts)
													Analysis.AnalyzeCharged( PolygonalBinClas,beam,target,electronRec,helicity,XMass,MissingM,3,xB,-vecQ2.mass2(), New_Pi, PionCounts,Counts_Phi, Helicity0, Helicity1,Electron_counts,vecW2.mass(),y_var,epsilon);

												} // Loop over invariant mass and other cuts
											}//  Chi 2 and Status Check 
										}// Loop over all Hadrons 
									} // Good Electric status Fiducial 
								} // DIS electron
							} // Filter 

						}// If I have the rows
						if (debug==true) {
							//if(good_e_count==200000 ) {
							if(good_e_count==goodeventsStop ) {
								//if(good_e_count==8000000 ) {
								// my photons are
								System.out.println("Two photons events with photon cuts " + twophotons);   
								System.out.println("All count"+ count+ " while " + good_e_count + " are FINAL good electrons");	
								System.out.println(all_electrons + " FINAL all electron filtered ");
								SimpleDateFormat dateFormat = new SimpleDateFormat("MMMdd_hh_mm");
								Date now = new Date();
								String time = dateFormat.format(now);
								String time2; 
								if (MC==true) {
									time2 = now+"MC";
									LUND_pi.SaveHistograms_MultiBins(time2,workdirout); // produce the MC histograms
									Histo.Print_ChargedHadron(time2); // this print a comparison Generated- Reconstructed Hadron in MC 
								}
								else  {time2 = time;}
								System.out.println("Printing histos ");
								Histo.PrintUnSkim_el(time2); Histo.PrintSkim_el(time2);Histo.PrintUnSkim_ph(time2); Histo.PrintSkim_ph(time2);
								SIDIS_Cuts.Print(time2,workdirout);
									       
							
								System.out.println(" --- I am analyzing Charged Pions  " );
								Analysis.CountPions(XMass, MissingM, PionCounts, Counts_Phi,Helicity0,Helicity1, Electron_counts,time2, MC,workdirout);
								System.out.println(" ============================================= ");
								System.out.println(" Events read " + count);
								System.out.println(" Elettroni DIS " +good_e_count);//printing the invariant mass for pi0 
								System.out.println(" All Electrons " +all_electrons);//printing the invariant mass for pi0
								System.out.println(" All Good Pions" +ContaiPioniBuoni);//printing the invariant mass for pi0
								System.out.println(" ============================================= ");
								System.out.println(" ------ BININNG USED IN Z AND PT ARE ----");
								System.out.println(" Z " + Z_Bins);
								System.out.println(" Pt " + PT_Bins);
								return;
							} // if number of events = X
						}   // if debug is true  
					}// Reader has next file			
				} // A good File 
			}// List of Files
	

		System.out.println("Two photons events with photon cuts " + twophotons);   
		System.out.println("All count"+ count+ " while " + good_e_count + " are FINAL good electrons");	
		System.out.println(all_electrons + " FINAL all electron filtered ");
		SimpleDateFormat dateFormat = new SimpleDateFormat("MMMdd_hh_mm");
		Date now = new Date();
		String time = dateFormat.format(now);
		String time2; 
		if (MC==true) {
			time2 = now+"MC";
			LUND_pi.SaveHistograms_MultiBins(time2,workdirout); // produce the MC histograms
		}
		else  {time2 = time;}
		System.out.println("Printing histos ");		
		if (pi0==true) {
			System.out.println(" --- I am analyzing Neutral Pions  " );
			Analysis.CountPi0s(1, 3,XMass,MissingM, InvariantMassBins,Helicity0,Helicity1, Electron_counts,time2, MC,workdirout);
		}
		else {
			System.out.println(" --- I am analyzing Charged Pions  " );
			Analysis.CountPions(XMass, MissingM, PionCounts,Counts_Phi, Helicity0,Helicity1, Electron_counts,time2, MC,workdirout);
		}

		System.out.println(" ============================================= ");
		System.out.println(" Events read " + count);
		System.out.println(" Elettroni DIS " +good_e_count);//printing the invariant mass for pi0 
		System.out.println(" All Electrons " +all_electrons);//printing the invariant mass for pi0
		System.out.println(" Good Pions are" + ContaiPioniBuoni);

		System.out.println(" ============================================= ");
		LUND_pi.getNr("z");
		LUND_pi.Histo(time2,workdirout);
		System.out.println(" ------ BININNG USED IN Z AND PT ARE ----");
		System.out.println(" Z " + Z_Bins);
		System.out.println(" Pt " + PT_Bins);



	}
	
/**
	 * @param fromBank the bank containing the index variable
	 * @param idxVarName the name of the index variable
	 * @return map with keys being the index in toBank and values the indices in fromBank
	 */
	private static  Map<Integer,List<Integer>> loadMapByIndex(	Bank rec_Calorimeter,String idxVarName) {
		Map< Integer,List<Integer> > map = new HashMap <Integer, List<Integer> >();
		if (rec_Calorimeter!=null) {
			for (int iFrom=0; iFrom<rec_Calorimeter.getRows(); iFrom++) {
				//System.out.println(" IFROM + " + iFrom);
				final int iTo = rec_Calorimeter.getInt(idxVarName,iFrom);
				if (!map.containsKey(iTo)) map.put(iTo,new ArrayList<Integer>()); 
				//	System.out.println("iTO e' " + iTo);
				map.get(iTo).add(iFrom);

			}
		}
		else {System.out.println(" Banca e' nulla? Non dovrebbe! ");
		};
		return map;
	}

	private static void InitializeMultiBins(boolean PolygonalBinClas, int bincount_q2,int bincount_x, int bincount_z, int bincount_pt, int binmass, int bincount_phi,
			double[] bin_q2,double[] bin_x, double[] bin_z, double[] bin_pt, double[] bin_mass, double[] bin_phi) {
		int Q2_Bins=bincount_q2;
		int Xb_Bins=bincount_x;
		int Z_Bins=bincount_z;
		int PT_Bins=bincount_pt;
		int Phi_Bins=bincount_phi;
		int Mass_Bins=binmass;
		InvariantMassBins = new MultiBins(Q2_Bins,Xb_Bins,Z_Bins,PT_Bins,Mass_Bins);
		InvariantMassBins.SetName_1stVariable("Q2");InvariantMassBins.SetName_2ndVariable("Xb");
		InvariantMassBins.SetName_3rdVariable("z");InvariantMassBins.SetName_4thVariable("Pt");
		InvariantMassBins.SetName_5thVariable("Invariant Mass");
		System.out.println( " Multidimensional Analysis code ");
		System.out.println(" Bin1: " + InvariantMassBins.GetName_1stVariable() + " Bin2: " + InvariantMassBins.GetName_2ndVariable());
		System.out.println(" Bin3: " + InvariantMassBins.GetName_3rdVariable() + " Bin4: " + InvariantMassBins.GetName_4thVariable());
		InvariantMassBins.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_mass);
		InvariantMassBins.GenerateHistograms("InvariantMass");
		if(PolygonalBinClas==true) {
			InvariantMassBins.InizializeClasPolygons();
		}
		PionCounts = new MultiBins(Q2_Bins,Xb_Bins,Z_Bins,PT_Bins,0);
		PionCounts.SetName_1stVariable("Q2");PionCounts.SetName_2ndVariable("Xb");
		PionCounts.SetName_3rdVariable("z");PionCounts.SetName_4thVariable("Pt");	 System.out.println( " Multidimensional Analysis code ");
		System.out.println(" Bin1: " + PionCounts.GetName_1stVariable() + " Bin2: " + PionCounts.GetName_2ndVariable());
		System.out.println(" Bin3: " + PionCounts.GetName_3rdVariable() + " Bin4: " + PionCounts.GetName_4thVariable());
		PionCounts.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_mass);
		PionCounts.GenerateHistograms("PionCounts");
		if(PolygonalBinClas==true) { PionCounts.InizializeClasPolygons();}
		// Counting particle in Phi Bins 
		Counts_Phi = new MultiBins(Q2_Bins,Xb_Bins,Z_Bins,PT_Bins,Phi_Bins);
		Counts_Phi.SetName_1stVariable("Q2");Counts_Phi.SetName_2ndVariable("Xb");
		Counts_Phi.SetName_3rdVariable("z");Counts_Phi.SetName_4thVariable("Pt"); Counts_Phi.SetName_5thVariable("Phi");
		Counts_Phi.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_phi);
		Counts_Phi.GenerateHistograms("Counts_Phi");
		if(PolygonalBinClas==true) { Counts_Phi.InizializeClasPolygons();}

		// Counting particle with Helicity0 as function of Phi bins. 
		Helicity0 = new MultiBins(Q2_Bins,Xb_Bins,Z_Bins,PT_Bins,Phi_Bins);
		Helicity1 = new MultiBins(Q2_Bins,Xb_Bins,Z_Bins,PT_Bins,Phi_Bins);
		Helicity0.SetName_1stVariable("Q2");Helicity0.SetName_2ndVariable("Xb");
		Helicity0.SetName_3rdVariable("z");Helicity0.SetName_4thVariable("Pt");
		Helicity1.SetName_1stVariable("Q2");Helicity1.SetName_2ndVariable("Xb");
		Helicity1.SetName_3rdVariable("z");Helicity1.SetName_4thVariable("Pt");
		Helicity0.SetName_5thVariable("Phi");Helicity1.SetName_5thVariable("Phi");
		Helicity0.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_phi);
		Helicity1.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_phi);
		Helicity0.GenerateHistograms("Helicity0"); Helicity1.GenerateHistograms("Helicity1");
		if(PolygonalBinClas==true) {Helicity0.InizializeClasPolygons(); Helicity1.InizializeClasPolygons();}

		// Counting particle with MissinM as function of Phi bins. 
		//it was 200 bins now it is just 3 because not used 
		MissingM = new MultiBins(Q2_Bins,Xb_Bins,Z_Bins,PT_Bins,3);
		MissingM.SetName_1stVariable("Q2");Helicity0.SetName_2ndVariable("Xb");
		MissingM.SetName_3rdVariable("z");Helicity0.SetName_4thVariable("Pt");
		MissingM.SetName_5thVariable("XMass");Helicity1.SetName_5thVariable("XMass");
		MissingM.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,  new double[] {0,1,2,3});
		MissingM.GenerateHistograms("MissingM");
		if(PolygonalBinClas==true) {MissingM.InizializeClasPolygons();}
		// Creating MultiDimensional Histograms for MC Analysis of Given Particle 

		MCParticles = new MultiBins(Q2_Bins,Xb_Bins,Z_Bins,PT_Bins,0);
		MCParticles.SetName_1stVariable("Q2");MCParticles.SetName_2ndVariable("Xb");
		MCParticles.SetName_3rdVariable("z");MCParticles.SetName_4thVariable("Pt");
		MCParticles.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_phi );
		MCParticles.GenerateHistograms("MCParticles");
		if(PolygonalBinClas==true) {MCParticles.InizializeClasPolygons();}
		// Creating an H2F for the information of the electron 

		Electron_counts = new MultiBins(Q2_Bins,Xb_Bins,0,0,0);
		Electron_counts.SetName_1stVariable("Q2");MCParticles.SetName_2ndVariable("Xb");
		Electron_counts.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_phi );
		Electron_counts.GenerateHistograms("Electron_counts");
		if(PolygonalBinClas==true) { Electron_counts.SetPolygonal(PolygonalBinClas); // it will compute properly the borders}
		}

	}
}
