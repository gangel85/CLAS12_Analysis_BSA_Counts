package GoodAnalysis;

import java.util.ArrayList;
import java.util.List;

import org.jlab.groot.data.H2F;
import org.jlab.jnp.physics.Vector3;

public class Sampling_Fractions implements DetectorCut{
	double min_cut[] = {0.16,0.16,0.16,0.16,0.16};
	double max_cut[] = {0.30,0.30,0.30,0.30,0.30};
	H2F SampliFraction_Overall_pre = new H2F("SampliFraction_Overall_pre",400,1,9, 400, 0,0.4);
	H2F SampliFraction_S1_pre = new H2F("SampliFraction_S1_pre",400,1,9, 400, 0,0.4);
	H2F SampliFraction_S2_pre = new H2F("SampliFraction_S2_pre",400,1,9, 400, 0,0.4);
	H2F SampliFraction_S3_pre = new H2F("SampliFraction_S3_pre",400,1,9, 400, 0,0.4);
	H2F SampliFraction_S4_pre = new H2F("SampliFraction_S4_pre",400,1,9, 400, 0,0.4);
	H2F SampliFraction_S5_pre = new H2F("SampliFraction_S5_pre",400,1,9, 400, 0,0.4);
	H2F SampliFraction_S6_pre = new H2F("SampliFraction_S6_pre",400,1,9, 400, 0,0.4);
	H2F PCAL_ECAL_Corr_pre = new H2F("PCAL_ECAL_Corr_pre",400,0,0.35, 400, 0,0.35);
	H2F SampliFraction_Overall_post = new H2F("SampliFraction_Overall_post",400,1,9, 400, 0,0.4);
	H2F SampliFraction_S1_post = new H2F("SampliFraction_S1_post",400,1,9, 400, 0,0.4);
	H2F SampliFraction_S2_post = new H2F("SampliFraction_S2_post",400,1,9, 400, 0,0.4);
	H2F SampliFraction_S3_post = new H2F("SampliFraction_S3_post",400,1,9, 400, 0,0.4);
	H2F SampliFraction_S4_post = new H2F("SampliFraction_S4_post",400,1,9, 400, 0,0.4);
	H2F SampliFraction_S5_post = new H2F("SampliFraction_S5_post",400,1,9, 400, 0,0.4);
	H2F SampliFraction_S6_post = new H2F("SampliFraction_S6_post",400,1,9, 400, 0,0.4);
	H2F PCAL_ECAL_Corr_post = new H2F("PCAL_ECAL_Corr_post",400,0,0.35, 400, 0,0.35);
    ArrayList<H2F> Histograms_Cut_Pre = new ArrayList<H2F>();
    ArrayList<H2F> Histograms_Cut_Aft = new ArrayList<H2F>();
     private String name;
	// TO CHANGE WITH THE REAL VALUES
     
     // 1 PCAL, 4 EC1 , 7 ECOUT
    // However transfered to my stuff":      
    private int PCALHIT=70; // I am taking the x,y,z,on PCAL 
	private int PCALID=71; //The index used to idenfity the PCAL
	private int ECAL1ID=74; //The index used to idenfity the ECAL1
	private int ECAL2ID=77; //The index used to idenfity the ECAL2 
	// ------------------------------------------------

	
	
// TO ADD sector by sector 
	Sampling_Fractions(){
		
	}
	

	public void setName(String nameDetector) {
		this.name=nameDetector;
	}
	public String toString() {
		return this.name;
	}
	
	/**
	 *  This function return the result of the cut 
	 *  case implemented: Electron and Photon 
	 */
	public boolean Status(ParticleREC particle) {
		//System.out.println("I am here in status EL EC CUT NO PID DEFINED");
		if (particle.pid()==11) { 
			return this.EC_sampling_fraction_cut(particle) ;
		//	EC_sampling_fraction_cut
		//	return this.ElectronCut(particle) ;
		}
		else if (particle.pid()==22) return this.PhotonCut(particle);
		// ADD PARTICLE 211, etc.
		else return false;
		
	}

	// --- Internal operation 
	private boolean PhotonCut(ParticleREC particle) {
		double EnergyFraction=0; 
		if(particle.hasDetector(PCALID)==true){ 
			EnergyFraction+= particle.getEnergy(PCALID);
		}
		if (particle.hasDetector(ECAL1ID) ==true) {
			EnergyFraction+= particle.getEnergy(ECAL1ID);
		}
		if (particle.hasDetector(ECAL2ID)== true) {
			EnergyFraction+= particle.getEnergy(ECAL2ID);
		}
		double Fraction = EnergyFraction/particle.p();
		SampliFraction_Overall_pre.fill(EnergyFraction,Fraction);
		// mettere settore per settore:
		//Hits=particle.getDetectorHits(PCALID);
		//int sector = getSector(Hits);
		
		if( Fraction >min_cut[0] && Fraction < max_cut[0] ) {
			
			SampliFraction_Overall_post.fill(EnergyFraction,Fraction);
			return true;}
		else return false;
	}

	
	boolean EC_sampling_fraction_cut(ParticleREC particle){
		
		double EnergyFraction=0; 
		double Epcal =0;
		double Einner=0;
		double Eout=0;
		int settore = particle.getSector(PCALID);
		if(particle.hasDetector(PCALID)==true){ 
			Epcal= particle.getEnergy(PCALID);
		}
		if (particle.hasDetector(ECAL1ID) ==true) {
			Einner=particle.getEnergy(ECAL1ID);
		}
		if (particle.hasDetector(ECAL2ID)== true) {
			Eout=particle.getEnergy(ECAL2ID);
		}
		
		double EnergyTotal=Epcal+Einner+Eout;

		  double [][]ecal_e_sampl_mu = {{  0.2531 ,  0.2550 ,  0.2514 ,  0.2494 ,  0.2528 ,  0.2521 },
		                                  { -0.6502 , -0.7472 , -0.7674 , -0.4913 , -0.3988 , -0.703  },
		                                  {  4.939  ,  5.350  ,  5.102  ,  6.440  ,  6.149  ,  4.957  }};

		  double [][] ecal_e_sampl_sigm = {{  2.726e-3 ,  4.157e-3 ,  5.222e-3 ,  5.398e-3 ,  8.453e-3 ,  6.533e-3 },
		                                    {  1.062    ,  0.859    ,  0.5564   ,  0.6576   ,  0.3242   ,  0.4423   },
		                                    { -4.089    , -3.318    , -2.078    , -2.565    , -0.8223   , -1.274    }};

		  double sigma_range = 3.5;

		  double mean = 0;
		  double sigma = 0;
		  double upper_lim_total = 0;
		  double lower_lim_total = 0;

		  for(int k = 0; k < 6; k++){  
		    if(settore-1 == k){
		      mean = ecal_e_sampl_mu[0][k] + ecal_e_sampl_mu[1][k]/1000*Math.pow(particle.p()-ecal_e_sampl_mu[2][k],2);
		      sigma = ecal_e_sampl_sigm[0][k] + ecal_e_sampl_sigm[1][k]/(10*(particle.p()-ecal_e_sampl_sigm[2][k]));
		      upper_lim_total = mean + sigma_range * sigma;
		      lower_lim_total = mean - sigma_range * sigma;
		    }
		  }

		  boolean pass_band = (EnergyTotal/particle.p()) <= upper_lim_total && (EnergyTotal/particle.p()) >= lower_lim_total;
		  boolean pass_triangle = false;

		  if(particle.p() < 4.5){ pass_triangle = true;}
		  else{pass_triangle = Einner/particle.p() > (0.2 -  Epcal/particle.p());}

		  if(pass_band && pass_triangle && Epcal>0.07) return true;
		 // if(pass_band && pass_triangle && Epcal>0.07) return true;
		  else return false;
		}
	
	
	
	
	
	
	private boolean ElectronCut(ParticleREC particle){
		
		double[][] ecal_e_sampl_mu = {{  0.2531 ,  0.2550 ,  0.2514 ,  0.2494 ,  0.2528 ,  0.2521 },
                { -0.6502 , -0.7472 , -0.7674 , -0.4913 , -0.3988 , -0.703  },
                {  4.939  ,  5.350  ,  5.102  ,  6.440  ,  6.149  ,  4.957  }};

       double[][] ecal_e_sampl_sigm = {{  2.726e-3 ,  4.157e-3 ,  5.222e-3 ,  5.398e-3 ,  8.453e-3 ,  6.533e-3 },
                  {  1.062    ,  0.859    ,  0.5564   ,  0.6576   ,  0.3242   ,  0.4423   },
                  { -4.089    , -3.318    , -2.078    , -2.565    , -0.8223   , -1.274    }};

double sigma_range = 3.5;
double mean = 0;
double sigma = 0;
double upper_lim_total = 0;
double lower_lim_total = 0;
//System.out.println("electron sector " + particle.getSector());

int sect =particle.getSector();
//int sect=  particle.getSector(PCALID);

//System.out.println("electron sector " + particle.getSector() + " dal hit PCALHIT "+ sect );


for(int k = 0; k < 6; k++){  
  if(sect -1  == k){
    mean = ecal_e_sampl_mu[0][k] + ecal_e_sampl_mu[1][k]/1000*Math.pow(particle.p()-ecal_e_sampl_mu[2][k],2);
    sigma = ecal_e_sampl_sigm[0][k] + ecal_e_sampl_sigm[1][k]/(10*(particle.p()-ecal_e_sampl_sigm[2][k]));
    upper_lim_total = mean + sigma_range * sigma;
    lower_lim_total = mean - sigma_range * sigma;
  }
}

		
		
		double EnergyFraction=0; 
		double Epcal =0;
		double Einner=0;
		double Eout=0;
		if(particle.hasDetector(PCALID)==true){ 
			EnergyFraction+= particle.getEnergy(PCALID);
			Epcal= particle.getEnergy(PCALID);
			//System.out.println("Epcal "+ Epcal);
		}
		if (particle.hasDetector(ECAL1ID) ==true) {
			EnergyFraction+= particle.getEnergy(ECAL1ID);
			Einner=particle.getEnergy(ECAL1ID);
		}
		if (particle.hasDetector(ECAL2ID)== true) {
			EnergyFraction+= particle.getEnergy(ECAL2ID);
			Eout=particle.getEnergy(ECAL2ID);
		}
		double Fraction = EnergyFraction/particle.p();
		SampliFraction_Overall_pre.fill(particle.p(),Fraction);
		SampliFraction_S1_pre.fill(particle.p(),Fraction);
		SampliFraction_S2_pre.fill(particle.p(),Fraction);
		SampliFraction_S3_pre.fill(particle.p(),Fraction);
		SampliFraction_S4_pre.fill(particle.p(),Fraction);
		
		boolean pass_band = false ;
		if (Fraction <= upper_lim_total && Fraction>= lower_lim_total) pass_band = true;
		

		// mettere settore per settore:
		//Hits=particle.getDetectorHits(PCALID);
		//int sector = getSector(Hits);
		double ratioPcal= Epcal/particle.p();
		double ratioEin= Einner/particle.p();
		double ratioEout = Eout/particle.p();
		PCAL_ECAL_Corr_pre.fill(ratioPcal, ratioEin);
		boolean pass_triangle=false; 
		if(particle.p() < 4.5){ pass_triangle = true;}
		  else{
			 if( pass_triangle = ratioEin > (0.2 - ratioPcal)) pass_triangle =true;
			  
		  }
		 if(pass_band &&  Epcal> 0.007 && (ratioPcal+ratioEin+ratioEout)>0.2  ) {
		//  if(pass_band && pass_triangle &&  Epcal> 0.007 && (ratioPcal+ratioEin+ratioEout)>0.2  ) {
			  SampliFraction_Overall_post.fill(particle.p(),Fraction);
				SampliFraction_S1_post.fill(particle.p(),Fraction);
				SampliFraction_S2_post.fill(particle.p(),Fraction);
				SampliFraction_S3_post.fill(particle.p(),Fraction);
				SampliFraction_S4_post.fill(particle.p(),Fraction);
				PCAL_ECAL_Corr_post.fill(ratioPcal, ratioEin);
			  return true;
		  }
		  else return false;
	}

	
	@Override
	/**
	 * This function return the sector of a given hit .
	 * @param Hit is a Vector3 representing x,y,z of a detector reconstructed cluster
	 **/
		public int getSector(Vector3 Hit) {
		//	System.out.println (" Calcolo dell' hit angolo" + Math.toDegrees(Hit.phi()));
			double myphi = Math.toDegrees(Hit.phi());
			int sec =0;
			if (myphi>0 && myphi<=30) sec= 1;
			else if (myphi>30 && myphi <=90) sec= 2;
			else if (myphi>90 && myphi <=150 ) sec= 3;
			else if (myphi>150 && myphi<= 180 ) sec=4;
			else if ( myphi>-30 && myphi<=0) sec=1;
			else if (myphi>-90 && myphi<=-30) sec=6;
			else if (myphi>-150 && myphi <= -90) sec=5;
			else if ( myphi>=-180 && myphi <=-150) sec= 4;
			return sec; 
			
			}

	
	public List<H2F> Histograms_Pre() {
		this.Histograms_Cut_Pre.add(SampliFraction_Overall_pre);
		this.Histograms_Cut_Pre.add(SampliFraction_S1_pre);
		this.Histograms_Cut_Pre.add(SampliFraction_S2_pre);
		this.Histograms_Cut_Pre.add(SampliFraction_S3_pre);
		this.Histograms_Cut_Pre.add(SampliFraction_S4_pre);
		return this.Histograms_Cut_Pre;
		// TODO Auto-generated method stub
	}
	public List<H2F> Histograms_Aft() {
		this.Histograms_Cut_Aft.add(SampliFraction_Overall_post);
		this.Histograms_Cut_Aft.add(SampliFraction_S1_post);
		this.Histograms_Cut_Aft.add(SampliFraction_S2_post);
		this.Histograms_Cut_Aft.add(SampliFraction_S3_post);
		this.Histograms_Cut_Aft.add(SampliFraction_S4_post);
		return this.Histograms_Cut_Aft;
		// TODO Auto-generated method stub
	}
	

	}

