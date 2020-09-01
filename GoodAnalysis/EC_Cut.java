package GoodAnalysis;

import java.util.ArrayList;
import java.util.List;

import org.jlab.groot.data.H2F;
import org.jlab.jnp.physics.Vector3;

public class EC_Cut implements DetectorCut{
	//Particle that doesn't survive the cuts
	H2F EC1Histo_Pre = new H2F("EC1_Pre",450,0,450,450,0,450);
	H2F EC2Histo_Pre = new H2F("EC2_Pre",450,0,450,450,0,450);
	H2F EC3Histo_Pre = new H2F("EC3_Pre",450,0,450,450,0,450);
	// Particle that does survive the cuts 
	H2F EC1Histo_Aft = new H2F("EC1_Aft",450,0,450,450,0,450);
	H2F EC2Histo_Aft = new H2F("EC2_Aft",450,0,450,450,0,450);
	H2F EC3Histo_Aft = new H2F("EC3_Aft",450,0,450,450,0,450);
    ArrayList<H2F> Histograms_Cut_Pre = new ArrayList<H2F>();
    ArrayList<H2F> Histograms_Cut_Aft = new ArrayList<H2F>();
     private String name;
	// THIS IS MY INTERNAL NOTATION USED IN PARTICLE REC 
    private int PCALHIT=70; // I am taking the x,y,z,on PCAL 
	private int PCALID=71; //The index used to idenfity the PCAL
	private int ECAL1ID=74; //The index used to idenfity the PCAL
	private int ECAL2ID=77; //The index used to idenfity the PCAL
	// ------------------------------------------------
	private double lu_cutMin,lv_cutMin,lw_cutMin,lu_cutMax,lv_cutMax,lw_cutMax;
	// Variables use internally
	private boolean PCAL_cut,EC1_cut,EC2_cut ;
	private Vector3 Hits;

	
// These are the definition of Stefan's Diehl
	  double min_v_tight_inb[] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
	  double min_v_med_inb[]   = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
	  double min_v_loose_inb[] = {9.0,  9.0,  9.0,  9.0,  9.0,  9.0 };
	  //
	  double max_v_tight_inb[] = {400, 400, 400, 400, 400, 400};
	  double max_v_med_inb[]   = {400, 400, 400, 400, 400, 400};
	  double max_v_loose_inb[] = {400, 400, 400, 400, 400, 400};
	  //
	  double min_w_tight_inb[] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
	  double min_w_med_inb[]   = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
	  double min_w_loose_inb[] = {9.0,  9.0,  9.0,  9.0,  9.0,  9.0 };
	  // 
	  double max_w_tight_inb[] = {400, 400, 400, 400, 400, 400};
	  double max_w_med_inb[]   = {400, 400, 400, 400, 400, 400};
	  double max_w_loose_inb[] = {400, 400, 400, 400, 400, 400};
	
// TO ADD sector by sector 
	EC_Cut(){
		//Implment to some ok values as inizialization
		this.lu_cutMin=20;this.lv_cutMin=20;this.lw_cutMin=20;this.lu_cutMax=400;this.lv_cutMax=400;this.lw_cutMax=400;
	}
	
	//Constructor with costum limits to apply to all calorimeters
	EC_Cut(double lum, double lvm, double lwm ,double luM, double lvM, double lwM){
		this.lu_cutMin=lum;this.lv_cutMin=lvm;this.lw_cutMin=lwm;this.lu_cutMax=luM;this.lv_cutMax=lvM;this.lw_cutMax=lwM;	
	}

	public void setName(String nameDetector) {
		this.name=nameDetector;
	}
	public String toString() {
		return this.name;
	}
	
	/**
	 * With this function the user can set new cuts for the variable lu,lv,lw
	 * @author gangelini
	 * @param new cuts value
	 * lu min, lv min , lw min, lu max, lv max ,lw max 
	 */
	public void SetNewLimits(double lu_min, double lv_min, double lw_min, double lu_max, double lv_max, double lw_max) {
		this.lu_cutMin=lu_min;this.lv_cutMin=lv_min;this.lw_cutMin=lw_min;this.lu_cutMax=lu_max;this.lv_cutMax=lv_max;this.lw_cutMax=lw_max;
	}
	/**
	 *  This function return the result of the cut 
	 *  case implemented: Electron and Photon
	 *  No other particles are implemented 
	 */
	public boolean Status(ParticleREC particle) {
		if (particle.pid()==11) { 
			return this.ElectronCut(particle) ;
		}
		else if (particle.pid()==22) return this.PhotonCut(particle);
		else return false;
		
	//	}
	}

	// --- Internal operation 
	private boolean PhotonCut(ParticleREC particle) {
		if(particle.hasDetector(PCALID)==true){ 
			Hits=particle.getDetectorHits(PCALID);
			return EC_cut(Hits.x(),Hits.y(),Hits.z());
		}
		else return false;
	}

	private boolean ElectronCut(ParticleREC particle){
		if(particle.hasDetector(PCALID)==true) {
			PCAL_cut=false;
			Hits=particle.getDetectorHits(PCALID);
			int sector = getSector(Hits);
		    //Setting Stefan's limits
			SetNewLimits(0,min_v_loose_inb[sector-1],min_w_loose_inb[sector-1],450,max_v_loose_inb[sector-1],max_w_loose_inb[sector-1]);
			if(particle.getDetectorHits(71).y()>min_v_loose_inb[sector-1] && max_v_loose_inb[sector-1]>particle.getDetectorHits(71).y()&& max_v_loose_inb[sector-1]>particle.getDetectorHits(71).z() && min_w_loose_inb[sector-1]<particle.getDetectorHits(71).z()	)
			{
					
return true;

		}
			else return false;
		}
		/*
	
		 if (particle.hasDetector(ECAL1ID)==true) {
			Hits=particle.getDetectorHits(ECAL1ID);
			int sector = getSector(Hits);
			SetNewLimits(0,0,0,450,450,450);
			//SetNewLimits(0,min_v_med_inb[sector-1],min_w_med_inb[sector-1],450,max_v_med_inb[sector-1],max_w_med_inb[sector-1]);
		//	SetNewLimits(0,min_v_med_inb[sector-1],min_w_med_inb[sector-1],450,max_v_med_inb[sector-1],max_w_med_inb[sector-1]);

			EC1_cut=EC_cut(Hits.x(),Hits.y(),Hits.z());
		}
		 if ( particle.hasDetector(ECAL2ID)==true) {
			Hits=particle.getDetectorHits(ECAL2ID);
			int sector = getSector(Hits);
			//SetNewLimits(0,min_v_med_inb[sector-1],min_w_med_inb[sector-1],450,max_v_med_inb[sector-1],max_w_med_inb[sector-1]);

			SetNewLimits(0,0,0,450,450,450);
			EC2_cut=EC_cut(Hits.x(),Hits.y(),Hits.z());
			//System.out.println("ECAL2: The hits x of the photon is " +Hits.get(0).x() +" The hits y of the photon is " +Hits.get(0).y()+" The hits z of the photon is " + Hits.get(0).z());
		}
		
		// The electron can have hits not only in PCAL
		if (particle.hasDetector(ECAL1ID)==false && particle.hasDetector(ECAL2ID)==false && PCAL_cut==true) return true; //hits only in PCAL
		else if (particle.hasDetector(ECAL2ID)==false && PCAL_cut==true && EC1_cut==true) return true; //hits in PCAL and EC1
		else if ( PCAL_cut==true && EC1_cut==true && EC2_cut==true) return true; // Hits in PCAL + EC1 +EC2 
		else if (particle.hasDetector(PCALID)==false&& EC1_cut==true && EC2_cut==true) return true;// Hits in EC1+EC2 but no signal in PCAL	
		else {
			//System.out.println("Non ho passato i tagli");
			return false;}
			*/
	//	if(PCAL_cut==true && particle.getDetectorHits(61).y()<14) System.out.println(" --> Siamo ubriachi");
		return PCAL_cut;
		
	}

	/**
	 * The class takes the lu,lv,lw coordinates from EC banks and apply the cut.
	 * @param d 
	 * @param e
	 * @param f
	 * @return boolean T if pass, F if not passed 
	 */
	private  boolean EC_cut (double d, double e, double f) {
		if (d>=this.lu_cutMin && e >= this.lv_cutMin && f >= this.lw_cutMin && d<= this.lu_cutMax && e<=this.lv_cutMax & f<=this.lw_cutMax) {
			//System.out.println("Il taglio era bono con minimo lv" + lv_cutMin);
			return true;
		}
		else {
			//System.out.println("Che schifo ! the velues of d, e and f are" + d + " " + e  + " " + f );
			return false ; 
		}
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
	@Override
	public List<H2F> Histograms_Pre() {
		this.Histograms_Cut_Pre.add(EC1Histo_Pre);
		this.Histograms_Cut_Pre.add(EC2Histo_Pre);
		this.Histograms_Cut_Pre.add(EC3Histo_Pre);
		return this.Histograms_Cut_Pre;
		// TODO Auto-generated method stub
	}
	public List<H2F> Histograms_Aft() {
		this.Histograms_Cut_Aft.add(EC1Histo_Aft);
		this.Histograms_Cut_Aft.add(EC2Histo_Aft);
		this.Histograms_Cut_Aft.add(EC3Histo_Aft);
		return this.Histograms_Cut_Aft;
		// TODO Auto-generated method stub
	}
	

	}

