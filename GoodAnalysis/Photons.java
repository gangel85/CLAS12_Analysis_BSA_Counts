package GoodAnalysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jlab.jnp.physics.LorentzVector;
import org.jlab.jnp.physics.Vector3;
//import org.jlab.clas.physics.Vector3;
//import org.jlab.geom.prim.Vector3D;
//import org.jlab.io.base.DataBank;
//import org.jlab.io.base.DataEvent;
import org.jlab.jnp.hipo.data.HipoEvent;
import org.jlab.jnp.hipo.data.HipoGroup;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;

public class Photons {
    private boolean Photon_Cuts=false;
	public float px,py,pz,ph_mom,ph_theta,ph_phi;
	public int PhotonsNr;
	public LorentzVector Vph;
	public double lu_1st,lv_1st,lw_1st;
	
	private ArrayList<LorentzVector> Photons_4Vect = new ArrayList<>();	
	public ArrayList <Float>  Photons_px= new ArrayList<Float>();
	public ArrayList <Float>  Photons_py= new ArrayList<Float>();
	public ArrayList <Float>  Photons_pz= new ArrayList<Float>();
	public ArrayList <Float>  Photons_mom= new ArrayList<Float>();
	public ArrayList <Float>  Photons_theta= new ArrayList<Float>();
	public ArrayList <Float>  Photons_phi= new ArrayList<Float>();
	public ArrayList <Float>  Photons_pxAC= new ArrayList<Float>();
	public ArrayList <Float>  Photons_pyAC= new ArrayList<Float>();
	public ArrayList <Float>  Photons_pzAC= new ArrayList<Float>();
	public ArrayList <Float>  Photons_momAC= new ArrayList<Float>();
	public ArrayList <Float>  Photons_thetaAC= new ArrayList<Float>();
	public ArrayList <Float>  Photons_phiAC= new ArrayList<Float>();
	public ArrayList <Float>  Photons_beta= new ArrayList<Float>();
	
	public double ph_ecal_XNC,ph_ecal_XAC,ph_ecal_YNC,ph_ecal_YAC,ph_ecal_Z_NC,ph_ecal_ZAC;
	
	
	public Photons(ParticleREC electron, Bank Rec_Part,  Bank Rec_Calorimeter, double startime)
	{  double vectAngles = 0; boolean EC_cut = false;

		PhotonsNr=0;
		Map<Integer,List<Integer>> caloMap = loadMapByIndex(Rec_Calorimeter,"pindex");
		int nrows = Rec_Part.getRows();
		for (int ipart =0 ; ipart<nrows; ipart++) {
		int  pid = Rec_Part.getInt("pid",ipart);
		int charge= (int) Rec_Part.getInt("charge",ipart);
		if(pid == 22){
			double time=0;
			double path=0;
			//Rec_Part.show();
		     px = Rec_Part.getFloat("px",ipart);
			 py = Rec_Part.getFloat("py",ipart);
			 pz = Rec_Part.getFloat("pz",ipart);
			 ph_mom = (float)Math.sqrt(px*px+py*py+pz*pz);
			 ph_theta = (float)Math.toDegrees(Math.acos(pz/ph_mom));
			 ph_phi = (float)Math.toDegrees(Math.atan2(py,px));		
			 Vph = new LorentzVector(px,py,pz,ph_mom);
			// System.out.println( " Phi Lorentz vecttor" + Math.toDegrees(Vph.phi()));
			 PhotonsNr++;
			 Photons_px.add(px);
			 Photons_py.add(py);
			 Photons_pz.add(pz);
			 Photons_mom.add(ph_mom);
			 Photons_theta.add(ph_theta);
			 Photons_phi.add(ph_phi);
			
			 //System.out.println(" Phi of my photon "+ ph_phi);
			//AnglePH_F.fill(ph_theta);
			//EnergiePH_F.fill(ph_mom1);

			//PH_angle =true;
			//System.out.println(" ======" );
			if(caloMap.get(ipart)!=null){
				 time = Rec_Calorimeter.getFloat("time",caloMap.get(ipart).get(0));
				 path = Rec_Calorimeter.getFloat("path",caloMap.get(ipart).get(0));
				int First_det = Rec_Calorimeter.getInt("layer", caloMap.get(ipart).get(0));
				//Rec_Calorimeter.show();
				int Second_det=0; int Third_det=0;
				if (caloMap.get(ipart).size()==2 )Second_det= Rec_Calorimeter.getInt("layer",caloMap.get(ipart).get(1));
				else if (caloMap.get(ipart).size()>2) {
					Second_det=  Rec_Calorimeter.getInt("layer",caloMap.get(ipart).get(1));
					Third_det= Rec_Calorimeter.getInt("layer",caloMap.get(ipart).get(2));
				//	System.out.println(" Sec " +Second_det+" Third " + Third_det);
				}
                  float ph_ecal_E=0,ph_ecal_X = 0,ph_ecal_Y = 0,ph_ecal_Z = 0;
                  boolean ph_EC1_cut = false,ph_EC2_cut = false,ph_EC3_cut = false;
                  int det=-5;
                  float lu,lv,lw =0;
				  for (int icalo : caloMap.get(ipart)) {
					
					         det =Rec_Calorimeter.getInt("layer",icalo);
					        //final int pind = Rec_Calorimeter.getInt("pindex",icalo);
							 ph_ecal_X = Rec_Calorimeter.getFloat("x",icalo);
							 ph_ecal_Y = Rec_Calorimeter.getFloat("y",icalo);
							 ph_ecal_Z = Rec_Calorimeter.getFloat("z",icalo);
							 ph_ecal_E += Rec_Calorimeter.getFloat("energy",icalo);
							 lu = Rec_Calorimeter.getFloat("lu",icalo);
							 lv = Rec_Calorimeter.getFloat("lv",icalo);
							 lw = Rec_Calorimeter.getFloat("lw",icalo);
							if(det==First_det){	 
								lu_1st=lu;
								lv_1st=lv;
								lw_1st=lw;
								Vector3 vector1 = new Vector3(electron.vector().vect().x(), electron.vector().vect().y(),electron.vector().vect().z());
								Vector3 vector2 = new Vector3(Vph.vect().x(),Vph.vect().y(),Vph.vect().z());
								//Vector3D vector1 = new Vector3D(electron.vector().vect().x(), electron.vector().vect().y(),electron.vector().vect().z());
								//Vector3D vector2 = new Vector3D(Vph.vect().x(),Vph.vect().y(),Vph.vect().z());
								 vectAngles = Math.toDegrees((vector1.dot(vector2))/vector1.mag()*vector2.mag());	
								// The electron vector is empty! Should study why 
							     ph_EC1_cut = new_EC_cut(lu,lv,lw);
							   }
							else if(det==Second_det &&det!=0) {
							ph_EC2_cut = new_EC_cut(lu,lv,lw);
							}
							else if(det==Third_det &&det!=0) {
							ph_EC3_cut = new_EC_cut(lu,lv,lw);
							}
							if(Second_det==0 && Third_det==0 && ph_EC1_cut==true)EC_cut=true;
							else if (Third_det==0 && ph_EC1_cut==true && ph_EC2_cut==true) EC_cut=true;
							else if (Second_det!=0&& Third_det!=0 && ph_EC1_cut==true && ph_EC2_cut==true && ph_EC3_cut==true) EC_cut=true;
				 }
		if(EC_cut==true) { ph_ecal_XAC=ph_ecal_X;ph_ecal_YAC=ph_ecal_Y;ph_ecal_ZAC=ph_ecal_Z;}
		else if(EC_cut!=true) {
			//System.out.println("This photon is wrong"); 
			ph_ecal_XNC=ph_ecal_X;ph_ecal_YNC=ph_ecal_Y;ph_ecal_Z_NC=ph_ecal_Z;}
		// It was 10 degree before 
		if (EC_cut==true && First_det!=7 && First_det!=4 && Math.toDegrees(Vph.vect().theta())>5 ) {
				//if (EC_cut==true && vectAngles>4 && First_det!=7 && First_det!=4 && Math.toDegrees(Vph.vect().theta())>10) {
					//System.out.println(" Qua non ci sono ");
				//if (EC_cut==true && vectAngles>2) {
					//System.out.println(" The detector is " + det);
					//H_ECal_xy.fill(ph_ecal_X, ph_ecal_Y);
			
			//System.out.println( " EC CUTS ARE TRUE ");
					Photon_Cuts=true;
					Photons_4Vect.add(Vph);   	
					Photons_pxAC.add(px);
					Photons_pyAC.add(py);
					Photons_pzAC.add(pz);
					Photons_momAC.add(ph_mom);
					Photons_thetaAC.add(ph_theta);
					Photons_phiAC.add(ph_phi);
				double beta=(path/(time-startime))/29.9792;
					Photons_beta.add((float) beta);
			
					}
	    }
				
	   //   }
	
	 }
	
	}
	}
	
	private  Map<Integer,List<Integer>> loadMapByIndex(	Bank rec_Calorimeter,String idxVarName) {
		Map< Integer,List<Integer> > map = new HashMap <Integer, List<Integer> >();
		if (rec_Calorimeter!=null) {
			for (int iFrom=0; iFrom<rec_Calorimeter.getRows(); iFrom++) {
				//System.out.println(iFrom +" index loop ");
				final int iTo = rec_Calorimeter.getInt(idxVarName,iFrom);
			//	System.out.println(iTo + " iTo value for string"+ idxVarName);
				if (!map.containsKey(iTo)) map.put(iTo,new ArrayList<Integer>()); 
				map.get(iTo).add(iFrom);
			}
		}
		return map;
	}
	
	public static boolean new_EC_cut (float lu, float lv, float lw) {
		if (lu>=34 && lv >= 14 && lw >= 14 && lu<= 405 && lv<=405 & lw<=405) return true;
		else return false ; 
	}
	
	public ArrayList<LorentzVector> getPhoton4Vects() {
		return Photons_4Vect;
	}
	
	public boolean getStatus(int i) {
		if(this.PhotonsNr>=i) {
		return Photon_Cuts;}
		else return false;
	}
	
	}
