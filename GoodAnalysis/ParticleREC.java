/**
 * 
 */
package GoodAnalysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jlab.jnp.hipo.data.HipoEvent;
import org.jlab.jnp.hipo.data.HipoGroup;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.physics.Particle;
import org.jlab.jnp.physics.Vector3;

/**
 * @author gangelini
 *Check how to get the sector without using the DC. Now I take it from the first DC hit when I go into the Fiducial cuts 
 */
public class ParticleREC extends Particle {
	// The map has for key a detector ID and for value a list of vector3 representing the position of the hits within that detector.
	// Note the EC vector is not composed by x,y,z but lu,lv,lw
	//Map<Integer, List<Vector3>> ParticleHits = new HashMap<Integer, List<Vector3>>();
	Map<Integer, Vector3> ParticleHits = new HashMap<Integer, Vector3>();
	Map<Integer, Double> ParticleEnergy = new HashMap<Integer, Double>();
	Map<Integer, Integer> ParticleSector = new HashMap<Integer, Integer>();
	//Info that I use for the printing 
	//public double vX=-9999, vY=-9999, vZ=-9999, Theta=-9999, Phi=-9999, Mom=-9999;
	// To discuss about EC with someone experts.
	public int sector=-1; // Sector 
	private int DetKey=99999; //this value means the key is unknown
	// ID DETCTORS:
	int DC1=61; int DC2=62; int DC3=63; //Detector 6: DC . Layers 1,2,3
	int PCAL=71; int EC1=74; int EC2=77;  // Detector 7: ECAL  . Layers 1,4,7
	int PCAL_layer=1, EC1_layer=4, EC2_layer=7;
	int DC_detector=6 ;
	// OLD:
//	int DC_layer1 =6 , DC_layer2=18, DC_layer3=30 ;  
	// NEW:
	int DC_layer1 =6 , DC_layer2=18, DC_layer3=36 ; 
	private double mom =0;
	private double Theta =0 ;
	private double Phi =0;
	private double vX=-999;
	private double vY= -999;
	private double vZ=-999;
	private boolean CalCheck=false;
	private boolean DCCheck=false;


	@Override
	public double p() {
		// TODO Auto-generated method stub
		return this.mom;
	}
	@Override
	public double theta() {
		// TODO Auto-generated method stub
		return this.Theta;
	}

	@Override
	public double phi() {
		// TODO Auto-generated method stub
		return this.Phi;
	}

	@Override
	public double vx() {
		// TODO Auto-generated method stub
		return this.vX;
	}
	@Override
	public double vy() {
		// TODO Auto-generated method stub
		return this.vY;
	}
	@Override
	public double vz() {
		// TODO Auto-generated method stub
		return this.vZ;
	}
	
	public ParticleREC(Particle particle,Bank Rec_Part,Bank Rec_Traj,Bank Rec_Calorimeter) {
		//System.out.println("Ciao I am in partREC with particle id" + particle.pid());
		ParticleHits.clear(); // Rest the list containing the hits of a particle, so that each new events will have an empty list.
		ParticleEnergy.clear(); // Same as particle hits, but for the energies deposited.
		ParticleSector.clear();
		// Filling standard info inherited from the Particle object. 
		this.vX=particle.vx(); this.vY=particle.vy(); this.vZ=particle.vz();
		this.mom=particle.p();this.Theta=particle.theta();this.Phi=particle.phi();
		int pid = particle.pid();		int charge = particle.charge();
		//System.out.println(" PID is " + pid);

		boolean PCAL_Hit = false;
		float lu,lv,lw; // Hits on the calorimeter
		float dx,dy,dz;// Hits on DC
		Map<Integer,List<Integer>> caloMap = loadMapByIndex(Rec_Calorimeter,"pindex");
		Map<Integer,List<Integer>> trajMap = loadMapByIndex(Rec_Traj,"pindex");
		//int nrows = Rec_Part.getNode("pid").getDataSize();
		//for (int ipart =0 ; ipart<nrows; ipart++) {
		
		// ----> Electron Specific Information 
		if (pid ==11) {
			if(caloMap.get(0)!=null){  
			//	Rec_Calorimeter.show();
				CalCheck=true;
				for (int icalo : caloMap.get(0)) {
					//System.out.println("Hello");
					//Rec_Calorimeter.show();
					int det = Rec_Calorimeter.getInt("layer", icalo);
					int sector = Rec_Calorimeter.getInt("sector", icalo);
					lu = Rec_Calorimeter.getFloat("lu",icalo);
					lv = Rec_Calorimeter.getFloat("lv",icalo);
					lw = Rec_Calorimeter.getFloat("lw",icalo);
					double energy = Rec_Calorimeter.getFloat("energy",icalo);
					//Devo fare qualcosa cosi qui come PCAL_Hit
					if(det==PCAL_layer){	 	
						PCAL_Hit=true;
						//if( lv<20)System.out.println("Hi");
						Vector3 HitUVW= new Vector3(lu,lv,lw);
						ParticleHits.put(PCAL,HitUVW);
						ParticleEnergy.put(PCAL,energy);	
						ParticleSector.put(PCAL,sector);
				}
					else if(det==EC1_layer) {
						Vector3 HitUVW= new Vector3(lu,lv,lw);
						ParticleHits.put(EC1,HitUVW);	
						ParticleEnergy.put(EC1,energy);	
						ParticleSector.put(EC1,sector);
					}
					else if(det==EC2_layer) {
						Vector3 HitUVW= new Vector3(lu,lv,lw);
						ParticleHits.put(EC2,HitUVW);
						ParticleEnergy.put(EC2,energy);
						
					}
			}
		}
			boolean checkHits=false;
			if(charge!=0 && trajMap.get(0)!=null ){
				for(int index =0 ; index<trajMap.get(0).size(); index ++ ) {
					DCCheck=true;
					//Rec_Traj.show();
					
					
					if( Rec_Traj.getInt("detector",trajMap.get(0).get(index)) == DC_detector && Rec_Traj.getInt("layer",trajMap.get(0).get(index))== DC_layer1){
						//System.out.println( " Layer e' " + Rec_Traj.getInt("layer",trajMap.get(ipart).get(index)));
						dx = Rec_Traj.getFloat("x",trajMap.get(0).get(index));
						dy = Rec_Traj.getFloat("y",trajMap.get(0).get(index));
						dz = Rec_Traj.getFloat("z",trajMap.get(0).get(index));
						Vector3 HitXYZ= new Vector3(dx,dy,dz);
						ParticleHits.put(DC1,HitXYZ);
						ParticleSector.put(DC1,determineSector(HitXYZ));
						checkHits=true;
					}
					else if( Rec_Traj.getInt("detector",trajMap.get(0).get(index)) == DC_detector && Rec_Traj.getInt("layer",trajMap.get(0).get(index))== DC_layer2){
						dx = Rec_Traj.getFloat("x",trajMap.get(0).get(index));
						dy = Rec_Traj.getFloat("y",trajMap.get(0).get(index));
						dz = Rec_Traj.getFloat("z",trajMap.get(0).get(index));
						Vector3 HitXYZ2= new Vector3(dx,dy,dz);
						ParticleHits.put(DC2,HitXYZ2);
						ParticleSector.put(DC2,determineSector(HitXYZ2));
						checkHits=true;
					}

					else if( Rec_Traj.getInt("detector",trajMap.get(0).get(index)) == DC_detector && Rec_Traj.getInt("layer",trajMap.get(0).get(index))== DC_layer3){
						dx = Rec_Traj.getFloat("x",trajMap.get(0).get(index));
						dy = Rec_Traj.getFloat("y",trajMap.get(0).get(index));
						dz = Rec_Traj.getFloat("z",trajMap.get(0).get(index));
						Vector3 HitXYZ3= new Vector3(dx,dy,dz);
						ParticleHits.put(DC3,HitXYZ3); 		
						ParticleSector.put(DC3,determineSector(HitXYZ3));
						checkHits=true;
					}
				}			
			}
			else {
				Vector3 Empty = new Vector3 (0,0,0);
				ParticleHits.put(DC1,Empty);
				ParticleHits.put(DC2,Empty);
				ParticleHits.put(DC3,Empty);
			}
		}
		//In the case I don't have an electron
		else {
			//System.out.println(" ------- >> Ao ma che fai? Questo non e' un elettrone");
		for (int ipart =0 ; ipart<Rec_Part.getRows(); ipart++) {
			if(pid == Rec_Part.getInt("pid", ipart)) {
				if(caloMap.get(ipart)!=null){  
					CalCheck=true;
					//Am I sure about this? 
					for (int icalo : caloMap.get(ipart)) {
					//	Rec_Calorimeter.show();
						int sector = Rec_Calorimeter.getInt("sector", icalo);
						int det = Rec_Calorimeter.getInt("layer", icalo);
						lu = Rec_Calorimeter.getFloat("lu",icalo);
						lv = Rec_Calorimeter.getFloat("lv",icalo);
						lw = Rec_Calorimeter.getFloat("lw",icalo);
						double energy = Rec_Calorimeter.getFloat("energy",icalo);
						if(det==PCAL_layer){	 	
							PCAL_Hit=true;
							Vector3 HitUVW= new Vector3(lu,lv,lw);
							ParticleHits.put(PCAL,HitUVW);
							ParticleEnergy.put(PCAL,energy);
							ParticleSector.put(PCAL,sector);

						}
						else if(det==EC1_layer) {
							Vector3 HitUVW= new Vector3(lu,lv,lw);
							ParticleHits.put(EC1,HitUVW);
							ParticleEnergy.put(EC1,energy);	
							ParticleSector.put(EC1,sector);

						   
						}
						else if(det==EC2_layer) {
							Vector3 HitUVW= new Vector3(lu,lv,lw);
							ParticleHits.put(EC2,HitUVW);
							ParticleEnergy.put(EC2,energy);	
							ParticleSector.put(EC2,sector);

							//ParticleHits.put(EC2,new ArrayList<Vector3>());
							//ParticleHits.get(EC2).add(HitUVW);   
						}
					}
				}
				boolean checkHits=false;
		
				if(charge!=0 && trajMap.get(ipart)!=null ){
					for(int index =0 ; index<trajMap.get(ipart).size(); index ++ ) {
						DCCheck=true;
						if( Rec_Traj.getInt("detector",trajMap.get(ipart).get(index)) == DC_detector && Rec_Traj.getInt("layer",trajMap.get(ipart).get(index))== DC_layer1){
							dx = Rec_Traj.getFloat("x",trajMap.get(ipart).get(index));
							dy = Rec_Traj.getFloat("y",trajMap.get(ipart).get(index));
							dz = Rec_Traj.getFloat("z",trajMap.get(ipart).get(index));
							Vector3 HitXYZ= new Vector3(dx,dy,dz);
							ParticleHits.put(DC1,HitXYZ);
							ParticleSector.put(DC1,determineSector(HitXYZ));
						//	ParticleHits.put(DC1,new ArrayList<Vector3>());
						//	ParticleHits.get(DC1).add(HitXYZ);
							checkHits=true;
							//System.out.println(" HIT " + HitXYZ.getXYZString());
						//	System.out.println(" HIT da ParticleHits " +ParticleHits.get(DC1).get(0).getXYZString());
						}
						else if( Rec_Traj.getInt("detector",trajMap.get(ipart).get(index)) == DC_detector && Rec_Traj.getInt("layer",trajMap.get(ipart).get(index))== DC_layer2){

							dx = Rec_Traj.getFloat("x",trajMap.get(ipart).get(index));
							dy = Rec_Traj.getFloat("y",trajMap.get(ipart).get(index));
							dz = Rec_Traj.getFloat("z",trajMap.get(ipart).get(index));
							Vector3 HitXYZ2= new Vector3(dx,dy,dz);
							ParticleHits.put(DC2,HitXYZ2);
							ParticleSector.put(DC2,determineSector(HitXYZ2));
							//ParticleHits.put(DC2,new ArrayList<Vector3>());
							//ParticleHits.get(DC2).add(HitXYZ2);
							checkHits=true;
						}

						else if( Rec_Traj.getInt("detector",trajMap.get(ipart).get(index)) == DC_detector && Rec_Traj.getInt("layer",trajMap.get(ipart).get(index))== DC_layer3){

							dx = Rec_Traj.getFloat("x",trajMap.get(ipart).get(index));
							dy = Rec_Traj.getFloat("y",trajMap.get(ipart).get(index));
							dz = Rec_Traj.getFloat("z",trajMap.get(ipart).get(index));
							Vector3 HitXYZ3= new Vector3(dx,dy,dz);
							ParticleHits.put(DC3,HitXYZ3); 
							ParticleSector.put(DC3,determineSector(HitXYZ3));
							//ParticleHits.put(DC3,new ArrayList<Vector3>());
							//ParticleHits.get(DC3).add(HitXYZ3); 
							checkHits=true;
						}
					
		/*			
				
else {
							System.out.println("ciao");
							Vector3 Empty = new Vector3 (0,0,0);
							ParticleHits.put(DC1,Empty);
							ParticleHits.put(DC2,Empty);
							ParticleHits.put(DC3,Empty);
							/*
							ParticleHits.put(DC1,new ArrayList<Vector3>());
							ParticleHits.get(DC1).add(Empty);   	
							ParticleHits.put(DC2,new ArrayList<Vector3>());
							ParticleHits.get(DC2).add(Empty);   	
							ParticleHits.put(DC3,new ArrayList<Vector3>());
							ParticleHits.get(DC3).add(Empty);   	
							
						}
						
*/
					}
				
					
				
				}
				else {
					//Rec_Traj.show();
					System.out.println("Cazzo! ho le banche DC am non ho hit??");
					System.out.println(" Charge is + " + charge );
					System.out.println( " Traj Map " + trajMap.get(ipart));
					Vector3 Empty = new Vector3 (0,0,0);
					ParticleHits.put(DC1,Empty);
					ParticleHits.put(DC2,Empty);
					ParticleHits.put(DC3,Empty);
				}
				/*
				if ( checkHits==false ) {System.out.println("Cazzo! ho le banche DC am non ho hit??");
				Vector3 Empty = new Vector3 (0,0,0);
				ParticleHits.put(DC1,Empty);
				ParticleHits.put(DC2,Empty);
				ParticleHits.put(DC3,Empty);
				//ParticleHits.put(DC1,new ArrayList<Vector3>());
				//ParticleHits.get(DC1).add(Empty);   	
				//ParticleHits.put(DC2,new ArrayList<Vector3>());
				//ParticleHits.get(DC2).add(Empty);   	
				//ParticleHits.put(DC3,new ArrayList<Vector3>());
				//ParticleHits.get(DC3).add(Empty);   	
				}
*/
			}
		}
		}
	}

	public Map getHitsMap() {
		return ParticleHits;
	}
	public Map getEnergyMap() {
		return ParticleEnergy;
	}
	public Map getSectorMap() {
		return ParticleSector;
	}
	public boolean checkDetectorID(int pCALID) {
		// TODO Auto-generated method stub
		return false;
	}

	// How can I improve this ? 
	// Ask GAgik
	public boolean hasDetector(int pCALID) {
		for (Integer key : ParticleHits.keySet()) {
			// To check what this stuff does 
			if (key==pCALID) { this.DetKey=key; return true;}
		}
		return false;
	}
/*
	public List<Vector3> getDetectorHits(int pCALID)
	{   
		if (this.DetKey==99999) {
			for (Integer key : this.ParticleHits.keySet()) { 
				if (key==pCALID) this.DetKey=key;
			}
		}
		int Key=DetKey;
		this.DetKey=99999; //restore the key to "unknown"
		//System.out.println("I am getting ParticleHits with get ID" + Key);
		return this.ParticleHits.get(Key);

	}
*/
	
	public Vector3 getDetectorHits(int pCALID)
	{   
		if (this.DetKey==99999) {
			for (Integer key : this.ParticleHits.keySet()) { 
				if (key==pCALID) this.DetKey=key;
			}
		}
		int Key=DetKey;
		this.DetKey=99999; //restore the key to "unknown"
		//System.out.println("I am getting ParticleHits with get ID" + Key);
		return this.ParticleHits.get(Key);

	}
	
	public double getEnergy(int pCALID) {
		if (this.DetKey==99999) {
			for (Integer key : this.ParticleEnergy.keySet()) { 
				if (key==pCALID) this.DetKey=key;
			}
		}
		int Key=DetKey;
		this.DetKey=99999; //restore the key to "unknown"
		//System.out.println("I am getting ParticleEnergy with get ID" + Key);
		return this.ParticleEnergy.get(Key);
		
	}
	
	public int getSector( int pCALID) {
		if (this.DetKey==99999) {
			for (Integer key : this.ParticleSector.keySet()) { 
				if (key==pCALID) this.DetKey=key;
			}
		}
		int Key=DetKey;
		this.DetKey=99999; //restore the key to "unknown"
		//System.out.println("I am getting ParticleEnergy with get ID" + Key);
		return this.ParticleSector.get(Key);
		
	}
	
	
	
	public void setSector(int sect)
	{
		this.sector=sect;
	}
 public  int getSector()
 {
	 return this.sector;
 }
	/**
	 * @param fromBank the bank containing the index variable
	 * @param idxVarName the name of the index variable
	 * @return map with keys being the index in toBank and values the indices in fromBank
	 */
	private  Map<Integer,List<Integer>> loadMapByIndex(	Bank rec_Calorimeter,String idxVarName) {
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

	public boolean Status()
	{
		if(this.CalCheck==true && this.DCCheck==true )	{	
			return true;
		}
		//else return false;
		else return false;
	}


	 int determineSector(Vector3 Hit){

	   double phi = 180 / Math.PI* Math.atan2(Hit.y() /
	Math.sqrt(Math.pow(Hit.x(), 2) + Math.pow(Hit.y(), 2) +
	Math.pow(Hit.z(), 2)), Hit.x() / Math.sqrt(Math.pow(Hit.x(), 2)
	                + Math.pow(Hit.y(), 2) + Math.pow(Hit.z(), 2)));

	   if(phi < 30 && phi >= -30){        return 1;}
	   else if(phi < 90 && phi >= 30){    return 2;}
	   else if(phi < 150 && phi >= 90){   return 3;}
	   else if(phi >= 150 || phi < -150){ return 4;}
	   else if(phi < -90 && phi >= -150){ return 5;}
	   else if(phi < -30 && phi >= -90){  return 6;}

	   return 0;
	 }
	
}
