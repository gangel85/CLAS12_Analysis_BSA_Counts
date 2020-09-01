package GoodAnalysis;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import org.jlab.jnp.physics.Particle;
import org.jlab.jnp.physics.LorentzVector;
import org.jlab.jnp.physics.Vector3;
import org.jlab.jnp.physics.reaction.TransMatrix;
import org.jlab.groot.data.H1F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.TGCanvas;


public class ChargedPi_Particle {
	private	double Sector1_B=-30,Sector1_E=30,Sector2_B=30, Sector2_E=90, Sector3_B=90, Sector3_E=150, Sector4_B=150, Sector4_E=-150, Sector5_B=-150,Sector5_E=-90, Sector6_B=-90, Sector6_E=-30;
	private int PhiB0= -180,PhiB1= -150; private int PhiB2= -120;  private int PhiB3= -90;  private int PhiB4= -60; private int PhiB5= -30; private int PhiB6= 0;private int PhiB7= 30; private int PhiB8= 60; private int PhiB9= 90; private int PhiB10= 120;  private int PhiB11= 150; private int PhiB12= 180; 
	private double minMass = 0.01;
	private double maxMass= 0.015;
	private double Phi_LAB=99;
	private int nrPi0=0;
	static ArrayList<LorentzVector> Photons4Vectors = new ArrayList<>();
	private LorentzVector Pi_4Vector= new LorentzVector();
	static double Pi_Z, Pi_Pt, Pi_Phi,Pi_Mass;
	static double Pi_PhiClas; 
	private double  Pi_XMass ;
	private double minPhi=999;
	private double maxPhi=-999;
	private double Pi_Energy; 
	private double Pi_Angle;
	private double Pi_Momentum;
    private int sector;
	 int Method =0 ; 
	 double EnergyCut=0;
	 public int N0_B1=0,N0_B2=0,N0_B3=0,N0_B4=0,N0_B5=0,N0_B6=0,N0_B7=0,N0_B8=0,N0_B9=0,N0_B10=0,N0_B11=0,N0_B0=0;
	 public int N1_B1=0,N1_B2=0,N1_B3=0,N1_B4=0,N1_B5=0,N1_B6=0,N1_B7=0,N1_B8=0,N1_B9=0,N1_B10=0,N1_B11=0,N1_B0=0;
	private double XF;

	public void setMethod(int k) {
		if (k<=2) {
			this.Method= k;}
		else this.Method=0;
	}
	public void setEnergyCut(double En) {
		this.EnergyCut= En;
	}

	public  ChargedPi_Particle() {
		this.Method=0;
		this.EnergyCut=0;
	}
	
	    public void getChargedPi(Particle pion, LorentzVector beam, LorentzVector target, LorentzVector vecE, LorentzVector vecQ2, LorentzVector vecW2, double xB, double Q2, double Degrees) {
		this.Phi_LAB=99;
		this.nrPi0=0;
		double Angl_Cut = Math.toRadians(Degrees);
		boolean AngleCheck = true; 
		double Pq = -target.px()*vecQ2.px()-target.py()*vecQ2.py()-target.pz()*vecQ2.pz()+target.e()*vecQ2.e();
		LorentzVector q = new LorentzVector(0,0,0,0);
		double modulus =  (Math.sqrt(vecQ2.px()*vecQ2.px()+vecQ2.py()*vecQ2.py()+vecQ2.pz()*vecQ2.pz()));
		q.setPxPyPzE(vecQ2.px()/modulus, vecQ2.py()/modulus,vecQ2.pz()/modulus, vecQ2.e());
			double sector1 =0, sector2=0;
			double PPh = 0, Z_var=0,Angolo=0, AngleAA=0,AngleSin=0;
			double Pi0M = 0,Pls = 0, Pt=0, Plod = 0;
			double amp2 = 0.88035;
				
				LorentzVector PHM = new LorentzVector();// PHM is the lorentz vector of the pion
				this.Pi_Momentum=Math.sqrt(pion.px()*pion.px()+pion.py()*pion.py()+pion.pz()*pion.pz());
				
			//	double PionEnergy=Math.sqrt(pion.px()*pion.px()+pion.py()*pion.py()+pion.pz()*pion.pz()+0.13957*0.13957);
			//	System.out.println( " Pion energy from particle rec e()" + pion.e()+ " pion energy from calculation "+ PionEnergy);
				PHM.setPxPyPzE(pion.px(), pion.py(), pion.pz(), pion.e());
				PPh = -target.px()*PHM.px()-target.py()*PHM.py()-target.pz()*PHM.pz()+target.e()*PHM.e();
				Z_var= PPh/Pq;
				LorentzVector Vphotons = new LorentzVector(0,0,0,0);
				LorentzVector Gammacm = new LorentzVector(0,0,0,0);
				LorentzVector Protoncm = new LorentzVector(0,0,0,0);
				LorentzVector PIcm = new LorentzVector(0,0,0,0); 
				LorentzVector Beamcm = new LorentzVector(0,0,0,0);
				LorentzVector X = new LorentzVector(0,0,0,0);X.add(target).add(beam).sub(vecE).sub(PHM);
				LorentzVector VGS = new LorentzVector(0,0,0,0); VGS.add(beam).sub(vecE); //THe proper virtual photon
				Vector3 VectorPh = VGS.vect(); // Virtual photon is beam-el scattered but adding the target doesent change tri momenta
				Vector3 VectorPi = PHM.vect();
				Vector3 VectorP = target.vect();
				Vector3 VectorB = beam.vect();
				this.sector=getSector(VectorPi);
				this.Pi_PhiClas=Math.toDegrees(VectorPi.phi());
				Vector3 VectorPhNorm = new Vector3(VectorPh.x()/VectorPh.mag(),VectorPh.y()/VectorPh.mag(),VectorPh.z()/VectorPh.mag());				
				Vector3 Angl1 = VectorPhNorm.cross(VectorB);
						Vector3 Angl2 = VectorPhNorm.cross(VectorPi);
						Angl1.setXYZ(Angl1.x()/Angl1.mag(), Angl1.y()/Angl1.mag(),Angl1.z()/Angl1.mag());
						Angl2.setXYZ(Angl2.x()/Angl2.mag(), Angl2.y()/Angl2.mag(),Angl2.z()/Angl2.mag());
						double AngoloCos=Angl1.dot(Angl2);
						Angolo= Math.toDegrees(Math.acos(AngoloCos));
						//this.Phi_LAB = Math.toDegrees(Math.acos(Angolo/(Angl1.mag()*Angl2.mag())));
						Vector3 AngleA = VectorB.cross(VectorPi);
						AngleAA = AngleA.dot(VectorPhNorm);
						Vector3 AngleB=VectorPhNorm.cross(VectorB);
						Vector3 AngleBB =VectorPhNorm.cross(VectorPi);
						AngleSin= AngleAA/(AngleB.mag()*AngleBB.mag());
						double AlternativePhi= Math.toDegrees(Math.asin(AngleSin));
						
						this.Phi_LAB = Angolo;
						
						if(Angolo>0 && AlternativePhi>0) this.Phi_LAB=Angolo;
						else if(Angolo<0 && AlternativePhi>0) this.Phi_LAB=Angolo;
						else if (Angolo<0  && AlternativePhi<0) this.Phi_LAB=Angolo+90;
						else if (Angolo>0 && AlternativePhi<0 ) this.Phi_LAB=360-Angolo;
						//System.out.println("Angolo Ph e' " + this.Phi_LAB );
						// Nuovo Angolo ORLANDO
						//double AnglePHI_Orlando1= (VectorPh.cross(VectorB)).dot(VectorPi);
						//AnglePHI_Orlando1 = AnglePHI_Orlando1/Math.abs(AnglePHI_Orlando1);
						//System.out.println( "Angolo Orlando " + (AnglePHI_Orlando1*Angolo));
						Vector3 rotationVector = new Vector3(0,0,0);
						TransMatrix rotationMatrix = new TransMatrix();
						//Capire come funziona sta roba!
						
						rotationMatrix.compose(VectorPh);
						Vector3 PhCM = rotationMatrix.mult(VectorPh);
		//System.out.println(" Hello all " + " x "+ PhCM.x()+" y "+PhCM.y()+ " z "+PhCM.z());
						
						Vector3 PiCM = rotationMatrix.mult(VectorPi);
						Vector3 PCM = rotationMatrix.mult(VectorP);
						Vector3 BCM= rotationMatrix.mult(VectorB);
						Gammacm.setPxPyPzE(PhCM.x(), PhCM.y(), PhCM.z(), VGS.e());
						//Check these energie! 
						Protoncm.setPxPyPzE(PCM.x(), PCM.y(), PCM.z(), target.e());
						//System.out.println("PhCM" + PCM.x() + " y " + PCM.y() + " z " + PCM.z());
						PIcm.setPxPyPzE(PiCM.x(), PiCM.y(), PiCM.z(), PHM.e());
						//System.out.println("PHM e "+PHM.e()+ " Calculated " +((0.135*0.135)+PiCM.mag2()));
						Beamcm.setPxPyPzE(BCM.x(), BCM.y(), BCM.z(), beam.e());
						//System.out.println(" Energia beam " + Beamcm.e());
						Vphotons.add(Gammacm).add(Protoncm);
						Vector3 boost = Vphotons.boostVector();
						//System.out.println(" Boost x" + boost.x() + " Boost y"+ boost.y()+" boost z"+ boost.z());
						//boost.negative();
						Vphotons.boost(boost);
						Gammacm.boost(boost);
						Protoncm.boost(boost);
						PIcm.boost(boost);
						Beamcm.boost(boost);
						
						 Pi0M= PHM.mass();
						 Pls = PIcm.vect().dot(Gammacm.vect())/Gammacm.vect().mag();
					
						
						 //this is pt2
						 Pt = PIcm.vect().mag2()-Pls*Pls; // QUesto e' PT^2 
						 Plod = PHM.vect().dot(VGS.vect())/VGS.vect().mag();
						 if(Math.toDegrees(PIcm.angle(Gammacm))>90) System.out.println( "Attenzione all'angolo " + Math.toDegrees(PIcm.angle(Gammacm)));
						
						
						 // Cose nuove:
						 LorentzVector Vettore= new LorentzVector();
						 LorentzVector NewGammaCm = new LorentzVector() ;
						 NewGammaCm.setPxPyPzE(VGS.px(), VGS.py(), VGS.pz(), VGS.e());
						 LorentzVector NewPionCm= new LorentzVector();
						 NewPionCm.setPxPyPzE(pion.px(), pion.py(), pion.pz(), pion.e());
						Vettore.add(VGS).add(target);
						Vector3 boostNew = Vettore.boostVector();
						NewGammaCm.boost(boostNew);
						NewPionCm.boost(boostNew);						
						Double New_Pls = NewPionCm.vect().dot(NewGammaCm.vect())/NewGammaCm.vect().mag();
						Double New_Pt  = NewPionCm.vect().mag2()-New_Pls*New_Pls;
						
						
						 Pi_Z= Z_var;
						 
						 // Misuro Pt
						 Pi_Pt=Math.sqrt(Pt);
						 Pi_Phi = Phi_LAB;
						 Pi_Mass=PHM.mass();
						 this.Pi_XMass=X.mass(); //invariant mass
						 Pi_Energy= pion.e();
						 Pi_Angle= Math.toDegrees(pion.theta());
						 Pi_4Vector=PHM;
						 double myXf= 2*Pls/vecW2.mass();
						 this.XF = myXf;
						
						/* Cosa di Harut 
						 * Not working so far*/
					      
						
					     double EHS=(vecW2.mass2()+PHM.mass2()-X.mass2())/(2*vecW2.mass());    // Energy of the Hadron in the center of mass	
					     //if (xB>0.4) System.out.println(" Energy CM Harut vs Boost" + EHS +" vs "+ (PIcm.e())) ; 
					       double P2=(0.25*Math.pow(vecW2.mass2()-(target.e()*target.e())-Q2, 2)+(target.e()*target.e())-Q2)/(vecW2.mass2()); // Momentum squared of Hadron center of mass
					      // if (xB>0.4) System.out.println(" Momentum squared hadron CM" + P2 + " VS  calculated by me " + PiCM.mag2());
					  
					       double PlS=( PHM.e()*beam.e()-EHS*Math.sqrt(P2+(target.e()*target.e())))/Math.sqrt(P2);
			          
					       //  System.out.println("Energy"+ PHM.e()+ " PLS " + PlS);
					       
					  	
					       double XfH=2*PlS/vecW2.mass(); //XF of the Hadron
					      /// if (xB>0.4) System.out.println("I am in charged partcile, W is " + vecW2.mass() + " PLS "+ PlS + " XF " +XfH);
					      // if (xB>0.4) System.out.println(" XF my new "+myXf );
					       //	this.XF = XfH;
					  	double PTH=(EHS*EHS)-PHM.mass2()-(PlS*PlS);
						//System.out.println("PTH 1 " +(EHS*EHS) +" Second piece "+ PHM.mass2());

		
	}
	    
	    public int getSector() {
	    	return this.sector;
	    }
	 public double getMomentum() {
		 return this.Pi_Momentum;
	 }

	 public double getXF() {
		 return this.XF;
	 }
	public double getPhiLab()
	{
		return this.Phi_LAB;
	}
	
	public double getTheta()
	{
		return this.Pi_Angle;
	}
	
	public double getEnergy()
	{
		return this.Pi_Energy;
	}
	
	public int getNrPi0s()
	{
		return this.nrPi0;
	}
	
public LorentzVector getLorentzVector(){
	return this.Pi_4Vector;
}

public double getMissMass(){
	return this.Pi_XMass;
}
public double getZ()
{
	return this.Pi_Z;
}
public double getPt()
{
	return this.Pi_Pt;
}

public double getPhiClas12() {
	return this.Pi_PhiClas;
}
/**
 * 
 * @return List of pi0s angle in the event in radiant measured from lab frame
 * */
public double getPhi()
{
	return this.Pi_Phi;
}

public double getMass()
{
	return this.Pi_Mass;
}


	public int getSector(Vector3 Hit) {
	//	System.out.println (" Calcolo dell' hit angolo" + Math.toDegrees(Hit.phi()));
		// remove 30 to put sector 1 on the O, add 180 to transform everything positively.
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
	
}
