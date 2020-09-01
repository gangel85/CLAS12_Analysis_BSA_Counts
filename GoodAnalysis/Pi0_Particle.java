package GoodAnalysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.jlab.jnp.physics.LorentzVector;
import org.jlab.jnp.physics.Vector3;
import org.jlab.jnp.physics.reaction.TransMatrix;
import org.jlab.groot.data.H1F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.TGCanvas;



public class Pi0_Particle {
	private	double Sector1_B=-30,Sector1_E=30,Sector2_B=30, Sector2_E=90, Sector3_B=90, Sector3_E=150, Sector4_B=150, Sector4_E=-150, Sector5_B=-150,Sector5_E=-90, Sector6_B=-90, Sector6_E=-30;

	private int PhiB0= -180,PhiB1= -150; private int PhiB2= -120;  private int PhiB3= -90;  private int PhiB4= -60; private int PhiB5= -30; private int PhiB6= 0;private int PhiB7= 30; private int PhiB8= 60; private int PhiB9= 90; private int PhiB10= 120;  private int PhiB11= 150; private int PhiB12= 180; 
private double minMass = 0.01;
private double maxMass= 0.015;
	private double Phi_LAB=99;
	private int nrPi0=0;
	
	private double Xf=0;
	
	static ArrayList<LorentzVector> Photons4Vectors = new ArrayList<>();
	static ArrayList<Double> Pi0_Z = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_Pt = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_Phi = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_Mass = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_XF = new ArrayList<Double> ();
	static ArrayList<LorentzVector> Pi0_4Vector = new ArrayList<LorentzVector> ();
	static ArrayList<LorentzVector> Pi0_X4Vector = new ArrayList<LorentzVector> ();
	
private double minPhi=999;
private double maxPhi=-999;

	int Method =0 ; 
	double EnergyCut=0;
	 public int N0_B1=0,N0_B2=0,N0_B3=0,N0_B4=0,N0_B5=0,N0_B6=0,N0_B7=0,N0_B8=0,N0_B9=0,N0_B10=0,N0_B11=0,N0_B0=0;
	 public int N1_B1=0,N1_B2=0,N1_B3=0,N1_B4=0,N1_B5=0,N1_B6=0,N1_B7=0,N1_B8=0,N1_B9=0,N1_B10=0,N1_B11=0,N1_B0=0;

	public void setMethod(int k) {
		if (k<=2) {
			this.Method= k;}
		else this.Method=0;
	}
	public void setEnergyCut(double En) {
		this.EnergyCut= En;
	}

	public  Pi0_Particle() {
		this.Method=0;
		this.EnergyCut=0;
	}

	//	public void Neutral_Pion_Analysis(ElectronSelection electron, double Degrees) {
	    public void getPi0s(LorentzVector beam, LorentzVector target, LorentzVector vecE,LorentzVector vecQ2, LorentzVector vecW2, double xB, double Q2, double Degrees) {
		this.Phi_LAB=99;
		this.nrPi0=0;
		this.Pi0_Z.clear(); this.Pi0_Pt.clear(); this.Pi0_Phi.clear(); this.Pi0_Mass.clear(); this.Pi0_4Vector.clear();
		double Angl_Cut = Math.toRadians(Degrees);
		boolean AngleCheck=true; 
		double Pq = -target.px()*vecQ2.px()-target.py()*vecQ2.py()-target.pz()*vecQ2.pz()+target.e()*vecQ2.e();
		LorentzVector q = new LorentzVector(0,0,0,0);
		double modulus =  (Math.sqrt(vecQ2.px()*vecQ2.px()+vecQ2.py()*vecQ2.py()+vecQ2.pz()*vecQ2.pz()));
		q.setPxPyPzE(vecQ2.px()/modulus, vecQ2.py()/modulus,vecQ2.pz()/modulus, vecQ2.e());
		if(Method==0 || EnergyCut==0 ) { System.out.println(" Attention the Pi0 extrapolation method should be set and the Energy min cut should be set");return;}
		else if( Method==1) {
			double sector1 =0, sector2=0;
			double PPh=0, Z_var=0,Angolo=0, AngleAA=0,AngleSin=0;
			double Pi0M=0,Pls = 0, Pt=0, Plod = 0;
			double amp2=0.88035;
			for( int kk=0 ; kk<Photons4Vectors.size()-1; kk ++)
			{
				sector1=getSector(Photons4Vectors.get(kk).vect()); 
				for(int qq=kk+1 ; qq< Photons4Vectors.size(); qq++) 
				{   
					//System.out.println(" I fotoni li ho ");
					sector2= getSector(Photons4Vectors.get(qq).vect()); // check the two photons are in the same sector
					//double angolo = Photons4Vectors.get(kk).vect().dot(Photons4Vectors.get(qq).vect()) ;
					//System.out.println("Angolo Cut " + Angl_Cut + " Angolo between photon case 2 : " + angolo );
					if(Photons4Vectors.get(kk).vect().dot(Photons4Vectors.get(qq).vect()) <=Angl_Cut ) {
						double angolo2 = Photons4Vectors.get(kk).vect().dot(Photons4Vectors.get(qq).vect()) ;
						System.out.println("Angolo Cut " + Angl_Cut + " Angolo between photon case 2 : " + angolo2 );
						AngleCheck=false;}
					else if(Photons4Vectors.get(qq).vect().dot(Photons4Vectors.get(kk).vect() )<= Angl_Cut){
						double angolo2 = Photons4Vectors.get(qq).vect().dot(Photons4Vectors.get(kk).vect());
						System.out.println("Angolo Cut " + Angl_Cut + " Angolo between photon case 2 : " + angolo2 );
						AngleCheck=false;
						}
					
					if(Photons4Vectors.get(kk).e()>EnergyCut  && Photons4Vectors.get(qq).e()>EnergyCut && AngleCheck==true && sector1== sector2 ) {
						//System.out.println("I am here !!! Hiii! ");
						LorentzVector PHM = new LorentzVector(0,0,0,0);
						PHM.add(Photons4Vectors.get(kk));
						PHM.add(Photons4Vectors.get(qq));
						PPh = -target.px()*PHM.px()-target.py()*PHM.py()-target.pz()*PHM.pz()+target.e()*PHM.e();
					//	double Z_new=PHM.e()/(beam.e()-vecE.e());
						Z_var= PPh/Pq;
						
						
						//System.out.println( " New Z " +Z_new + " old Z "+ Z_var);
					
						LorentzVector Vphotons = new LorentzVector(0,0,0,0);
						LorentzVector Gammacm = new LorentzVector(0,0,0,0);
						LorentzVector Protoncm = new LorentzVector(0,0,0,0);
						LorentzVector PIcm = new LorentzVector(0,0,0,0); 
						LorentzVector Beamcm = new LorentzVector(0,0,0,0);
						LorentzVector X= new LorentzVector(0,0,0,0);X.add(target).add(beam).sub(vecE).sub(PHM);
						LorentzVector VGS = new LorentzVector(0,0,0,0); VGS.add(beam).sub(vecE); //THe proper virtual photon
						Vector3 VectorPh = VGS.vect(); // Virtual photon is beam-el scattered but adding the target doesent change tri momenta
						Vector3 VectorPi = PHM.vect();
						Vector3 VectorP = target.vect();
						Vector3 VectorB=beam.vect();
						Vector3 VectorPhNorm = new Vector3(VectorPh.x()/VectorPh.mag(),VectorPh.y()/VectorPh.mag(),VectorPh.z()/VectorPh.mag());				
						Vector3 Angl1= VectorPhNorm.cross(VectorB);
						Vector3 Angl2= VectorPhNorm.cross(VectorPi);
						Angl1.setXYZ(Angl1.x()/Angl1.mag(), Angl1.y()/Angl1.mag(),Angl1.z()/Angl1.mag());
						Angl2.setXYZ(Angl2.x()/Angl2.mag(), Angl2.y()/Angl2.mag(),Angl2.z()/Angl2.mag());
						double AngoloCos=Angl1.dot(Angl2);
						Angolo= Math.toDegrees(Math.acos(AngoloCos));
						//this.Phi_LAB = Math.toDegrees(Math.acos(Angolo/(Angl1.mag()*Angl2.mag())));
						
						Vector3 AngleA = VectorB.cross(VectorPi);
						AngleAA= AngleA.dot(VectorPhNorm);
						Vector3 AngleB=VectorPhNorm.cross(VectorB);
						Vector3 AngleBB =VectorPhNorm.cross(VectorPi);
						AngleSin= AngleAA/(AngleB.mag()*AngleBB.mag());
						double AlternativePhi= Math.toDegrees(Math.asin(AngleSin));
						//System.out.println("----------");
						//System.out.println(" Phi cos " +Angl1.dot(Angl2) + " sin " + AngleSin );
						this.Phi_LAB = Angolo;
						//if(AngleSin<0) this.Phi_LAB=-this.Phi_LAB; 
						//System.out.println("Angolo LAB [degrees] " + Angolo);
						//if(AlternativePhi<0) AlternativePhi=-AlternativePhi; 
						//System.out.println(" Angolo LAB Alternativo " + AlternativePhi);
					
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
						
						//NUOVO ANGOLO
						//System.out.println( " Angolo nuovo " + Math.toDegrees(Math.atan2(AngleSin,AngoloCos)) + " versus " + this.Phi_LAB);
						
						//System.out.println("ProtonCM" +Protoncm.vect().x() + " y " + Protoncm.vect().y() + " z " + Protoncm.vect().z());
						//System.out.println("Photon CM " +Gammacm.vect().x() + " y " + Gammacm.vect().y() + " z " + Gammacm.vect().z());

						 Pi0M= PHM.mass();
						 Pls = PIcm.vect().dot(Gammacm.vect())/Gammacm.vect().mag();
						 //this is pt2 
						 Pt=PIcm.vect().mag2()-Pls*Pls;
						 Plod = PHM.vect().dot(VGS.vect())/VGS.vect().mag();
						// Maybe wrong
						 double Xf = (2 * Pls) / vecW2.mass();
						
						 
						 if ( Pls>0) {
						Pi0_Z.add(Z_var);
						Pi0_Pt.add(Pt);
						Pi0_Phi.add(this.Phi_LAB);
						Pi0_4Vector.add(PHM);
						Pi0_Mass.add(PHM.mass());
						Pi0_X4Vector.add(X);
						Pi0_XF.add(Xf);
						this.nrPi0++;					
						 }
						 else if (Pls<0) System.out.println(" PLS is negative and XF is " + Xf);
						/* Cosa di Harut 
						 * Not working so far*/
					      
						
					     double EHS=(vecW2.mass2()+PHM.mass2()-X.mass2())/(2*vecW2.mass());    // Energy of the Hadron in the center of mass	
					     // System.out.println(" Energy CM Harut vs Boost" + EHS +" vs "+ (PIcm.e())) ; 
					       double P2=(0.25*Math.pow(vecW2.mass2()-(target.e()*target.e())-Q2, 2)+(target.e()*target.e())-Q2)/(vecW2.mass2()); // Momentum squared of Hadron center of mass
			             // System.out.println(" Momentum squared hadron CM" + P2 + " VS  calculated by me " + PiCM.mag2());
					       double PlS=( PHM.e()*beam.e()-EHS*Math.sqrt(P2+(target.e()*target.e())))/Math.sqrt(P2);
			            //  System.out.println("Energy"+ PHM.e()+ " PLS " + PlS);
					  //	double XfH=2*PlS/vecW2.mass();
					  	double PTH=(EHS*EHS)-PHM.mass2()-(PlS*PlS);
						//System.out.println("PTH 1 " +(EHS*EHS) +" Second piece "+ PHM.mass2());
					  	/*
							if(X.mass()>1.3 && Z_var>0.2 && Z_var<0.9 ){
							Z_Dist_Overall.fill(Z_var);
							XB_Dist_Overall.fill(xB);
							Pt_Dist_Overall.fill(Pt);
							if( Z_var <= 0.1 ) { 
							InvMassPlots.get(0).fill(PHM.mass());
							Z_Dist_1.fill(Z_var);
							XB_Dist_1.fill(xB);
							Pt_Dist_1.fill(Pt);
							}
							else if (Z_var>0.1 && Z_var<=0.2) {InvMassPlots.get(1).fill(PHM.mass());
							Z_Dist_2.fill(Z_var);
							XB_Dist_2.fill(xB);
							Pt_Dist_2.fill(Pt);
							}
							else if (Z_var>0.2 && Z_var<=0.3) {InvMassPlots.get(2).fill(PHM.mass()); 
							Z_Dist_3.fill(Z_var);
							XB_Dist_3.fill(xB);
							Pt_Dist_3.fill(Pt);
							}
							else if (Z_var>0.3 && Z_var<=0.4) {InvMassPlots.get(3).fill(PHM.mass());
							Z_Dist_4.fill(Z_var);
							XB_Dist_4.fill(xB);
							Pt_Dist_4.fill(Pt);
							}
							else if (Z_var>0.4 && Z_var<=0.5) {
							InvMassPlots.get(4).fill(PHM.mass()); 
							Z_Dist_5.fill(Z_var);
							XB_Dist_5.fill(xB);
							Pt_Dist_5.fill(Pt);
							}
							else if (Z_var>0.5 && Z_var<=0.6) {InvMassPlots.get(5).fill(PHM.mass());
							Z_Dist_6.fill(Z_var);
							XB_Dist_6.fill(xB);
							Pt_Dist_6.fill(Pt);
							}
							else if (Z_var>0.6 && Z_var<=0.7) {InvMassPlots.get(6).fill(PHM.mass()); 
							Z_Dist_7.fill(Z_var);
							XB_Dist_7.fill(xB);
							Pt_Dist_7.fill(Pt);
							}
							else if (Z_var>0.7 && Z_var<=0.8) {InvMassPlots.get(7).fill(PHM.mass()); 
							Z_Dist_8.fill(Z_var);
							XB_Dist_8.fill(xB);
							Pt_Dist_8.fill(Pt);
							}
							else if (Z_var>0.8 && Z_var<=0.9) {InvMassPlots.get(8).fill(PHM.mass()); 
							Z_Dist_9.fill(Z_var);
							XB_Dist_9.fill(xB);
							Pt_Dist_9.fill(Pt);
							}
							else if (Z_var>0.9 && Z_var<=1)   {InvMassPlots.get(9).fill(PHM.mass());
							Z_Dist_10.fill(Z_var);
							XB_Dist_10.fill(xB);
							Pt_Dist_10.fill(Pt);
							}
						}
							*/
					}					
				}
			}
		}
		
		else if( Method==2) {		
			double Energy_max=0;
			LorentzVector Reference = new LorentzVector();
			LorentzVector MaxPh= new LorentzVector();
			LorentzVector SndPh= new LorentzVector();
			int Indice=-1;	
			int F1_sect=0;
			int F2_sect=0;
			// No vertex correction 
			for (int ll=0;ll < Photons4Vectors.size();ll++)
			{
				if(Photons4Vectors.get(ll).e() > Energy_max) { 
					Energy_max=Photons4Vectors.get(ll).e();
					MaxPh.copy(Photons4Vectors.get(ll));
					Indice=ll;
				}
			}
			Energy_max=0;
			for (int ii=0;ii<Photons4Vectors.size();ii++)
			{
				if(ii!=Indice)
				{	
					if( Photons4Vectors.get(ii).vect().dot(MaxPh.vect()) <= Angl_Cut ) AngleCheck=false;
					else if( MaxPh.vect().dot(Photons4Vectors.get(ii).vect()) <= Angl_Cut ) AngleCheck=false;
					if(Photons4Vectors.get(ii).e() > Energy_max &&Photons4Vectors.get(ii)!=Photons4Vectors.get(Indice) && AngleCheck==true) { 
						Energy_max=Photons4Vectors.get(ii).e();
						SndPh.copy(Photons4Vectors.get(ii));
					}
				}
			}
			Vector3 fotone1 = new Vector3();
			Vector3 fotone2 = new Vector3();
			// fotone1.setXYZ(MaxPh.vect().x(), MaxPh.vect().y(), MaxPh.vect().z()-e_vz);
			//  fotone2.setXYZ(SndPh.vect().x(), SndPh.vect().y(), SndPh.vect().z()-e_vz);
			fotone1.setXYZ(MaxPh.vect().x(), MaxPh.vect().y(), MaxPh.vect().z());
			fotone2.setXYZ(SndPh.vect().x(), SndPh.vect().y(), SndPh.vect().z());
			double Angle_ph=(fotone1.dot(fotone2))/(fotone1.mag()*fotone2.mag());
			//double Angle_ph=((MaxPh.vect()).dot(SndPh.vect()))/(MaxPh.vect().mag()*SndPh.vect().mag());
			LorentzVector PHM = new LorentzVector(0,0,0,0);
			PHM.add(MaxPh);
			PHM.add(SndPh);		
			double PPh = -target.px()*PHM.px()-target.py()*PHM.py()-target.pz()*PHM.pz()+target.e()*PHM.e();
			double Z_var= PPh/Pq;
			LorentzVector Vphotons = new LorentzVector(0,0,0,0);

			//check normalization+
			//Vphotons.setPxPyPzE(vecW2.px()/vecW2.mass(), vecW2.py()/vecW2.mass(),vecW2.pz()/vecW2.mass(), vecW2.e()/vecW2.mass());
			double Scalar_old = PHM.vect().dot(vecW2.vect())/vecW2.vect().mag();
			double Pt_old=Math.sqrt(PHM.vect().mag2()-Scalar_old*Scalar_old);

			//Z_Plot.fill(Z_var);

			if( MaxPh.e()> EnergyCut && SndPh.e()>EnergyCut ) {
				LorentzVector PHM2 = new LorentzVector(0,0,0,0);
				PHM2.add(MaxPh);
				PHM2.add(SndPh);		
				LorentzVector Vphotons2 = new LorentzVector(0,0,0,0);
				LorentzVector Gammacm = new LorentzVector(0,0,0,0);
				LorentzVector Protoncm = new LorentzVector(0,0,0,0);
				LorentzVector PIcm = new LorentzVector(0,0,0,0); 
				LorentzVector X= new LorentzVector(0,0,0,0);X.add(target).add(beam).sub(vecE).sub(PHM);
				LorentzVector VGS = new LorentzVector(0,0,0,0); VGS.add(beam).sub(vecE); //THe proper virtual photon
				//System.out.println("VecW2 x,y,z: " +vecW2.vect().x()+" . "+vecW2.vect().y()+" , "+  vecW2.vect().z());
				//	System.out.println("VecVGS x,y,z: " +VGS.vect().x()+" . "+VGS.vect().y()+" , "+  VGS.vect().z());
				Vector3 VectorPh = VGS.vect(); // Virtual photon is beam-el scattered but adding the target doesent change tri momenta
				Vector3 VectorPi=PHM2.vect();
				Vector3 VectorP = target.vect();
				Vector3 VectorB=beam.vect();
				Vector3 VectorPhNorm = new Vector3(VectorPh.x()/VectorPh.mag(),VectorPh.y()/VectorPh.mag(),VectorPh.z()/VectorPh.mag());
				Vector3 Angl1=  VectorPhNorm.cross(VectorB)	;
				Vector3 Angl2=  VectorPhNorm.cross(VectorPi)	;	
				double Angolo=  Angl1.dot(Angl2);
				this.Phi_LAB = Math.toDegrees(Math.acos(Angolo/(Angl1.mag()*Angl2.mag())));

				Vector3 rotationVector = new Vector3(0,0,0);
				TransMatrix rotationMatrix = new TransMatrix();
				rotationMatrix.compose(VectorPh);
				Vector3  PhCM= rotationMatrix.mult(VectorPh);
				Vector3 PiCM= rotationMatrix.mult(VectorPi);
				Vector3 PCM = rotationMatrix.mult(VectorP);
				Gammacm.setPxPyPzE(PhCM.x(), PhCM.y(), PhCM.z(), VGS.e());
				Protoncm.setPxPyPzE(PCM.x(), PCM.y(), PCM.z(), target.e());
				PIcm.setPxPyPzE(PiCM.x(), PiCM.y(), PiCM.z(), PHM.e());
				Vphotons.add(Gammacm).add(Protoncm);
				Vector3 boost = Vphotons.boostVector();
				boost.negative();
				Vphotons.boost(boost);
				Gammacm.boost(boost);
				Protoncm.boost(boost);
				PIcm.boost(boost);
				double Pi0M= PHM2.mass();
				double Pls = PIcm.vect().dot(Gammacm.vect())/Gammacm.vect().mag();
				double Pt=Math.sqrt(PIcm.vect().mag2()-Pls*Pls);
				double Plod = PHM2.vect().dot(VGS.vect())/VGS.vect().mag();
				Pi0_Z.add(Z_var);
				Pi0_Pt.add(Pt);
				Pi0_Phi.add(Angolo);
				Pi0_4Vector.add(PHM);
				Pi0_Mass.add(PHM.mass());
				
				this.nrPi0++;
				//Pi_0_M2_4_BK.fill(PHM2.mass());
				//Pi_0_M2_4.fill(PHM2.mass());
				if(Math.toDegrees(fotone1.phi())<Sector1_E && Math.toDegrees(fotone1.phi())>=Sector1_B) F1_sect=1;
				else if(Math.toDegrees(fotone1.phi())<Sector2_E && Math.toDegrees(fotone1.phi())>=Sector2_B) F1_sect=2;
				else if(Math.toDegrees(fotone1.phi())<Sector3_E && Math.toDegrees(fotone1.phi())>=Sector3_B) F1_sect=3;
				else if(Math.toDegrees(fotone1.phi())<Sector4_E || Math.toDegrees(fotone1.phi())>=Sector4_B) F1_sect=4;
				else if(Math.toDegrees(fotone1.phi())<Sector5_E && Math.toDegrees(fotone1.phi())>=Sector5_B) F1_sect=5;
				else if(Math.toDegrees(fotone1.phi())<Sector6_E && Math.toDegrees(fotone1.phi())>=Sector6_B) F1_sect=6;    
				if(Math.toDegrees(fotone2.phi())<Sector1_E && Math.toDegrees(fotone2.phi())>=Sector1_B) F2_sect=1;
				else if(Math.toDegrees(fotone2.phi())<Sector2_E && Math.toDegrees(fotone2.phi())>=Sector2_B) F2_sect=2;
				else if(Math.toDegrees(fotone2.phi())<Sector3_E && Math.toDegrees(fotone2.phi())>=Sector3_B) F2_sect=3;
				else if(Math.toDegrees(fotone2.phi())<Sector4_E || Math.toDegrees(fotone2.phi())>=Sector4_B) F2_sect=4;
				else if(Math.toDegrees(fotone2.phi())<Sector5_E && Math.toDegrees(fotone2.phi())>=Sector5_B) F2_sect=5;
				else if(Math.toDegrees(fotone2.phi())<Sector6_E && Math.toDegrees(fotone2.phi())>=Sector6_B) F2_sect=6;    
				//	if(F1_sect==1&F2_sect==1) Pi0M_S1.fill(PHM2.mass());
				//	if(F1_sect==2&F2_sect==2) Pi0M_S2.fill(PHM2.mass());
				//	if(F1_sect==3&F2_sect==3) Pi0M_S3.fill(PHM2.mass());
				//	if(F1_sect==4&F2_sect==4) Pi0M_S4.fill(PHM2.mass());
				//	if(F1_sect==5&F2_sect==5) Pi0M_S5.fill(PHM2.mass());
				//	if(F1_sect==6&F2_sect==6) Pi0M_S6.fill(PHM2.mass());
				if(X.mass()>1.3){
					
					double VPr= (PHM2.vect().dot(q.vect())/q.vect().mag());
					Vector3 P_l = new Vector3(VPr*q.vect().x(),VPr*q.vect().y(),VPr*q.vect().z() ) ;
					Vector3 P_t = new Vector3( PHM2.px()-P_l.x(),PHM2.py()-P_l.y(),PHM2.pz()-P_l.z());
					//  System.out.println(P_t.mag());	
					//	H_p_Pt.fill(P_t.mag());
					//H_p_Z.fill(Z_var);
				}
			}
		}
	}

	public double getPhiLab()
	{
		return this.Phi_LAB;
	}
	public void add(Photons photons)
	{
		for(int i =0  ; i< photons.getPhoton4Vects().size();i++) {
			Photons4Vectors.add(photons.getPhoton4Vects().get(i));
		}
	}
	public int getNrPi0s()
	{
		return this.nrPi0;
	}
	public int getNrPhotons()
	{
		return this.Photons4Vectors.size();
	}
public ArrayList<LorentzVector> getLorentzVector(){
	return this.Pi0_4Vector;
}

public ArrayList<LorentzVector> getInvMassLVector(){
	return this.Pi0_X4Vector;
}
public ArrayList<Double> getZ()
{
	return this.Pi0_Z;
}
public ArrayList<Double> getPt()
{
	return this.Pi0_Pt;
}
 
public ArrayList<Double> getXF()
{
	return this.Pi0_XF;
}

public int GetIndex_MostEnergy()
{
	int indice = -9;
	double Energy =0;
	
	for (int i =0 ; i< this.Pi0_4Vector.size(); i++)
	{
		if( this.Pi0_4Vector.get(i).e()> Energy) {
			Energy=this.Pi0_4Vector.get(i).e();
			indice=i; 
		}
	}
	return indice;
}

/**
 * 
 * @return List of pi0s angle in the event in radiant measured from lab frame
 * */
public ArrayList<Double> getPhi()
{
	return this.Pi0_Phi;
}

public ArrayList<Double> getMass()
{
	return this.Pi0_Mass;
}
	

	public void resetPhotons()
	{
		this.Photons4Vectors.clear();
	}
	public void printsHCounts() {
		System.out.println(" Angle min and max" + this.minPhi  +" " +this.maxPhi );
		System.out.println("H0 B1 " + N0_B1 + " H1 B1 " + N1_B1);
		System.out.println("H0 B2 " + N0_B2 + " H1 B2 " + N1_B2);
		System.out.println("H0 B3 " + N0_B3 + " H1 B3 " + N1_B3);
		System.out.println("H0 B4 " + N0_B4 + " H1 B4 " + N1_B4);
		System.out.println("H0 B5 " + N0_B5 + " H1 B5 " + N1_B5);
		System.out.println("H0 B6 " + N0_B6 + " H1 B6 " + N1_B6);
		System.out.println("H0 B7 " + N0_B7 + " H1 B7 " + N1_B7);
		System.out.println("H0 B8 " + N0_B8 + " H1 B8 " + N1_B8);
		System.out.println("H0 B9 " + N0_B9 + " H1 B9 " + N1_B9);
		System.out.println("H0 B10 " + N0_B10 + " H1 B10 " + N1_B10);
		System.out.println("H0 B11 " + N0_B11 + " H1 B11 " + N1_B11);
		// TODO Auto-generated method stub
		
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
