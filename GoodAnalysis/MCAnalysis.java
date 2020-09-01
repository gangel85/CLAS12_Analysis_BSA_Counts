package GoodAnalysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.H3F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.jnp.hipo.data.HipoEvent;
import org.jlab.jnp.hipo.data.HipoGroup;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.physics.LorentzVector;
import org.jlab.jnp.physics.Vector3;
import org.jlab.jnp.physics.reaction.TransMatrix;

public class MCAnalysis {
	double z1,z2,z3,z4,z5,z6,z7,z8,z9,z10;
	double pt1,pt2,pt3,pt4,pt5,pt6,pt7,pt8,pt9,pt10;
	private int nr_particle;
	double beamEnergy = 10.6;
	private LorentzVector VB   = new LorentzVector(0.0,0.0,beamEnergy,beamEnergy); // Getting the Lorentz Vector associated to the beam
	private LorentzVector VT = new LorentzVector(0.0,0.0,0.0,0.93827); // Defining Lorentz Vector for the proton target

	private int parentID;
	private FileWriter fileWriter;
	private  double min_pion= 0.05; private  double max_pion=0.25;	private int Mass_Bins=100;
	private int Z_Bins=0; private double Z_min=0.2; private double Z_max=0.8; 
	private int PT_Bins=0;private double PT_min=0.1; private double PT_max=1; 
	private H2F InvMass,H_e_vx_vy;
	private H1F H_e_vz,H_e_Mom, H_e_theta, H_e_phi;
	private H2F ElBinMigration;
	private double Q2_el_MC= 0;
	private double XB_el_MC=0; 
	private double y_el_MC=0;
	private double W_el_MC=0;
	
	private String FirstVariable, SecondVariable, ThirdVariable, ForthVariable  ;
	int FirstBins, SecondBins, ThirdBins, ForthBins;
	double FirstMin, FirstMax,SecondMin,SecondMax,ThirdMin, ThirdMax, ForthMin, ForthMax;

	
	// Variables of the electron:
	
	public double e_xB, e_Q2, e_W, e_theta, e_phi, e_p;
	
   private int CountsID=0;
   private int PID=0;

	
	
	private MultiBins Hadron_Lund,Hadron_Z,Hadron_PT, Hadron_Phi;
	
	public ArrayList <MultiBins>  Particles_LUND= new ArrayList<MultiBins>();
	
	private ArrayList<LorentzVector> Particle_4Vector = new ArrayList<>();	
	public ArrayList <Float>  Hadron_en= new ArrayList<Float>();
	public ArrayList <Float>  Hadron_mom= new ArrayList<Float>();
	public ArrayList <Float>  Hadron_theta= new ArrayList<Float>();
	public ArrayList <Float>  Hadron_phi= new ArrayList<Float>();
	public ArrayList <Float> Hadron_phiTrento= new ArrayList<Float>();
	public ArrayList <Float>  Hadron_z= new ArrayList<Float>();
	public ArrayList <Float>  Hadron_pt= new ArrayList<Float>();
	public ArrayList <Integer> Hadron_parent = new ArrayList<Integer>();
	public ArrayList <Integer> Hadron_pid = new ArrayList<Integer>();

	public MCAnalysis( MultiBins ParticleMC, int zbin, double zmin, double zmax, int ptbin, double ptmin, double ptmax, int phibins)
	{   

		
		this.Z_Bins=zbin; this.PT_Bins=ptbin;
		InvMass = new H2F("INVMASS",zbin, zmin, zmax, ptbin, ptmin, ptmax);
		H_e_vz = new H1F("H_e_vz","H_e_vz",200,-20,20);
		H_e_vz.setTitle("LUND Electron longitudinal vertex ");
		H_e_vz.setTitleX("v_{z} (cm)");
		H_e_vz.setFillColor(2);
		H_e_vx_vy = new H2F("H_e_vx_vy ", "H_e_vx_vy ", 100, -20,20, 100, -20, 20 );
		H_e_vx_vy.setTitle("LUND Electron vertex x vs y ");
		H_e_vx_vy.setTitleX(" El. x vertex [cm]");
		H_e_vx_vy.setTitleY(" El. y vertex [cm]");
		H_e_Mom = new H1F("H_e_Mom","H_e_Mom",100,0,9);
		H_e_Mom.setTitle(" LUND Momentum of electrons After Cuts ");
		H_e_Mom.setTitleX("Momentum [GeV/c]");
		H_e_Mom.setFillColor(5);
		H_e_theta = new H1F("H_e_theta","H_e_theta",100,0,80);
		H_e_theta.setTitle("LUND Electron Theta angle");
		H_e_theta.setTitleX("Theta angle [degree]");
		H_e_theta.setFillColor(3);
		H_e_phi = new H1F("H_e_phi","H_e_phi",100,-180,180);
		H_e_phi.setTitle("LUND Electron Phi angle");
		H_e_phi.setTitleX("Azimuthal Angle Phi [degree]");
		H_e_phi.setFillColor(4);
		int MigrationBin=9;
		ElBinMigration = new H2F("ElBinMigration", " ElBinMigration", 9, -4.5,4.5, 9, -4.5,4.5);
		ElBinMigration.setTitle("Migration of Electrons bins");
		ElBinMigration.setTitleX(ParticleMC.GetName_1stVariable());
		ElBinMigration.setTitleY(ParticleMC.GetName_2ndVariable());
		FirstVariable = ParticleMC.GetName_1stVariable();
		SecondVariable = ParticleMC.GetName_2ndVariable();
		ThirdVariable = ParticleMC.GetName_3rdVariable();
		ForthVariable = ParticleMC.GetName_4thVariable();
		FirstBins= ParticleMC.FstBin; FirstMin=ParticleMC.Fst_Min; FirstMax=ParticleMC.Fst_Max;
		SecondBins= ParticleMC.SndBin; SecondMin=ParticleMC.Snd_Min; SecondMax=ParticleMC.Snd_Max;
		ThirdBins= ParticleMC.TrdBin; ThirdMin=ParticleMC.Trd_Min; ThirdMax= ParticleMC.Trd_Max;
		ForthBins=ParticleMC.FrtBin; ForthMin= ParticleMC.Frt_Min; ForthMax= ParticleMC.Frt_Max;
		System.out.println(" ----------------------------------");
		System.out.println(" ----- MC Muilti Bins Analysis Inizialized ---");
		System.out.println(" Bins name: " + FirstVariable + " | " + SecondVariable + " | " + ThirdVariable + " | "+ ForthVariable);
		System.out.println(" Bins are: " +  FirstBins+ " -  "+ SecondBins+ " -  " + ThirdBins+ " -  " + ForthBins);
		System.out.println(" ----------------------------------");

		
		Hadron_Lund = new MultiBins(FirstBins,SecondBins,ThirdBins,ForthBins,0);
		Hadron_Lund.SetName_1stVariable(this.FirstVariable);Hadron_Lund.SetName_2ndVariable(this.SecondVariable);
		Hadron_Lund.SetName_3rdVariable(this.ThirdVariable);Hadron_Lund.SetName_4thVariable(this.ForthVariable);
		Hadron_Lund.SetBins(this.FirstMin, this.FirstMax, this.SecondMin, this.SecondMax, this.ThirdMin, this.ThirdMax, this.ForthMin, this.ForthMax, 0, 0);
		Hadron_Lund.GenerateHistograms("Hadron_Lund");
		
		Hadron_Z = new MultiBins(FirstBins,SecondBins,100,0,0);
		Hadron_Z.SetName_1stVariable(this.FirstVariable);Hadron_Lund.SetName_2ndVariable(this.SecondVariable);
		Hadron_Z.SetName_3rdVariable(this.ThirdVariable);
		Hadron_Z.SetBins(this.FirstMin, this.FirstMax, this.SecondMin, this.SecondMax, this.ThirdMin, this.ThirdMax, 0, 0, 0, 0);
		Hadron_Z.GenerateHistograms("Hadron_Z");
		
		Hadron_PT = new MultiBins(FirstBins,SecondBins,100,0,0);
		Hadron_PT.SetName_1stVariable(this.FirstVariable);Hadron_Lund.SetName_2ndVariable(this.SecondVariable);
		Hadron_PT.SetName_3rdVariable(this.ForthVariable);
		Hadron_PT.SetBins(this.FirstMin, this.FirstMax, this.SecondMin, this.SecondMax, this.ForthMin, this.ForthMax, 0, 0, 0, 0);
		Hadron_PT.GenerateHistograms("Haron_PT");
				
		Hadron_Phi = new MultiBins(FirstBins,SecondBins,ThirdBins,ForthBins,phibins);
		Hadron_Phi.SetName_1stVariable(this.FirstVariable);Hadron_Lund.SetName_2ndVariable(this.SecondVariable);
		Hadron_Phi.SetName_3rdVariable(this.ThirdVariable); Hadron_Phi.SetName_4thVariable(this.ForthVariable); Hadron_Phi.SetName_5thVariable("Phi cm");
		Hadron_Phi.SetBins(this.FirstMin, this.FirstMax, this.SecondMin, this.SecondMax, this.ThirdMin, this.ThirdMax, this.ForthMin,this.ForthMax, 0, 360);
		Hadron_Phi.GenerateHistograms("Haron_Phi"); 
		       
		       
		MultiBins ZPt_Migration = new MultiBins(ParticleMC.FrtBin,ParticleMC.SndBin,MigrationBin,MigrationBin,0);
		ZPt_Migration.SetName_1stVariable(ParticleMC.GetName_1stVariable());ZPt_Migration.SetName_2ndVariable(ParticleMC.GetName_2ndVariable());
		ZPt_Migration.SetName_3rdVariable(ParticleMC.GetName_3rdVariable()+" Migration");ZPt_Migration.SetName_4thVariable(ParticleMC.GetName_4thVariable() + " Migration ");

	}


	/**
	 * Process another event 
	 * @param event
	 */
	public void add(Bank MCParticle, Bank bankMC, int Particle ) {
		Hadron_en.clear();Hadron_mom.clear(); Hadron_theta.clear();Hadron_phi.clear();
		 Hadron_z.clear(); Hadron_pt.clear();Hadron_phiTrento.clear();
		Particles_LUND.clear();Hadron_parent.clear();Hadron_pid.clear();

		CountsID=0;
		// TODO Auto-generated method stub
		LorentzVector Ve = new LorentzVector ( 0,0,0,0); 
		// HipoGroup MCParticle = event.getGroup("MC::Particle");
		// Read the electrons and photons 
		// bankMC.show();
		//MCParticle.show();
		int nrows = MCParticle.getRows();
		float px=-100,py=-100,pz=-100,e_vx=-100,e_vy=-100,e_vz=-100;
		double e_mom=0, e_theta=0, e_phi=0;
		double Q2=0, xB=0;
		Vector3 ParticleVector= new Vector3(0,0,0);
		//MCParticle.show();
		//bankMC.show();
		LorentzVector VGS = new LorentzVector(0,0,0,0);
		LorentzVector vecW2= new LorentzVector(0,0,0,0);
		double Pq=0;
		double Pl=0;
		int nrows2 = bankMC.getRows();
	//	MCParticle.show();
	//	System.out.println(" EVENTO MC ---------------------------------- ");
		for (int ipart =0 ; ipart<nrows; ipart++) {
			//System.out.println(MCParticle.getInt("pid",ipart));
			if(11 == MCParticle.getInt("pid",ipart)) {
				px = MCParticle.getFloat("px",ipart);
				py = MCParticle.getFloat("py",ipart);
				pz = MCParticle.getFloat("pz",ipart);
				e_vx = MCParticle.getFloat("vx",ipart);
				e_vy = MCParticle.getFloat("vy",ipart);
				e_vz = MCParticle.getFloat("vz",ipart);
				e_mom = Math.sqrt(px*px+py*py+pz*pz);
				Ve.setPxPyPzE(px,py,pz,e_mom);
				H_e_vx_vy.fill(e_vx, e_vy);;
				H_e_vz.fill(e_vz);
				H_e_Mom.fill(e_mom);
				this.e_p=e_mom;
	//			System.out.println(" Elettrone px " + px);
				//THIS IS WHAT MAKES IT CRASH: ParticleVector	
			   // ParticleVector.setXYZ(px,py,pz);
			    double myphi = Math.toDegrees(Math.atan(py/px));
			    double mytheta = Math.toDegrees(Math.acos(pz/e_mom));
			   
			    if (px<0 && py>0) myphi+=180;
			    else if (px<0 && py<0) myphi-=180;
			    this.e_phi=myphi;
			    this.e_theta=mytheta;
				//e_theta= Math.toDegrees(ParticleVector.theta());
				e_phi= Math.toDegrees(ParticleVector.phi());
			   // System.out.println("Theta from ParticleVector " + e_theta + " Phi " + e_phi + " Mine calculation is theta " + mytheta+" Phi "+myphi);
				H_e_theta.fill(e_theta);
				H_e_phi.fill(e_phi);
				break;
			}
		}
	
		VGS.add(this.VB);
		VGS.sub(Ve);
		Q2=-VGS.mass2(); 
		this.e_Q2=Q2;
		Pq = -this.VT.px()*VGS.px()-this.VT.py()*VGS.py()-this.VT.pz()*VGS.pz()+this.VT.e()*VGS.e();
		vecW2.add(this.VT);vecW2.add(this.VB);vecW2.sub(Ve);
		Pl= -this.VT.px()*this.VB.px()-this.VT.py()*this.VB.py()-this.VT.pz()*this.VB.pz()+this.VT.e()*this.VB.e();
		xB=Q2/(2*Pq);
		this.e_xB=xB;
		y_el_MC = Pq/Pl;
		XB_el_MC = xB; 
		Q2_el_MC=-VGS.mass2();
		W_el_MC= vecW2.mass();
		this.e_W=W_el_MC;
		//HipoGroup bankMC = event.getGroup("MC::Lund");
		//bankMC.show();
	

		float p_px=-10,p_py=-10,p_pz=-10, p_E=0;
		  
		    
		    LorentzVector PHM = new LorentzVector(0,0,0,0);
		    LorentzVector Gammacm = new LorentzVector(0,0,0,0);
		    LorentzVector Protoncm = new LorentzVector(0,0,0,0);
		    LorentzVector PIcm = new LorentzVector(0,0,0,0);
		    LorentzVector Vphotons = new LorentzVector(0,0,0,0);
		    // Modify to look over the whole LUND structure
		    
		
		for(int k = 0; k < nrows2; k++){
			if(bankMC.getInt("pid",k)==Particle) { 
				this.CountsID++;
				//bankMC.show();
				//System.out.println(Q2);
			
				Hadron_pid.add(bankMC.getInt("pid",k));
				p_px = bankMC.getFloat("px",k);
				p_py = bankMC.getFloat("py",k);
				p_pz = bankMC.getFloat("pz",k);
				p_E = bankMC.getFloat("energy",k);
				int parentID= bankMC.getInt("parent", k);
				PHM.setPxPyPzE(p_px,p_py,p_pz,p_E);
				//LorentzVector PHM = new LorentzVector(p_px,p_py,p_pz,p_E);
				double PPh = -this.VT.px()*PHM.px()-this.VT.py()*PHM.py()-this.VT.pz()*PHM.pz()+this.VT.e()*PHM.e();
				Vector3 VectorPh = VGS.vect(); // Virtual photon is beam-el scattered but adding the target doesent change tri momenta
				Vector3 VectorPi=PHM.vect();
				Vector3 VectorP = this.VT.vect();
				Vector3 rotationVector = new Vector3(0,0,0);
				TransMatrix rotationMatrix = new TransMatrix();
				rotationMatrix.compose(VectorPh);
				Vector3 PhCM = rotationMatrix.mult(VectorPh);
				Vector3 PiCM = rotationMatrix.mult(VectorPi);
				Vector3 PCM = rotationMatrix.mult(VectorP);
				Gammacm.setPxPyPzE(PhCM.x(), PhCM.y(), PhCM.z(), VGS.e());
				Protoncm.setPxPyPzE(PCM.x(), PCM.y(), PCM.z(), this.VT.e());
				PIcm.setPxPyPzE(PiCM.x(), PiCM.y(), PiCM.z(), PHM.e());
				Vphotons.add(Gammacm).add(Protoncm);
				Vector3 boost = Vphotons.boostVector();
				boost.negative();
				Vphotons.boost(boost);
				Gammacm.boost(boost);
				Protoncm.boost(boost);
				PIcm.boost(boost);
				double Pi0M= PHM.mass();
				double Pls = PIcm.vect().dot(Gammacm.vect())/Gammacm.vect().mag();
				double PT=PIcm.vect().mag2()-Pls*Pls;
				double Z_var= PPh/Pq;
	
				Vector3 VectorB = VB.vect();
				Vector3 VectorPhNorm = new Vector3(VectorPh.x()/VectorPh.mag(),VectorPh.y()/VectorPh.mag(),VectorPh.z()/VectorPh.mag());				
				Vector3 Angl1 = VectorPhNorm.cross(VectorB);
				Vector3 Angl2 = VectorPhNorm.cross(VectorPi);
				Angl1.setXYZ(Angl1.x()/Angl1.mag(), Angl1.y()/Angl1.mag(),Angl1.z()/Angl1.mag());
				Angl2.setXYZ(Angl2.x()/Angl2.mag(), Angl2.y()/Angl2.mag(),Angl2.z()/Angl2.mag());
				double AngoloCos=Angl1.dot(Angl2);
				double Angolo= Math.toDegrees(Math.acos(AngoloCos));
				//this.Phi_LAB = Math.toDegrees(Math.acos(Angolo/(Angl1.mag()*Angl2.mag())));
				Vector3 AngleA = VectorB.cross(VectorPi);
				double AngleAA = AngleA.dot(VectorPhNorm);
				Vector3 AngleB=VectorPhNorm.cross(VectorB);
				Vector3 AngleBB =VectorPhNorm.cross(VectorPi);
				double AngleSin= AngleAA/(AngleB.mag()*AngleBB.mag());
				double AlternativePhi= Math.toDegrees(Math.asin(AngleSin));
				//System.out.println("----------");
				//System.out.println(" Phi cos " +Angl1.dot(Angl2) + " sin " + AngleSin );
				double Phi_LAB = Angolo;
				//if(AngleSin<0) this.Phi_LAB=-this.Phi_LAB; 
				//System.out.println("Angolo LAB [degrees] " + Angolo);
				//if(AlternativePhi<0) AlternativePhi=-AlternativePhi; 
				//System.out.println(" Angolo LAB Alternativo " + AlternativePhi);
				if(Angolo>0 && AlternativePhi>0) Phi_LAB=Angolo;
				else if(Angolo<0 && AlternativePhi>0) Phi_LAB=Angolo;
				else if (Angolo<0  && AlternativePhi<0) Phi_LAB=Angolo+90;
				else if (Angolo>0 && AlternativePhi<0 ) Phi_LAB=360-Angolo;
				//System.out.println("PID " + Particle +  " My generted  Phi  is " + Phi_LAB + "  Particle Mass " + Pi0M + " The electron  momentum is (x,y,z):  " + px  +" " + py +" "+pz );
				
		    	this.InvMass.fill(Z_var, PT);	
		    	Hadron_en.add(p_E);
		    	Hadron_mom.add((float) Math.sqrt(px*px+py*py+pz*pz));
		    	Hadron_theta.add((float) Math.toDegrees(PHM.vect().theta()));
		        Hadron_phi.add((float) Math.toDegrees(PHM.vect().phi()));
		        Hadron_phiTrento.add((float) Phi_LAB);
		    	Hadron_z.add((float) Z_var);
		    	Hadron_pt.add((float) PT);
		    	Hadron_parent.add(parentID);
				//System.out.println("Bins " + Hadron_Lund.FstBin + " | "+ Hadron_Lund.SndBin +" |  " + Hadron_Lund.TrdBin + "|  "+ Hadron_Lund.FrtBin);
				//System.out.println("First Min "  + Hadron_Lund.Fst_Min + " First max "+ Hadron_Lund.Fst_Max + " Second min "+ Hadron_Lund.Snd_Min + " Second Max" + Hadron_Lund.Snd_Max);
				//System.out.println(" Third Min " + Hadron_Lund.Trd_Min + " Third_Max " + Hadron_Lund.Trd_Max + " Forth Min " + Hadron_Lund.Frt_Min + " Forth Max " + Hadron_Lund.Frt_Max);
				
				int BinQ2=Hadron_Lund.ComputeBin(1, Q2);
				int BinXB=Hadron_Lund.ComputeBin(2, xB);
				
				//System.out.println(" The bins are " + BinQ2+" with a Q2 " + Q2 + "and " + BinXB+ " with XB "+ xB);
			  
				//Hadron_Lund.GetH2F(BinQ2,BinXB).getXAxis().setTitle("Ciao");
			  // To readd 
				if (BinQ2>=0 && BinXB>=0) {
				Hadron_Lund.GetH2F(BinQ2, BinXB).fill(Z_var,PT);
				Hadron_Z.GetH1F(BinQ2, BinXB).fill(Z_var);
				Hadron_PT.GetH1F(BinQ2, BinXB).fill(PT);
				Hadron_Phi.GetH3F(BinQ2, BinXB).fill(Z_var,PT,Phi_LAB);
			    }
				Particles_LUND.add(Hadron_Lund);
			/*
				//Here i have to inizialize a new multibins and put it into the array.
				if( Z_var <= 0.1 ) this.z1++;
				else if (Z_var>0.1 && Z_var<=0.2)this.z2++;
				else if (Z_var>0.2 && Z_var<=0.3)this.z3++;
				else if (Z_var>0.3 && Z_var<=0.4)this.z4++;
				else if (Z_var>0.4 && Z_var<=0.5)this.z5++;
				else if (Z_var>0.5 && Z_var<=0.6)this.z6++;
				else if (Z_var>0.6 && Z_var<=0.7)this.z7++;
				else if (Z_var>0.7 && Z_var<=0.8)this.z8++;
				else if (Z_var>0.8 && Z_var<=0.9)this.z9++;
				else if (Z_var>0.9 && Z_var<=1)this.z10++;
				if (PT<=0.2) this.pt1++;
				else if (PT>0.2 &&PT<=0.3) this.pt2++;
				else if (PT>0.3 &&PT<=0.4) this.pt3++;
				else if (PT>0.4 &&PT<=0.6) this.pt4++;
				else if (PT>0.6 &&PT<=0.8) this.pt5++;
				else if (PT>0.8 &&PT<=1) this.pt6++;
				else if (PT>1 &&PT<=1.5) this.pt7++;
*/		
			}
		}
		
	}

	/**
	 *  Set the Parent ID for searching PID
	 * @param i
	 */
	public void setParent(int i) {
		// Check if i has been already used in the class
		this.parentID=i;
	}

	/**
	 * Count how many particles in the events
	 * @param String of the bin : z or pt
	 * @return the nr of particles with the given id in the events analyzided
	 * @throws FileNotFoundException 
	 * @throws IOException 
	 */
	public void getNr(String s) throws FileNotFoundException {
		File File_txt = new File("/Users/gangelini/Desktop/Plots/Lund.txt");
		if(File_txt.exists()) {
			System.out.println("File found" );
			@SuppressWarnings("resource")
			PrintWriter outP = new PrintWriter(File_txt);
			for(int ii=0;ii<Z_Bins; ii++) {
				for (int jj=0;jj<PT_Bins; jj++) {
					System.out.println(" BIN Z: " +ii+ " BIN PT: " +jj);
					System.out.println(InvMass.getBinContent(ii, jj));
					outP.println(" BIN Z: " +ii+ " BIN PT: " +jj+" Counts Pi0s" + InvMass.getBinContent(ii, jj) );
				}
			} 
		}
		System.out.print(" Z COUNTS" );

		// TODO Auto-generated method stub
		if (s=="z") {
			System.out.println(" --- Z ---- " );
			System.out.println("MC counts Z<=0.1 " + this.z1);
			System.out.println("MC count 0.1<Z bin <=0.2 " + this.z2);
			System.out.println("MC count 0.2<Z bin <=0.3 " + this.z3);
			System.out.println("MC count 0.3<Z bin <=0.4 " + this.z4);
			System.out.println("MC count 0.4<Z bin <=0.5 "+ this.z5);
			System.out.println("MC count 0.5<Z bin <=0.6 " + this.z6);
			System.out.println("MC count 0.6<Z bin <=0.7 "+ this.z7);
			System.out.println("MC count 0.7<Z bin <=0.8 "+this.z8);
			System.out.println("MC count 0.8<Z bin <=0.9 "+ this.z9);
			System.out.println("MC count 0.9<Z bin <=1 " + this.z10);

		}
		else if ( s=="pt") {
			System.out.println(" --- PT ---- " );
			System.out.println("MC counts PT<=0.2 " + this.pt1);
			System.out.println("MC count 0.2<PT <=0.3 " + this.pt2);
			System.out.println("MC count 0.3<PT <=0.4 " + this.pt3);
			System.out.println("MC count 0.4<PT  <=0.6 " + this.pt4);
			System.out.println("MC count 0.6<PT <=0.8 "+ this.pt5);
			System.out.println("MC count 0.8<PT  <=1 " + this.pt6);
			System.out.println("MC count 1<PT  <=1.5 "+ this.pt7);

		}

	}

	public void setTarget(LorentzVector target) {
		// TODO Auto-generated method stub
		this.VT=target;

	}

	public void setBeam(LorentzVector beam) {
		// TODO Auto-generated method stub
		this.VB=beam;
	}

	/**
	 * 
	 * @return the MC electron Q2 for the event 
	 */
	public double getQ2()
	{
		return this.Q2_el_MC;
	}

	/**
	 * 
	 * @return the MC electron y for the event 
	 */
	public double getY()
	{
		return this.y_el_MC;
	}

	/**	 
	 * @return the MC electron Q2 for the event 
	 */
	public double getW()
	{
		return this.W_el_MC;
	}		

	public double getXB() {
		return this.XB_el_MC;
	}

	public void Histo(String time, String workdirout) {
		String s1=workdirout;
		new File(s1+"/Plots").mkdir();
		String PathFolder = s1+"/Plots/"+time;
		File dir = new File(PathFolder);
		dir.mkdir();
		PathFolder=s1+"/Plots/"+time+"/LUND";
		new File (PathFolder).mkdir();


		System.out.println(" (Pi) - > I am Plotting LUND information");

		EmbeddedCanvas lund = new EmbeddedCanvas();
		lund.setSize(1600,1000);
		lund.divide(3,2);
		lund.setAxisTitleSize(20);
		lund.setAxisFontSize(20);
		lund.setTitleSize(20);
		// I plot z from  from 0.1 to 0.9 because first and last bin are useless
		lund.cd(0);lund.draw(H_e_vz);
		lund.cd(1);lund.draw(H_e_vx_vy);
		lund.cd(2);lund.draw(H_e_Mom);
		lund.cd(3);lund.draw(H_e_theta);
		lund.cd(4); lund.draw(H_e_phi);
		String strg0 = String.format("%s/Lund.png",PathFolder);
		System.out.println("Saving plots in "+ PathFolder);
		lund.save(strg0);	

		//I could save the invairiant mass somewhere.


	}

public int CountParticles()
{
	return CountsID;
}
/**
 * Return the MultiBins for the MC 
 * @param i
 * @return
 */
public MultiBins GetParticleBins (int i)
{
	return Hadron_Lund;
}

/**
 * Save Histograms 
 * @param workdirout 
 */
public void SaveHistograms_MultiBins(String time, String workdirout) {
	// Producing the H1F for the Getting counts as function of Phi 
	MultiBins PhiCounts = new MultiBins(Hadron_Phi.TrdBin,Hadron_Phi.FrtBin,Hadron_Phi.FthBin, 0 ,0);
	PhiCounts.SetBins(Hadron_Phi.Trd_Min, Hadron_Phi.Trd_Max, Hadron_Phi.Frt_Min, Hadron_Phi.Frt_Max,  Hadron_Phi.Fth_Min, Hadron_Phi.Fth_Max, 0, 0, 0, 0);
	PhiCounts.GenerateHistograms("PhiCounts");
	
	String s1=workdirout;
	
	String PathFolder = s1+"/Plots/"+time;
	File dir = new File(PathFolder);dir.mkdir();
	PathFolder=s1+"/Plots/"+time+"/MultiBins_Generated";
	new File (PathFolder).mkdir();
	for(int x=0 ; x <(Hadron_Lund.FstBin);x++) { //Q2bins
		String number = Integer.toString(x);
		 PathFolder=s1+"/Plots/"+time+"/MultiBins_Generated/"+Hadron_Lund.GetName_1stVariable()+ "_"+number;
		File dirbin = new File(PathFolder);dirbin.mkdir();
		for (int y=0 ; y<(Hadron_Lund.SndBin); y++) { //Xb bins
			
			TDirectory directory = new TDirectory();
			directory.mkdir("mc"); directory.cd("mc");
			
			String number2=Integer.toString(y);
    		PathFolder=s1+"/Plots/"+time+"/MultiBins_Generated/"+Hadron_Lund.GetName_1stVariable()+ "_"+number+"/"+Hadron_Lund.GetName_2ndVariable()+ "_"+number2;
        	File dirbin2 = new File(PathFolder);dirbin2.mkdir();
			File File_txt = new File(PathFolder+"/results.txt");
			try {
			File_txt.createNewFile();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if(File_txt.exists()) {
				System.out.println(" Output File Found" );
				PrintWriter outP= null;
				try {
					outP = new PrintWriter(File_txt);
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				Hadron_Lund.GetH2F(x, y).setTitle("Generated Hadron Z vs PT");
				Hadron_Z.GetH1F(x, y).setTitle(" Generated Hadron Z distributuion");
				Hadron_PT.GetH1F(x, y).setTitle(" Generated Hadron PT distribution");
				Hadron_Lund.GetH2F(x, y).setTitleX(" Z"); Hadron_Lund.GetH2F(x, y).setTitleX(" PT"); 
				Hadron_Z.GetH1F(x, y).setTitleX(" Z");Hadron_PT.GetH1F(x, y).setTitleX(" PT");
				directory.addDataSet(Hadron_Lund.GetH2F(x, y));
				directory.addDataSet(Hadron_PT.GetH1F(x, y));
				directory.addDataSet(Hadron_Z.GetH1F(x, y));
				
				for(int i = 0; i < (Hadron_Phi.TrdBin); i++) { //Z bins
					for(int k = 0;k < (Hadron_Phi.FrtBin); k++) {  // pt
						directory.mkdir("/Z"+i+"PT"+k+"/"); directory.cd("/Z"+i+"PT"+k+"/");
						System.out.println(" I just got ino Z and PT folder I made ");
						//PhiCounts.GetH1F(i, k).setTitle("Counts Bin Z: " + i + " Bin Pt: "+ k);
						//PhiCounts.GetH1F(i, k).setTitleX(" Phi [degree]");
						for(int ll=0; ll<Hadron_Phi.FthBin; ll++) {
							System.out.println(" Hadron Phi bin  " + ll);
							
							double binsize=(Hadron_Phi.Fth_Max-Hadron_Phi.Fth_Min)/Hadron_Phi.FthBin;
							System.out.println("Bin Size"+binsize);
							double binPhi = (ll)*(binsize)+(binsize/2)+Hadron_Phi.Frt_Min;
							//System.out.println(" The bin position " + binPhi);
							double value= Hadron_Phi.GetH3F(x, y).getBinContent(i, k, ll);	
							PhiCounts.GetH1F(i, k).fill(binPhi,value);
						}// Phi bins 
						System.out.println("I am adding the Phi Counts to my dataset");
						directory.addDataSet(PhiCounts.GetH1F(i, k));
					}// End Pt
				} /// End Z 
				//directory.addDataSet(Hadron_Phi.GetH3F(x, y));
				directory.writeFile(PathFolder+"/Histograms_SIDIS_MC.hipo"); // Maybe outside the loop? 
				EmbeddedCanvas LundC = new EmbeddedCanvas();
				LundC.setSize(1600,1000);
				LundC.divide(2,2);
				LundC.setAxisTitleSize(22);
				LundC.setAxisFontSize(22); 
				LundC.setTitleSize(22); 	    		
				LundC.cd(0); 
				LundC.draw(Hadron_Lund.GetH2F(x, y));
				LundC.cd(1);
				LundC.draw(Hadron_Z.GetH1F(x, y));
				LundC.cd(2);
				LundC.draw(Hadron_PT.GetH1F(x, y));
				//LundC.cd(3);
				//LundC.draw(Hadron_Phi.GetH1F(x, y));
			    String stringa = String.format("%s/LUND_Counts.png",PathFolder);
			    LundC.save(stringa);
	    	
			    for(int i = 0; i < (Hadron_Lund.TrdBin); i++) { //Z bins
				for(int k = 0;k < (Hadron_Lund.FrtBin); k++) {
					
					outP.println(" Bin Z  "+ i +" Bin PT  " + k   +"  -> LUND COUNTS " + Hadron_Lund.GetH2F(x, y).getBinContent(i, k) );
					System.out.println(" LUND Bin Z  "+ i +" Bin PT  " +  k   +"  -> LUND COUNTS " + Hadron_Lund.GetH2F(x, y).getBinContent(i, k)); 

				}//Z PT
				}//Pt bin
            outP.close();						    // Printing stuff: 	
			} //File exits
	} //Xb bins 
 } //Q2 bins
	
}

public  ArrayList<Float> Get_H_En() {
	return Hadron_en;
}

public  ArrayList<Float> Get_H_mom() {
	return Hadron_mom;
}

public  ArrayList<Float> Get_H_theta() {
	return Hadron_theta;
}
public  ArrayList<Float> Get_H_phi() {
	return Hadron_phi;
}

public ArrayList<Float> Get_H_phiTrento(){
	return Hadron_phiTrento;
}
public  ArrayList<Float> Get_H_z() {
	return Hadron_z;
}
public  ArrayList<Float> Get_H_PT() {
	return Hadron_pt;
}

public ArrayList<Integer> Get_H_PID(){
	return Hadron_pid;
}

public ArrayList<Integer> Get_H_Parent(){
	return Hadron_parent;
}


public void ResetHadrons () {
	
	Hadron_pt.clear();Hadron_z.clear();Hadron_phi.clear();Hadron_theta.clear();Hadron_mom.clear();Hadron_en.clear();
	Hadron_phiTrento.clear();Hadron_parent.clear();Hadron_pid.clear();
	
}

}
