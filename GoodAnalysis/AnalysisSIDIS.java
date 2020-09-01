package GoodAnalysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.H3F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.TGCanvas;
import org.jlab.jnp.physics.LorentzVector;

public class AnalysisSIDIS {
	private double pi_peak=0.135; // this quantity is used to find the center of the peak 
	private double EC_res=0.012; //sigma of the pi0 peak in the calorimeter
private double Polarization = 0.863;
	public static H1F Phi_Dist_Overall,Z_Dist_Overall,XB_Dist_Overall,Pt_Dist_Overall, BSA_Phi;
	static List<H1F> InvMassPlots = new ArrayList<H1F>();
	static List<H1F> N0_Phi = new ArrayList<H1F>();
	static List<H1F> N1_Phi = new ArrayList<H1F>();
	static List<H1F> Pt_Dist = new ArrayList<H1F>();
	static List<H1F> Z_Dist = new ArrayList<H1F>();
	static List<H1F> XB_Dist = new ArrayList<H1F>();


	private H1F H_MissingMass = new H1F("H_MissingMass",200,0,3);

	private int zbins =0;
	private int PTbins=0;
	private int PT_Bins=0;
	private int Z_Bins=0;
	private int PhiBins=10;
	private double Phi_Min=-180, Phi_Max=180;
	private MultiBins PionCounts;
	private H3F InvMass;
	private H3F N0_Count;
	private H3F N1_Count;
	private boolean SingleHadron;


	private MultiBins Hadron_Z,Hadron_PT,Pion_Energy, Pion_Theta, Pion_Momentum, Pion_Z, Pion_Phi, Pion_PT,Pion_PhiClas;
	private MultiBins Dist_Q2,Dist_xB,Dist_W,Dist_y,Dist_eps,Dist_Pt,Dist_z;

	public AnalysisSIDIS(boolean PolygonalBinClas, MultiBins MultiB,int Z_bins, double Z_min, double Z_max, int PT_bins, double PT_min, double PT_max, boolean HadronBoolean)
	{
		double []bin100En = new double[100];
		for(int i=0; i<100; i++) {
			//double binsize=0.1;
			double binsize =8.0/99 ;
			double minvalue=0.0;
			double conto=minvalue+binsize*i;
			bin100En[i]=conto;
			//System.out.println("Ao"+conto);
		}
		Pion_Energy = new MultiBins(MultiB.FstBin,MultiB.SndBin,100,0,0);
		Pion_Energy.SetUnevenBins(MultiB.GetUnevenBins(1), MultiB.GetUnevenBins(2),bin100En, new double[]{0,0},new double[]{0,0});
		Pion_Energy.GenerateHistograms("Pion_Energy");
		if(PolygonalBinClas==true) {Pion_Energy.InizializeClasPolygons(); }
		
		
		double []bin100MoM = new double[100];;
		for(int i=0; i<100; i++) {
			double binsize =10.0/99 ;
			double minvalue=0.0;
			double conto=minvalue+binsize*i;
			bin100MoM[i]=conto;
		}
		Pion_Momentum = new MultiBins(MultiB.FstBin,MultiB.SndBin,100,0,0);
		Pion_Momentum.SetUnevenBins(MultiB.GetUnevenBins(1), MultiB.GetUnevenBins(2),bin100MoM, new double[]{0,0},new double[]{0,0});
		Pion_Momentum.GenerateHistograms("Pion_Momentum");
		if(PolygonalBinClas==true) {Pion_Momentum.InizializeClasPolygons(); }
		
		double []bin100Theta= new double[100];
		for(int i=0; i<100; i++) {
			double binsize =50.0/99 ;
			double minvalue=0.0;
			double conto=minvalue+binsize*i;
			bin100Theta[i]=conto;
		}
		Pion_Theta = new MultiBins(MultiB.FstBin,MultiB.SndBin,100,0,0);
		Pion_Theta.SetUnevenBins(MultiB.GetUnevenBins(1), MultiB.GetUnevenBins(2),bin100Theta, new double[]{ 0,0},new double[]{ 0,0});
		Pion_Theta.GenerateHistograms("Pion_Theta");
		if(PolygonalBinClas==true) {Pion_Theta.InizializeClasPolygons(); }
		
		double []bin100Z= new double[100];
		for(int i=0; i<100; i++) {
			double binsize =1.0/99 ;
			double minvalue=0.0;
			double conto=minvalue+binsize*i;
			bin100Z[i]=conto;
		}
		Pion_Z = new MultiBins(MultiB.FstBin,MultiB.SndBin,100,0,0);
		Pion_Z.SetUnevenBins(MultiB.GetUnevenBins(1), MultiB.GetUnevenBins(2),bin100Z, new double[]{ 0,0},new double[]{ 0,0});
		Pion_Z.GenerateHistograms("Pion_Z");
		if(PolygonalBinClas==true) {Pion_Z.InizializeClasPolygons(); }
		
		double []bin100PT= new double[100];
		for(int i=0; i<100; i++) {
			double binsize =2.0/99 ;
			double minvalue=0.0;
			double conto=minvalue+binsize*i;
			bin100PT[i]=conto;
		}
		Pion_PT = new MultiBins(MultiB.FstBin,MultiB.SndBin,100,0,0);
		Pion_PT.SetUnevenBins(MultiB.GetUnevenBins(1), MultiB.GetUnevenBins(2),bin100PT, new double[]{ 0,0},new double[]{ 0,0});
		Pion_PT.GenerateHistograms("Pion_PT");
		if(PolygonalBinClas==true) {Pion_PT.InizializeClasPolygons(); }
		
		double []bin100Phi= new double[100];
		for(int i=0; i<100; i++) {
			double binsize =360.0/99 ;
			double minvalue=0.0;
			double conto=minvalue+binsize*i;
			bin100Phi[i]=conto;
		}
		Pion_Phi = new MultiBins(MultiB.FstBin,MultiB.SndBin,100,0,0);
		Pion_Phi.SetUnevenBins(MultiB.GetUnevenBins(1), MultiB.GetUnevenBins(2),bin100Phi, new double[]{0,0},new double[]{ 0,0});
		Pion_Phi.GenerateHistograms("Pion_Phi");
		if(PolygonalBinClas==true) {Pion_Phi.InizializeClasPolygons(); }
		
		double []bin100PhiClas= new double[100];
		for(int i=0; i<100; i++) {
			double binsize =360.0/99 ;
			double minvalue=-180.0;
			double conto=minvalue+binsize*i;
			bin100PhiClas[i]=conto;
		}
		Pion_PhiClas =  new MultiBins(MultiB.FstBin,MultiB.SndBin,100,0,0);
		Pion_PhiClas.SetUnevenBins(MultiB.GetUnevenBins(1), MultiB.GetUnevenBins(2),bin100PhiClas, new double[]{ 0,0},new double[]{ 0,0});
		Pion_PhiClas.GenerateHistograms("Pion_PhiClas");
		if(PolygonalBinClas==true) {Pion_PhiClas.InizializeClasPolygons(); }
		
		
		double []binDistQ2= new double[400];
		for(int i=0; i<400; i++) {
			double binsize =12.0/399 ;
			double minvalue=1.0;
			double conto=minvalue+binsize*i;
			binDistQ2[i]=conto;
		}
		Dist_Q2 = new MultiBins(MultiB.FstBin,MultiB.SndBin,400,0,0);
		Dist_Q2.SetUnevenBins(MultiB.GetUnevenBins(1), MultiB.GetUnevenBins(2),binDistQ2, new double[]{0,0},new double[]{ 0,0});
		Dist_Q2.GenerateHistograms("Dist_Q2");
		if(PolygonalBinClas==true) {Dist_Q2.InizializeClasPolygons(); }
		
		double []binDistxB= new double[400];
		for(int i=0; i<400; i++) {
			double binsize = 0.9/399;
			double minvalue=0;
			double conto=minvalue+binsize*i;
			binDistxB[i]=conto;
		}
		Dist_xB = new MultiBins(MultiB.FstBin,MultiB.SndBin,400,0,0);
		Dist_xB.SetUnevenBins(MultiB.GetUnevenBins(1), MultiB.GetUnevenBins(2),binDistxB, new double[]{ 0,0},new double[]{ 0,0});
		Dist_xB.GenerateHistograms("Dist_xB");
		if(PolygonalBinClas==true) {Dist_xB.InizializeClasPolygons(); }
		
		
		double []binDistW= new double[400];
		for(int i=0; i<400; i++) {
			double binsize =3.5/399 ;
			double minvalue=1.5;
			double conto=minvalue+binsize*i;
			binDistW[i]=conto;
		}
		Dist_W = new MultiBins(MultiB.FstBin,MultiB.SndBin,400,0,0);
		Dist_W.SetUnevenBins(MultiB.GetUnevenBins(1), MultiB.GetUnevenBins(2),binDistW, new double[]{ 0,0},new double[]{ 0,0});
		Dist_W.GenerateHistograms("Dist_W");
		if(PolygonalBinClas==true) {Dist_W.InizializeClasPolygons(); }
		

		double []binDisty= new double[400];
		for(int i=0; i<400; i++) {
			double binsize =1.0/399 ;
			double minvalue=0;
			double conto=minvalue+binsize*i;
			binDisty[i]=conto;
		}
		Dist_y = new MultiBins(MultiB.FstBin,MultiB.SndBin,400,0,0);
		Dist_y.SetUnevenBins(MultiB.GetUnevenBins(1), MultiB.GetUnevenBins(2),binDisty, new double[]{ 0,0},new double[]{ 0,0});
		Dist_y.GenerateHistograms("Dist_y");
		if(PolygonalBinClas==true) {Dist_y.InizializeClasPolygons(); }
		
		
		double []binDisteps= new double[400];
		for(int i=0; i<400; i++) {
			double binsize =2.0/399 ;
			double minvalue=0;
			double conto=minvalue+binsize*i;
			binDisteps[i]=conto;
		}
		Dist_eps = new MultiBins(MultiB.FstBin,MultiB.SndBin,400,0,0);
		Dist_eps.SetUnevenBins(MultiB.GetUnevenBins(1), MultiB.GetUnevenBins(2),binDisteps, new double[]{ 0,0},new double[]{ 0,0});
		Dist_eps.GenerateHistograms("Dist_eps");
		if(PolygonalBinClas==true) {Dist_eps.InizializeClasPolygons(); }
		
		double []binDistz= new double[400];
		for(int i=0; i<400; i++) {
			double binsize =1.0/399 ;
			double minvalue=0;
			double conto=minvalue+binsize*i;
			binDistz[i]=conto;
		}
	//	System.out.println("Valore asse x z bin 0  "+ MultiB.GetUnevenBins(3)[0] );
	//	System.out.println("Valore asse x z bin 1  "+ MultiB.GetUnevenBins(3)[1] );
		//System.out.println("Valore asse y pt bin 0  "+ MultiB.GetUnevenBins(4)[0] );
	//	System.out.println("Valore asse pt y bin 1  "+ MultiB.GetUnevenBins(4)[1] );
		Dist_z = new MultiBins(MultiB.TrdBin,MultiB.FrtBin,400,0,0);
		Dist_z.SetUnevenBins(MultiB.GetUnevenBins(3), MultiB.GetUnevenBins(4),binDistz, new double[]{ 0,0},new double[]{ 0,0});
		Dist_z.GenerateHistograms("Dist_z");
		
		double []binDistPt= new double[400];
		for(int i=0; i<400; i++) {
			double binsize =2.0/399 ;
			double minvalue=0;
			double conto=minvalue+binsize*i;
			binDistPt[i]=conto;
		}
	//	System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" );
		//System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" );
	//	System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" );
	//	System.out.println("----->> Valore asse pt y bin 1  "+ MultiB.GetUnevenBins(4)[1] );
	//	System.out.println("---->> Valore asse y pt bin 0  "+ MultiB.GetUnevenBins(4)[0] );
	//	System.out.println("----->> Valore asse pt y bin 8  "+ MultiB.GetUnevenBins(4)[8] );
	//	System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" );
	//	System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" );
	//	System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" );
		Dist_Pt = new MultiBins(MultiB.TrdBin,MultiB.FrtBin,400,0,0);
		Dist_Pt.SetUnevenBins(MultiB.GetUnevenBins(3), MultiB.GetUnevenBins(4),binDistPt, new double[]{0,0},new double[]{0,0});
		Dist_Pt.GenerateHistograms("Dist_Pt");
		
		
		
		
		this.SingleHadron=HadronBoolean;
		this.Z_Bins=Z_bins;
		this.PT_Bins=PT_bins;
		//InvMass = new H3F(Z_bins, Z_min, Z_max, PT_bins, PT_min, PT_max,Mass_Bins,min_pion,max_pion);
		N0_Count = new H3F(Z_bins, Z_min,Z_max, PT_bins, PT_min, PT_max,12,0,360);
		N1_Count = new H3F(Z_bins, Z_min, Z_max, PT_bins, PT_min, PT_max,12,0,360);	
		//	Plots();
		Z_Dist_Overall = new H1F("Z_Dist_Overall","Z_Dist_Overall",100,0,1);
		Z_Dist_Overall.setTitle("Z distribution integrated");
		Z_Dist_Overall.setTitleX("Z");

		XB_Dist_Overall = new H1F("Xb_Dist_Overall","Xb_Dist_Overall",100,0,1);
		XB_Dist_Overall.setTitle("Xb distribution integrated in Z");
		XB_Dist_Overall.setTitleX("Xb");

		Pt_Dist_Overall = new H1F("Pt_Dist_Overall","Pt_Dist_Overall",100,0,2);
		Pt_Dist_Overall.setTitle("Pt distribution integrated in Z");
		Pt_Dist_Overall.setTitleX("Pt ");

		Phi_Dist_Overall = new H1F("Phi_Dist_Overall","Phi_Dist_Overall",12,0,360);
		Phi_Dist_Overall.setTitle("Phi distribution integrated in Z");
		Phi_Dist_Overall.setTitleX("Phi");

		H_MissingMass.setTitle("H_MissingMass");
		H_MissingMass.setTitleX("Missing Mass [GeV]");

	}

	public void Analyze(boolean PolygonalBinClas, LorentzVector beam, LorentzVector target, ParticleREC electronRec, int helicity, H1F xMass, MultiBins missingM, int sig, double xB, double Q2, Pi0_Particle new_Pi0s, MultiBins invariantMassBins, MultiBins H0_Counts, MultiBins H1_Counts, double W, double y, double eps  ) {
		// TODO Auto-generated method stub
		//InvMass =invariantMassBins;
		// Counts of Helicity 0 
		/*
		MultiBins H0_Counts = new MultiBins(invariantMassBins.FstBin,invariantMassBins.SndBin,invariantMassBins.TrdBin,invariantMassBins.FrtBin,this.PhiBins);
		H0_Counts.SetBins(invariantMassBins.Fst_Min, invariantMassBins.Fst_Max, invariantMassBins.Snd_Min, invariantMassBins.Snd_Max, invariantMassBins.Trd_Min, invariantMassBins.Trd_Max, invariantMassBins.Frt_Min, invariantMassBins.Frt_Max, this.Phi_Min, this.Phi_Max);
		H0_Counts.SetName_1stVariable(invariantMassBins.GetName_1stVariable()); H0_Counts.SetName_2ndVariable(invariantMassBins.GetName_2ndVariable()); H0_Counts.SetName_3rdVariable(invariantMassBins.GetName_3rdVariable());
		H0_Counts.SetName_4thVariable(invariantMassBins.GetName_4thVariable()); H0_Counts.SetName_5thVariable(" Phi [Degrees] ");
		H0_Counts.GenerateHistograms();

		//Counts of Helicity 1 
		MultiBins H1_Counts = new MultiBins(invariantMassBins.FstBin,invariantMassBins.SndBin,invariantMassBins.TrdBin,invariantMassBins.FrtBin,0);
		H1_Counts.SetBins(invariantMassBins.Fst_Min, invariantMassBins.Fst_Max, invariantMassBins.Snd_Min, invariantMassBins.Snd_Max, invariantMassBins.Trd_Min, invariantMassBins.Trd_Max, invariantMassBins.Frt_Min, invariantMassBins.Frt_Max, this.Phi_Min, this.Phi_Max);
		H1_Counts.SetName_1stVariable(invariantMassBins.GetName_1stVariable()); H1_Counts.SetName_2ndVariable(invariantMassBins.GetName_2ndVariable()); H1_Counts.SetName_3rdVariable(invariantMassBins.GetName_3rdVariable());
		H1_Counts.SetName_4thVariable(invariantMassBins.GetName_4thVariable()); H1_Counts.SetName_5thVariable(" Phi [Degrees] ");
		H1_Counts.GenerateHistograms();
		 */

		//Electron_counts.fill( Q2, xB);
				String FirstVariable = invariantMassBins.GetName_1stVariable();
				String SecondVariable = invariantMassBins.GetName_2ndVariable();
				String ThirdVariable = invariantMassBins.GetName_3rdVariable();
				String ForthVariable = invariantMassBins.GetName_4thVariable();
				int FirstBins= invariantMassBins.FstBin; double FirstMin=invariantMassBins.Fst_Min; double FirstMax=invariantMassBins.Fst_Max;
				int SecondBins= invariantMassBins.SndBin; double SecondMin=invariantMassBins.Snd_Min; double SecondMax=invariantMassBins.Snd_Max;
				int ThirdBins= invariantMassBins.TrdBin; double ThirdMin=invariantMassBins.Trd_Min; double ThirdMax= invariantMassBins.Trd_Max;
				int ForthBins=invariantMassBins.FrtBin; double ForthMin= invariantMassBins.Frt_Min; double ForthMax= invariantMassBins.Frt_Max;

				Hadron_Z = new MultiBins(FirstBins,SecondBins,100,0,0);
				Hadron_Z.SetName_1stVariable(FirstVariable);
				Hadron_Z.SetName_3rdVariable(ThirdVariable);
				Hadron_Z.SetBins(FirstMin, FirstMax, SecondMin, SecondMax, ThirdMin, ThirdMax, 0, 0, 0, 0);
				Hadron_Z.GenerateHistograms("Hadron_Z");
		        Hadron_Z.InizializeClasPolygons();
				Hadron_PT = new MultiBins(FirstBins,SecondBins,100,0,0);
				Hadron_PT.SetName_1stVariable(FirstVariable);
				Hadron_PT.SetName_3rdVariable(ForthVariable);
				Hadron_PT.SetBins(FirstMin, FirstMax,SecondMin, SecondMax, ForthMin, ForthMax, 0, 0, 0, 0);
				Hadron_PT.GenerateHistograms("Haron_PT"); 
				Hadron_PT.InizializeClasPolygons();

		if(this.SingleHadron==true) {
			int k= new_Pi0s.GetIndex_MostEnergy();

			Z_Dist_Overall.fill(new_Pi0s.getZ().get(k));
			XB_Dist_Overall.fill(xB);
			Pt_Dist_Overall.fill(new_Pi0s.getPt().get(k));

			
			LorentzVector MissM= new LorentzVector();
			LorentzVector electronLV = new LorentzVector();
			electronLV.setPxPyPzE(electronRec.px(), electronRec.py(),electronRec.pz(), electronRec.e());
			MissM.add(beam).add(target).sub(electronLV).sub(new_Pi0s.getLorentzVector().get(k));
			double MissingMass= MissM.mass();
			xMass.fill(MissingMass);
			int Q2bin,XBbin;
			if(PolygonalBinClas==true) {
				 Q2bin= PionCounts.ComputePolygonalBins(Q2, xB).get(0);
				 XBbin = PionCounts.ComputePolygonalBins(Q2, xB).get(1);
				}
				else {
				 Q2bin = PionCounts.ComputeBin(1, Q2);
				 XBbin = PionCounts.ComputeBin(2, xB);
				}
			
			

			missingM.GetH3F(Q2bin, XBbin).fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), MissingMass);
	
				if(Q2bin>=0 && XBbin>=0 ) {
					invariantMassBins.GetH3F(Q2bin,XBbin).fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getMass().get(k));
					if(helicity==-1)	{
						N0_Count.fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getPhi().get(k)); 
						if (new_Pi0s.getMass().get(k) <(pi_peak+sig*EC_res) && new_Pi0s.getMass().get(k)> (pi_peak-sig*EC_res)) {
							H0_Counts.GetH3F(Q2bin, XBbin).fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getPhi().get(k));
							
						}
					}
					else if(helicity==1) 
					{
						N1_Count.fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getPhi().get(k));
						if (new_Pi0s.getMass().get(k) <(pi_peak+sig*EC_res) && new_Pi0s.getMass().get(k)> (pi_peak-sig*EC_res)) {
							H1_Counts.GetH3F(Q2bin, XBbin).fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getPhi().get(k));
						}
					}
				}
		}
		else if (this.SingleHadron==false){		
			int Q2bin,XBbin;
			if(PolygonalBinClas==true) {
				 Q2bin= PionCounts.ComputePolygonalBins(Q2, xB).get(0);
				 XBbin = PionCounts.ComputePolygonalBins(Q2, xB).get(1);
				}
				else {
				 Q2bin = PionCounts.ComputeBin(1, Q2);
				 XBbin = PionCounts.ComputeBin(2, xB);
				}
			
			Dist_Q2.GetH1F(Q2bin, XBbin).fill(Q2);
			Dist_xB.GetH1F(Q2bin, XBbin).fill(xB);
			Dist_W.GetH1F(Q2bin, XBbin).fill(W);
			Dist_y.GetH1F(Q2bin, XBbin).fill(y);
			Dist_eps.GetH1F(Q2bin, XBbin).fill(eps);
			
			for(int k=0;k<new_Pi0s.getNrPi0s(); k++) {
				Z_Dist_Overall.fill(new_Pi0s.getZ().get(k));
				XB_Dist_Overall.fill(xB);
				Pt_Dist_Overall.fill(new_Pi0s.getPt().get(k));
			
				LorentzVector MissM= new LorentzVector();
				LorentzVector electronLV = new LorentzVector();
				electronLV.setPxPyPzE(electronRec.px(), electronRec.py(),electronRec.pz(), electronRec.e());
				MissM.add(beam).add(target).sub(electronLV).sub(new_Pi0s.getLorentzVector().get(k));
				double MissingMass= MissM.mass();
				xMass.fill(MissingMass);
				
			
				missingM.GetH3F(Q2bin, XBbin).fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), MissingMass);

				
					if(Q2bin>=0 && XBbin>=0 ) {
						invariantMassBins.GetH3F(Q2bin,XBbin).fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getMass().get(k));
						
						if(helicity==-1)	{
							N0_Count.fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getPhi().get(k)); 
							if (new_Pi0s.getMass().get(k) <(pi_peak+sig*EC_res) && new_Pi0s.getMass().get(k)> (pi_peak-sig*EC_res)) {
								H0_Counts.GetH3F(Q2bin, XBbin).fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getPhi().get(k));
								
							}
						}

						else if(helicity==1) 
						{
							
							N1_Count.fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getPhi().get(k));
							if (new_Pi0s.getMass().get(k) <(pi_peak+sig*EC_res) && new_Pi0s.getMass().get(k)> (pi_peak-sig*EC_res)) {
								H1_Counts.GetH3F(Q2bin, XBbin).fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getPhi().get(k));
								
							}
						}
					}
					
			}


		}
	}

	public void AnalyzeCharged(boolean PolygonalBinClas, LorentzVector beam, LorentzVector target, ParticleREC electronRec, int helicity, H1F xMass, MultiBins missingM, int sig, double xB, double Q2, ChargedPi_Particle new_Pi0s, MultiBins PionCounts, MultiBins Counts_Phi, MultiBins H0_Counts, MultiBins H1_Counts , MultiBins Electron_counts, double W, double y, double eps) {
		Z_Dist_Overall.fill(new_Pi0s.getZ());
		XB_Dist_Overall.fill(xB);
		Pt_Dist_Overall.fill(new_Pi0s.getPt());			
		
	int Q2bin,XBbin;
if(PolygonalBinClas==true) {
		 Q2bin= PionCounts.ComputePolygonalBins(Q2, xB).get(0);
		 XBbin = PionCounts.ComputePolygonalBins(Q2, xB).get(1);
		}
		else {
		 Q2bin = PionCounts.ComputeBin(1, Q2);
		 XBbin = PionCounts.ComputeBin(2, xB);
		 //System.out.println("My Q2"+ Q2 + " my Q2 bin "+ Q2bin + " my xB" + xB + " my xB bin "+XBbin );
		}
		

if(Q2bin>=0 && XBbin>=0) {
		Dist_Q2.GetH1F(Q2bin, XBbin).fill(Q2);
		Dist_xB.GetH1F(Q2bin, XBbin).fill(xB);
		Dist_W.GetH1F(Q2bin, XBbin).fill(W);
		Dist_y.GetH1F(Q2bin, XBbin).fill(y);
		Dist_eps.GetH1F(Q2bin, XBbin).fill(eps);

		
		int Zbin = PionCounts.ComputeBin(3, new_Pi0s.getZ());
		int Ptbin= PionCounts.ComputeBin(4, new_Pi0s.getPt());

		if(Zbin>=0 && Ptbin >=0) {
		Dist_z.GetH1F(Zbin, Ptbin).fill(new_Pi0s.getZ());
		Dist_Pt.GetH1F(Zbin, Ptbin).fill(new_Pi0s.getPt());
		}
		xMass.fill(new_Pi0s.getMissMass());
		missingM.GetH3F(Q2bin, XBbin).fill(new_Pi0s.getZ(), new_Pi0s.getPt(), new_Pi0s.getMissMass());
	
			Pion_Energy.GetH1F(Q2bin, XBbin).fill(new_Pi0s.getEnergy());
			Pion_Theta.GetH1F(Q2bin, XBbin).fill(new_Pi0s.getTheta());
			Pion_Momentum.GetH1F(Q2bin, XBbin).fill(new_Pi0s.getMomentum());
			Pion_Z.GetH1F(Q2bin, XBbin).fill(new_Pi0s.getZ());
			Pion_Phi.GetH1F(Q2bin, XBbin).fill(new_Pi0s.getPhi());
			Pion_PT.GetH1F(Q2bin, XBbin).fill(new_Pi0s.getPt());
			Pion_PhiClas.GetH1F(Q2bin, XBbin).fill(new_Pi0s.getPhiClas12());
			
			
	
				PionCounts.GetH2F(Q2bin,XBbin).fill(new_Pi0s.getZ(), new_Pi0s.getPt());
				Counts_Phi.GetH3F(Q2bin, XBbin).fill(new_Pi0s.getZ(), new_Pi0s.getPt(), new_Pi0s.getPhi());
				
				if(helicity==-1)		{
					N0_Count.fill(new_Pi0s.getZ(), new_Pi0s.getPt(), new_Pi0s.getPhi()); 
					H0_Counts.GetH3F(Q2bin, XBbin).fill(new_Pi0s.getZ(), new_Pi0s.getPt(), new_Pi0s.getPhi());
				}
				else if(helicity==1) 
				{
					N1_Count.fill(new_Pi0s.getZ(), new_Pi0s.getPt(), new_Pi0s.getPhi());
					H1_Counts.GetH3F(Q2bin, XBbin).fill(new_Pi0s.getZ(), new_Pi0s.getPt(), new_Pi0s.getPhi());
				}
			

			} //Q2 and Xb Bins >0
		
	}




	public void Plots()
	{


		BSA_Phi = new H1F("BSA_Phi","BSA_Phi",12,0,360);
		BSA_Phi.setTitle("BSA Phi Dist");
		BSA_Phi.setTitleX("Phi [rad]");

		Z_Dist_Overall = new H1F("Z_Dist_Overall","Z_Dist_Overall",100,0,1);
		Z_Dist_Overall.setTitle("Z distribution integrated");
		Z_Dist_Overall.setTitleX("Z");

		XB_Dist_Overall = new H1F("Xb_Dist_Overall","Xb_Dist_Overall",100,0,1);
		XB_Dist_Overall.setTitle("Xb distribution integrated in Z");
		XB_Dist_Overall.setTitleX("Xb");

		Pt_Dist_Overall = new H1F("Pt_Dist_Overall","Pt_Dist_Overall",100,0,2);
		Pt_Dist_Overall.setTitle("Pt^2 distribution integrated in Z");
		Pt_Dist_Overall.setTitleX("Pt^2 ");

		Phi_Dist_Overall = new H1F("Phi_Dist_Overall","Phi_Dist_Overall",12,0,360);
		Phi_Dist_Overall.setTitle("Phi distribution integrated in Z");
		Phi_Dist_Overall.setTitleX("Phi");


	}


	public void CountPi0s( int fit, int sig, H1F xMass, MultiBins MissingM, MultiBins invariantMassBins, MultiBins H0_Counts, MultiBins H1_Counts, MultiBins Electron_counts, String time, boolean mC, String workdirout) throws FileNotFoundException {
		System.out.println("*************************");
		System.out.println(" Bin content" + Electron_counts.GetH2().getBinContent(0, 0));
		System.out.println(" Bin content +1 " + Electron_counts.GetH2().getBinContent(1, 1));

		System.out.println(" Bin data" + Electron_counts.GetH2().getData(0, 0));
		System.out.println(" Bin data +1 " + Electron_counts.GetH2().getData(1, 1));


		List<F1D> functions = new ArrayList<F1D>(); 

		//System.out.println(" Valore Histogramma " + this.InvMass.getSlicesZ().get(20).projectionX().getEntries());
		double totEntry=0;

		// Producing the H2F containing the pion counts 
		MultiBins GaussianCounts = new MultiBins(invariantMassBins.FstBin,invariantMassBins.SndBin,invariantMassBins.TrdBin,invariantMassBins.FrtBin,0);
		GaussianCounts.SetBins(invariantMassBins.Fst_Min, invariantMassBins.Fst_Max, invariantMassBins.Snd_Min, invariantMassBins.Snd_Max, invariantMassBins.Trd_Min, invariantMassBins.Trd_Max, invariantMassBins.Frt_Min, invariantMassBins.Frt_Max, 0, 0);
		GaussianCounts.GenerateHistograms("GaussianCounts");


		MultiBins BinCounts = new MultiBins(invariantMassBins.FstBin,invariantMassBins.SndBin,invariantMassBins.TrdBin,invariantMassBins.FrtBin,0);
		BinCounts.SetBins(invariantMassBins.Fst_Min, invariantMassBins.Fst_Max, invariantMassBins.Snd_Min, invariantMassBins.Snd_Max, invariantMassBins.Trd_Min, invariantMassBins.Trd_Max, invariantMassBins.Frt_Min, invariantMassBins.Frt_Max, 0, 0);
		BinCounts.GenerateHistograms("BinCounts");




		int Mass_Bins=invariantMassBins.FthBin;
		double min_pion=invariantMassBins.Fth_Min;
		double max_pion= invariantMassBins.Fth_Max;


		int Phi_Bins= H0_Counts.FthBin;
		double min_phi = 0.0;
		double max_phi = 360.0;

		System.out.println(" Phi_Bins " + Phi_Bins);
		System.out.println(" Min Phi " + min_phi + " Max Phi "+ max_phi);

	


		String s1=workdirout;
		new File(s1+"/Plots").mkdir();
		TDirectory directory = new TDirectory();
		TDirectory directoryMain = new TDirectory();	
		String PathFolder = s1+"/Plots/"+time;
		File dir = new File(PathFolder);dir.mkdir();
		PathFolder=s1+"/Plots/"+time+"/MultiBins";
		new File (PathFolder).mkdir();
		System.out.println(" (Pi) - > I am Plotting Multi Bins Histograms");
		// Looping over all bins and fill the 3Dimensional Invariant Mass bins (z,pt,mass) by using the 5dimensional information of InvariantMassBins
		for(int x=0 ; x <(invariantMassBins.FstBin);x++) { //Q2bins
			String number = Integer.toString(x);
			PathFolder=s1+"/Plots/"+time+"/MultiBins/"+invariantMassBins.GetName_1stVariable()+ "_"+number;
			File dirbin = new File(PathFolder);dirbin.mkdir();
			//	directory.mkdir(PathFolder); directory.cd("PathFolder");
			for (int y=0 ; y<(invariantMassBins.SndBin); y++) { //Xb bins

				// Producing the H1F for the invariant mass computation 
				MultiBins PionMass = new MultiBins(invariantMassBins.TrdBin,invariantMassBins.FrtBin,invariantMassBins.FthBin, 0 ,0);
				PionMass.SetBins(invariantMassBins.Trd_Min, invariantMassBins.Trd_Max, invariantMassBins.Frt_Min, invariantMassBins.Frt_Max, invariantMassBins.Fth_Min, invariantMassBins.Fth_Max, 0, 0, 0, 0);
				PionMass.GenerateHistograms("PionMass");

				// Producing the H1F for the BSA mass computation 
				MultiBins PhiBSA = new MultiBins(invariantMassBins.TrdBin,invariantMassBins.FrtBin,H0_Counts.FthBin, 0 ,0);
				PhiBSA.SetBins(invariantMassBins.Trd_Min, invariantMassBins.Trd_Max, invariantMassBins.Frt_Min, invariantMassBins.Frt_Max, H0_Counts.Fth_Min, H0_Counts.Fth_Max, 0, 0, 0, 0);
				PhiBSA.GenerateHistograms("PhiBSA");

				// here I have to count the stuff:
				// Producing the H1F for the Helicity mass computation 
				MultiBins HelicityPositive = new MultiBins(invariantMassBins.TrdBin,invariantMassBins.FrtBin,H0_Counts.FthBin, 0 ,0);
				HelicityPositive.SetBins(invariantMassBins.Trd_Min, invariantMassBins.Trd_Max, invariantMassBins.Frt_Min, invariantMassBins.Frt_Max, H0_Counts.Fth_Min, H0_Counts.Fth_Max, 0, 0, 0, 0);
				HelicityPositive.GenerateHistograms("HelicityPositive");

				MultiBins HelicityNegative = new MultiBins(invariantMassBins.TrdBin,invariantMassBins.FrtBin,H0_Counts.FthBin, 0 ,0);
				HelicityNegative.SetBins(invariantMassBins.Trd_Min, invariantMassBins.Trd_Max, invariantMassBins.Frt_Min, invariantMassBins.Frt_Max, H0_Counts.Fth_Min, H0_Counts.Fth_Max, 0, 0, 0, 0);
				HelicityNegative.GenerateHistograms("HelicityNegative");


				MultiBins MissingMassX = new MultiBins(MissingM.TrdBin,MissingM.FrtBin,MissingM.FthBin,0,0);
				MissingMassX.SetBins(MissingM.Trd_Min, MissingM.Trd_Max, MissingM.Frt_Min, MissingM.Frt_Max, MissingM.Fth_Min, MissingM.Fth_Max, 0, 0, 0, 0);
				MissingMassX.GenerateHistograms("MissingMass");

				String number2=Integer.toString(y);
				PathFolder=s1+"/Plots/"+time+"/MultiBins/"+invariantMassBins.GetName_1stVariable()+ "_"+number+"/"+invariantMassBins.GetName_2ndVariable()+ "_"+number2;
				File dirbin2 = new File(PathFolder);dirbin2.mkdir();
				File File_txt = new File(PathFolder+"/results.txt");
				directoryMain.mkdir("/main/"); directoryMain.cd("/main/");

				//directoryMain.addDataSet(GaussianCounts.GetH2F(, k));

				try {
					File_txt.createNewFile();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(File_txt.exists()) {
					System.out.println(" Output File Found" );
					PrintWriter outP = new PrintWriter(File_txt);
					EmbeddedCanvas photons = new EmbeddedCanvas();
					EmbeddedCanvas BSA = new EmbeddedCanvas();
					EmbeddedCanvas CountsH1 = new EmbeddedCanvas();
					photons.setSize(1600,1000);BSA.setSize(1600,1000);CountsH1.setSize(1600,1000);
					photons.divide(invariantMassBins.TrdBin,invariantMassBins.FrtBin);
					BSA.divide(invariantMassBins.TrdBin,invariantMassBins.FrtBin);
					CountsH1.divide(invariantMassBins.TrdBin,invariantMassBins.FrtBin);
					photons.setAxisTitleSize(20); BSA.setAxisTitleSize(20); CountsH1.setAxisTitleSize(20);
					photons.setAxisFontSize(20); BSA.setAxisFontSize(20);CountsH1.setAxisTitleSize(20);
					photons.setTitleSize(20); BSA.setTitleSize(20);CountsH1.setAxisTitleSize(20);	
					int myindex=0;
					int badZ = 0;
					int badPT= 0; 
					boolean badFit=false;
					for(int i = 0; i < (invariantMassBins.TrdBin); i++) { //Z bins
						for(int k = 0;k < (invariantMassBins.FrtBin); k++) { // PT bins 
							directory.mkdir("/Z"+i+"PT"+k+"/"); directory.cd("/Z"+i+"PT"+k+"/");
							PionMass.GetH1F(i, k).setTitle("InvMass Bin Z:" +i+ " Bin Pt:" + k);
							PionMass.GetH1F(i, k).setTitleX("Invariant Mass [GeV]");
							PionMass.GetH1F(i, k).setTitleY(" Counts ");
							PhiBSA.GetH1F(i, k).setTitle("BSA Bin Z: " + i + " Bin Pt: "+ k);
							PhiBSA.GetH1F(i, k).setTitleX(" Phi [degree]");

							for(int ii=0; ii < MissingM.FthBin; ii++) {		
								MissingMassX.GetH1F(i, k).fill(MissingM.GetH3F(x, y).getBinContent(i, k, ii));
								System.out.println("Missing mass Value "+MissingM.GetH3F(x, y).getBinContent(i, k, ii) );
							}

							for(int ll=0; ll< Mass_Bins ; ll++) {
								//System.out.println(" I am making histo photons" );
								double binMass=(ll)*((max_pion-min_pion)/Mass_Bins)+(((max_pion-min_pion)/Mass_Bins)/2)+min_pion;
								PionMass.GetH1F(i, k).fill(binMass,invariantMassBins.GetH3F(x, y).getBinContent(i, k, ll));
								totEntry+=invariantMassBins.GetH3F(x, y).getBinContent(i, k, ll);
							}
							if (mC==false) {
								for(int ll=0; ll<Phi_Bins; ll++) {
									//	System.out.println(" Min Bin " + min_phi);
									System.out.println("Bin Size"+(max_phi-min_phi)/Phi_Bins);
									double binsize=(max_phi-min_phi)/Phi_Bins;
									double binPhi = (ll)*(binsize)+(binsize/2)+min_phi;
									//System.out.println(" The bin position " + binPhi);
									System.out.println("The bin content for H0 is "+ H0_Counts.GetH3F(x, y).getBinContent(i, k, ll));
									System.out.println("The bin content for H1 is "+ H1_Counts.GetH3F(x, y).getBinContent(i, k, ll));
									double num= H1_Counts.GetH3F(x, y).getBinContent(i, k, ll) - H0_Counts.GetH3F(x, y).getBinContent(i, k, ll);
									double den= H0_Counts.GetH3F(x, y).getBinContent(i, k, ll) + H1_Counts.GetH3F(x, y).getBinContent(i, k, ll);
									double Bsa= (1/Polarization)*(num/den);
									System.out.println( "NUM value" + num +"DEN value" + den+ " for phi bin "+ ll + " BSA " +Bsa);
									PhiBSA.GetH1F(i, k).fill(binPhi,Bsa);
									System.out.println( "H0 Counts are " + H0_Counts.GetH3F(x, y).getBinContent(i, k, ll));
									System.out.println( "H1 Counts are " + H1_Counts.GetH3F(x, y).getBinContent(i, k, ll));
									System.out.println( "BSA Phi are " + 	PhiBSA.GetH1F(i, k).getBinContent(ll));
									HelicityNegative.GetH1F(i, k).fill(binPhi,H0_Counts.GetH3F(x, y).getBinContent(i, k, ll));
									HelicityPositive.GetH1F(i, k).fill(binPhi,H1_Counts.GetH3F(x, y).getBinContent(i, k, ll));
								}
							}
							double max_pion_fit=0.24;
							F1D fun = null ;
							F1D fun2= null;
							DataFitter fitto = new DataFitter();
							if (fit ==1)  { 
								fun = new F1D("f1d","[p0]+[p1]*x +[p2]*x*x+[p3]*x*x*x+[a]*gaus(x,[m],[s])",min_pion,max_pion_fit);
								fun.setParameter(0, totEntry/Mass_Bins);
								fun.setParameter(1, 10);
								fun.setParameter(2, 10);
								fun.setParameter(3, -10); // Try negative next time 
								fun.setParameter(4, totEntry/Mass_Bins);
								fun.setParameter(5, 0.135);
								fun.setParameter(6, 0.013);
								fun.setLineColor(2);
								fun.setLineWidth(5);
								DataFitter.fit(fun, PionMass.GetH1F(i, k), "Q");
								
								System.out.println(" FITTING RESULTS ");
								System.out.println(" Par 0 " +fun.getParameter(0));
								System.out.println(" Par 1 " +fun.getParameter(1));
								System.out.println(" Par 2 " +fun.getParameter(2));
								System.out.println(" Par 3 " +fun.getParameter(3));
								System.out.println(" Par 4 " +fun.getParameter(4));
								System.out.println(" Par 5 " +fun.getParameter(5));
								System.out.println(" Par 6 " +fun.getParameter(6));
								//System.out.println(fun.getParameterEstimate());
								
								//Try to fix the fit:
								double valuep5= fun.getParameter(5); 
								double valuep6 =fun.getParameter(6);
								double valuep4= fun.getParameter(4);
						
								int counter=0;
								while(valuep5<0.131 || valuep5>0.139 || valuep6<0.01 || valuep6>0.02 || valuep4<0   ) {
									double sign =0;
									fun = new F1D("f1d","[p0]+[p1]*x +[p2]*x*x+[p3]*x*x*x+[a]*gaus(x,[m],[s])",min_pion+counter*(0.00002),max_pion_fit-counter*(0.00002));
									fun.setLineColor(4);
									fun.setLineWidth(5);
									fun.setParameter(0,(totEntry/Mass_Bins) * Math.random());
									if (Math.random()>0.5)sign=-1.0;
									else sign=+1;
									fun.setParameter(1, 10 *sign* Math.random());
									if (Math.random()>0.5)sign=-1.0;
									else sign=+1;
									fun.setParameter(2,10*sign* Math.random());
									if (Math.random()>0.5)sign=-1.0;
									else sign=+1;
									fun.setParameter(3,-10 *sign* Math.random());
									fun.setParameter(4, 100 * Math.random());
									DataFitter.fit(fun, PionMass.GetH1F(i, k), "Q");
									valuep4=fun.getParameter(4);
									valuep5=fun.getParameter(5);
									valuep6=fun.getParameter(6);
									//System.out.println("I AM HEReee! Parameter 4 " +valuep4 + " Parameter 5 " +  valuep5 + " Parameter 6" + valuep6);
									fun.setParameter(5, 0.135);
									fun.setParameter(6, 0.013);
									
									counter++;
									if (counter>1000) {
										badZ =i ;
										badPT=k ; 
										 badFit=true;
										break;
									
									}
								}
								
								functions.add(fun);	
								
							}
							else if (fit==2) {
								fun2 = new F1D("f1d","[p0]+[p1]*x +[p2]*x*x+[p3]*x*x*x+[p4]*x*x*x*x+[a]*gaus(x,[m],[s])",min_pion,max_pion_fit);
								fun2.setParameter(0, totEntry/Mass_Bins);
								fun2.setParameter(0, 100);
								fun2.setParameter(1, -1);
								fun2.setParameter(2, 1);
								fun2.setParameter(3, 2);
								fun2.setParameter(4, 0.5);
								fun2.setParameter(5,  totEntry/Mass_Bins);
								fun2.setParameter(6, 0.133);
								fun2.setParameter(7, 0.013);
								fun2.setLineColor(4);
								fun2.setLineWidth(5);
								DataFitter.fit(fun2, PionMass.GetH1F(i, k), "Q");
								functions.add(fun2);
							}

							// Here I need to Save my function somewhere on paper and add into the histos
							// System.out.println("  -> I have fitted");
							int ZetaCounts=0; // counter for nr of pi0 in each bin after BKG subtration
							int GCounts =0 ;// counter for nr of pi0 from Gaussian Fit 
							double interval = ((max_pion-min_pion)/Mass_Bins);
							if(fit==1) {
								for(int qq=0 ; qq < Mass_Bins; qq++ ){
									double bin_shift = qq*(interval+(max_pion-min_pion)/Mass_Bins)/2+min_pion;
									if( bin_shift <(pi_peak+sig*EC_res) && bin_shift > (pi_peak-sig*EC_res)  ){
										double Background=fun.getParameter(0)+fun.getParameter(1)*bin_shift+fun.getParameter(2)*Math.pow(bin_shift,2)+fun.getParameter(3)*Math.pow(bin_shift,3);
										ZetaCounts= (int) (ZetaCounts+ PionMass.GetH1F(i, k).getBinContent(qq)-Background);
									}
								}
								GCounts = (int) Math.abs(fun.getParameter(4)*(fun.getParameter(6)/interval)*Math.sqrt(2*Math.PI));
							}

							else if(fit==2) {
								for(int qq=0 ; qq < Mass_Bins; qq++ ){
									double bin_shift = qq*(interval+(max_pion-min_pion)/Mass_Bins)/2+min_pion;
									if( bin_shift <(pi_peak+sig*EC_res) && bin_shift > (pi_peak-sig*EC_res)  ){
										double Background=fun2.getParameter(0)+fun2.getParameter(1)*bin_shift+fun2.getParameter(2)*Math.pow(bin_shift,2)+fun2.getParameter(3)*Math.pow(bin_shift,3)+fun2.getParameter(4)*Math.pow(bin_shift,4);
										ZetaCounts= (int) (ZetaCounts+ PionMass.GetH1F(i, k).getBinContent(qq)-Background);
									}
								}
								GCounts = (int) Math.abs(fun2.getParameter(5)*(fun2.getParameter(7)/interval)*Math.sqrt(2*Math.PI));
							}
							GaussianCounts.GetH2F(x,y).setBinContent(i, k, GCounts);
							BinCounts.GetH2F(x, y).setBinContent(i, k, ZetaCounts);
							directory.addDataSet(PionMass.GetH1F(i, k)); // I can use this one to get back my miss fitting.
							directory.addDataSet(xMass);
							directory.addDataSet(MissingMassX.GetH1F(i, k)); //MissingMass
							if (mC==false) {
								directory.addDataSet(PhiBSA.GetH1F(i, k));
								directory.addDataSet(HelicityPositive.GetH1F(i, k));
								directory.addDataSet(HelicityNegative.GetH1F(i, k));
								directory.addDataSet(Dist_z.GetH1F(i, k));
								directory.addDataSet(Dist_Pt.GetH1F(i, k));
							}

							// directory write file was here 
							System.out.println("electrons "+ Electron_counts.GetH2().getBinContent(x, y));
							System.out.println(" Bin Z  "+ i +" Bin PT  " +  k   +"  -> Pi0 from Gauss: " + GCounts + " Pi0 from Histo-BKG: " + ZetaCounts ); 
							if(badFit==true) {
							outP.println(" - > The bad fits occurs at bins " + " Z " +i + "PT " +k+ " Please use the PionMass histoigram in Z/PT folders to refit");
							badFit=false; // reset the badFit boolean 
							}
							outP.println("electrons "+ Electron_counts.GetH2().getBinContent(x, y));
							outP.println(" Bin Z  "+ i +" Bin PT  " + k   +"  -> Pi0 from Gauss: " + GCounts + " Pi0 from Histo-BKG: " + ZetaCounts );; 					
							photons.cd(myindex); BSA.cd(myindex);
							photons.draw(PionMass.GetH1F(i, k)); 
							if (mC==false)		BSA.draw(PhiBSA.GetH1F(i, k));
							String strg0 = String.format("%s/InvMass_Pi0s.png",PathFolder);
							String strg1 = String.format("%s/BSA_Pi0s.png",PathFolder);
							photons.save(strg0); 
							if (mC==false)BSA.save(strg1); 
							myindex++;
						}//Z PT
					}//PT bin
					directory.writeFile(PathFolder+"/Histograms_SIDIS.hipo");
					//Add the distribution of the Pion
					directoryMain.addDataSet(Electron_counts.GetH2());
					directoryMain.addDataSet(GaussianCounts.GetH2F(x,y));
					directoryMain.addDataSet(BinCounts.GetH2F(x,y));
					directoryMain.addDataSet(Pion_Energy.GetH1F(x, y));
					directoryMain.addDataSet(Dist_Q2.GetH1F(x, y));
					directoryMain.addDataSet(Dist_xB.GetH1F(x, y));
					directoryMain.addDataSet(Dist_W.GetH1F(x, y));
					directoryMain.addDataSet(Dist_y.GetH1F(x, y));
					directoryMain.addDataSet(Dist_eps.GetH1F(x, y));
					directoryMain.writeFile(PathFolder+"/HistogramsCounts.hipo");
					
					outP.close();						    // Printing stuff: 	
				} //File exitst

			} //Xb bins 
		} //Q2 bins


		return;
	}


	public void CountPions( H1F xMass, MultiBins MissingM, MultiBins PionCounts,MultiBins Counts_Phi, MultiBins H0_Counts, MultiBins H1_Counts, MultiBins Electron_counts, String time, boolean mC, String workdirout) throws FileNotFoundException {

		
		List<F1D> functions = new ArrayList<F1D>(); 
		//System.out.println(" Valore Histogramma " + this.InvMass.getSlicesZ().get(20).projectionX().getEntries());
		double totEntry=0;
System.out.println(" ========================= SAVING ANALYSIS RECON COUNTS ================");
		
		int Phi_Bins= H0_Counts.FthBin;
		double min_phi = 0.0;
		double max_phi = 360.0;

		System.out.println(" Phi_Bins " + Phi_Bins);
		System.out.println(" Min Phi " + min_phi + " Max Phi "+ max_phi);

		String s1=workdirout;
		new File(s1+"/Plots").mkdir();
		// TDirectory was here 
		
		String PathFolder = s1+"/Plots/"+time;
		File dir = new File(PathFolder);dir.mkdir();
		PathFolder=s1+"/Plots/"+time+"/MultiBins";
		new File (PathFolder).mkdir();
		System.out.println(" (Pi) - > I am Plotting Multi Bins Histograms");
		// Looping over all bins and fill the 3Dimensional Invariant Mass bins (z,pt,mass) by using the 5dimensional information of InvariantMassBins
		for(int x=0 ; x <(PionCounts.FstBin);x++) { //Q2bins
			String number = Integer.toString(x);
			PathFolder=s1+"/Plots/"+time+"/MultiBins/"+PionCounts.GetName_1stVariable()+ "_"+number;
			File dirbin = new File(PathFolder);dirbin.mkdir();
			for (int y=0 ; y<(PionCounts.SndBin); y++) { //Xb bins
				
				TDirectory directory = new TDirectory();
				TDirectory directoryMain = new TDirectory();
//to put it before it 
		
			
				
// Count particle Phi counts
				MultiBins PhiCounts = new MultiBins(Counts_Phi.TrdBin,Counts_Phi.FrtBin,Counts_Phi.FthBin, 0 ,0);
				PhiCounts.SetUnevenBins(Counts_Phi.GetUnevenBins(3), Counts_Phi.GetUnevenBins(4),Counts_Phi.GetUnevenBins(5), new double[]{0,0},new double[]{ 0,0});
			//	PhiCounts.SetBins(Counts_Phi.Trd_Min, Counts_Phi.Trd_Max, Counts_Phi.Frt_Min, Counts_Phi.Frt_Max,  Counts_Phi.Fth_Min, Counts_Phi.Fth_Max, 0, 0, 0, 0);
				PhiCounts.GenerateHistograms("PhiCounts");
				
				
				// Producing the H1F for the invariant mass computation 
				MultiBins PhiBSA = new MultiBins(H0_Counts.TrdBin,H0_Counts.FrtBin,H0_Counts.FthBin, 0 ,0);
				PhiBSA.SetUnevenBins(H0_Counts.GetUnevenBins(3), H0_Counts.GetUnevenBins(4),H0_Counts.GetUnevenBins(5), new double[]{0,0},new double[]{ 0,0});
				//PhiBSA.SetBins(H0_Counts.Trd_Min, H0_Counts.Trd_Max, H0_Counts.Frt_Min, H0_Counts.Frt_Max, H0_Counts.Fth_Min, H0_Counts.Fth_Max, 0, 0, 0, 0);
				PhiBSA.GenerateHistograms("PhiBSA");

				// here I have to count the stuff:
				// Producing the H1F for the invariant mass computation 
				MultiBins HelicityPositive = new MultiBins(H1_Counts.TrdBin,H1_Counts.FrtBin,H1_Counts.FthBin, 0 ,0);
				HelicityPositive.SetUnevenBins(H1_Counts.GetUnevenBins(3), H1_Counts.GetUnevenBins(4),H1_Counts.GetUnevenBins(5), new double[]{0,0},new double[]{ 0,0});
				//HelicityPositive.SetBins(H1_Counts.Trd_Min, H1_Counts.Trd_Max, H1_Counts.Frt_Min, H1_Counts.Frt_Max, H1_Counts.Fth_Min, H1_Counts.Fth_Max, 0, 0, 0, 0);
				HelicityPositive.GenerateHistograms("HelicityPositive");

			
				MultiBins HelicityNegative = new MultiBins(H0_Counts.TrdBin,H0_Counts.FrtBin,H0_Counts.FthBin, 0 ,0);
				HelicityNegative.SetUnevenBins(H0_Counts.GetUnevenBins(3), H0_Counts.GetUnevenBins(4),H0_Counts.GetUnevenBins(5), new double[]{0,0},new double[]{ 0,0});
				//HelicityNegative.SetBins(H0_Counts.Trd_Min, H0_Counts.Trd_Max, H0_Counts.Frt_Min, H0_Counts.Frt_Max, H0_Counts.Fth_Min, H0_Counts.Fth_Max, 0, 0, 0, 0);
				HelicityNegative.GenerateHistograms("HelicityNegative");

				
				MultiBins MissingMassX = new MultiBins(MissingM.TrdBin,MissingM.FrtBin,MissingM.FthBin,0,0);
				MissingMassX.SetUnevenBins(MissingM.GetUnevenBins(3), MissingM.GetUnevenBins(4),MissingM.GetUnevenBins(5), new double[]{0,0},new double[]{ 0,0});
				//MissingMassX.SetBins(MissingM.Trd_Min, MissingM.Trd_Max, MissingM.Frt_Min, MissingM.Frt_Max, MissingM.Fth_Min, MissingM.Fth_Max, 0, 0, 0, 0);
				MissingMassX.GenerateHistograms("MissingMass");


				String number2=Integer.toString(y);
				PathFolder=s1+"/Plots/"+time+"/MultiBins/"+PionCounts.GetName_1stVariable()+ "_"+number+"/"+PionCounts.GetName_2ndVariable()+ "_"+number2;
				File dirbin2 = new File(PathFolder);dirbin2.mkdir();
				File File_txt = new File(PathFolder+"/results.txt");
				directoryMain.mkdir("/main/"); directoryMain.cd("/main/");
				try {
					File_txt.createNewFile();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(File_txt.exists()) {
					System.out.println(" Output File Found" );
					PrintWriter outP = new PrintWriter(File_txt);
					EmbeddedCanvas BSA = new EmbeddedCanvas();
					EmbeddedCanvas CountsH1 = new EmbeddedCanvas();
					BSA.divide(H1_Counts.TrdBin,H1_Counts.FrtBin);
					CountsH1.divide(H1_Counts.TrdBin,H1_Counts.FrtBin);
					int myindex=0;

					for(int i = 0; i < (PionCounts.TrdBin); i++) { //Z bins
						for(int k = 0;k < (PionCounts.FrtBin); k++) {  // pt
							
							directory.mkdir("/Z"+i+"PT"+k+"/"); directory.cd("/Z"+i+"PT"+k+"/");
							PhiBSA.GetH1F(i, k).setTitle("BSA Bin Z: " + i + " Bin Pt: "+ k);
							PhiBSA.GetH1F(i, k).setTitleX("Phi [degree]");

							for(int ii=0; ii < MissingM.FthBin; ii++) {

								MissingMassX.GetH1F(i, k).fill(MissingM.GetH3F(x, y).getBinContent(i, k, ii));
								//System.out.println("Missing mass Value "+MissingM.GetH3F(x, y).getBinContent(i, k, ii) );
							}
							double Integrated_HP =0;
							double Integrated_HM=0;
							if (mC==false) {
								for(int ll=0; ll<Phi_Bins; ll++) {
									
									double binsize=(max_phi-min_phi)/Phi_Bins;
									double binPhi = (ll)*(binsize)+(binsize/2)+min_phi;
									
							
									double num= H1_Counts.GetH3F(x, y).getBinContent(i, k, ll) - H0_Counts.GetH3F(x, y).getBinContent(i, k, ll);
									double den= H0_Counts.GetH3F(x, y).getBinContent(i, k, ll) + H1_Counts.GetH3F(x, y).getBinContent(i, k, ll);
									double Bsa= num/den;
									PhiBSA.GetH1F(i, k).fill(binPhi,Bsa);
									HelicityNegative.GetH1F(i, k).fill(binPhi,H0_Counts.GetH3F(x, y).getBinContent(i, k, ll));
									HelicityPositive.GetH1F(i, k).fill(binPhi,H1_Counts.GetH3F(x, y).getBinContent(i, k, ll));
									Integrated_HP+=H1_Counts.GetH3F(x, y).getBinContent(i, k, ll);
									Integrated_HM+=H0_Counts.GetH3F(x, y).getBinContent(i, k, ll);
								}// Phi bins 
								
								directory.addDataSet(PhiBSA.GetH1F(i, k));
								directory.addDataSet(HelicityPositive.GetH1F(i, k));
								directory.addDataSet(HelicityNegative.GetH1F(i, k));
								directory.addDataSet(MissingMassX.GetH1F(i, k));//MissingMass
								directory.addDataSet(Dist_z.GetH1F(i, k));
								directory.addDataSet(Dist_Pt.GetH1F(i, k));


							}
							else {
								// For MC == true; 
								for(int ll=0; ll<Phi_Bins; ll++) {
									double binsize=(max_phi-min_phi)/Phi_Bins;
									double binPhi = (ll)*(binsize)+(binsize/2)+min_phi;
									PhiCounts.GetH1F(i, k).fill(binPhi,Counts_Phi.GetH3F(x, y).getBinContent(i, k, ll));
									
								}
								directory.addDataSet(PhiCounts.GetH1F(i, k));
							}
							outP.println("electrons "+ Electron_counts.GetH2().getBinContent(x, y));
							outP.println(" Bin Z  "+ i +" Bin PT  " + k   +"  -> Charged Pions from EB: " + PionCounts.GetH2F(x, y).getBinContent(i, k) );			
						  outP.println(" BinZ " + i + " Bin PT " +k + " - > Phi integrated Helicity P " +Integrated_HP + " Phi Integrated Helicity M" + Integrated_HM);
							BSA.cd(myindex); CountsH1.cd(myindex);
							
							myindex++;
						}// PT binz
					}//Z bins
					//It was here
					System.out.println("----- ");
					System.out.println(" Pion_nergy naame" + Pion_Energy.GetNameHistos());
					directory.writeFile(PathFolder+"/Histograms_SIDIS.hipo");
					directoryMain.addDataSet(Electron_counts.GetH2());
					directoryMain.addDataSet(PionCounts.GetH2F(x,y));
					directoryMain.addDataSet(xMass);
					directoryMain.addDataSet(Pion_Energy.GetH1F(x, y));
					directoryMain.addDataSet(Pion_Theta.GetH1F(x, y));
					directoryMain.addDataSet(Pion_Momentum.GetH1F(x, y));
					directoryMain.addDataSet(Pion_Z.GetH1F(x, y));
					directoryMain.addDataSet(Pion_Phi.GetH1F(x, y));
					directoryMain.addDataSet(Pion_PT.GetH1F(x, y));
					directoryMain.addDataSet(Pion_PhiClas.GetH1F(x, y));
					directoryMain.addDataSet(Dist_Q2.GetH1F(x, y));
					directoryMain.addDataSet(Dist_xB.GetH1F(x, y));
					directoryMain.addDataSet(Dist_W.GetH1F(x, y));
					directoryMain.addDataSet(Dist_y.GetH1F(x, y));
					directoryMain.addDataSet(Dist_eps.GetH1F(x, y));
					//System.out.println("The pion Energy entries are " +Pion_Energy.GetH1F(x, y).getEntries());
					directoryMain.writeFile(PathFolder+"/HistogramsCounts.hipo");
					outP.close();						    // Printing stuff: 	
				} //File exitst



			} //Xb bins 
		} //Q2 bins


		return;
	}





	/**
	 * Save multi dimensional BSA
	 */
	public void Save_BSA() {


	}
	public void Histo (String time)
	{

		String s1="/Users/gangelini/Desktop/";
		new File(s1+"/Plots").mkdir();
		String PathFolder = s1+"/Plots/"+time;
		File dir = new File(PathFolder);
		dir.mkdir();
		PathFolder=s1+"/Plots/"+time+"/MultiBins";
		new File (PathFolder).mkdir();
		//new File("/Users/gangelini/Desktop/Plots/Inv_Mass").mkdir();


		System.out.println(" (Pi) - > I am Plotting Invariant Mass Histograms");

		EmbeddedCanvas photons = new EmbeddedCanvas();
		photons.setSize(1600,1000);
		photons.divide(PT_Bins,Z_Bins);
		photons.setAxisTitleSize(20);
		photons.setAxisFontSize(20);
		photons.setTitleSize(20);
		// I plot z from  from 0.1 to 0.9 because first and last bin are useless
		for(int i=0 ; i<(Z_Bins*PT_Bins); i++) {
			photons.cd(i);photons.draw(InvMassPlots.get(i));
		}
		for (int ii=0; ii<Z_Bins;ii++) {
			System.out.println(" Z bins:");
			System.out.println(this.InvMass.getXAxis().getBinCenter(ii));
		}
		for (int jj=0; jj<PT_Bins;jj++) {
			System.out.println(" PT bins:");
			System.out.println(this.InvMass.getYAxis().getBinCenter(jj));
		}


		String strg0 = String.format("%s/Pi0s.png",PathFolder);
		System.out.println("Saving plots in "+ PathFolder);
		photons.save(strg0);	

		PathFolder=s1+"/Plots/"+time+"/Distribution";
		new File (PathFolder).mkdir();
		//new File("/Users/gangelini/Desktop/Plots/Inv_Mass").mkdir();


		System.out.println(" (Pi) - > I am Plotting Invariant Mass Histograms");

		EmbeddedCanvas distribution = new EmbeddedCanvas();
		distribution.setSize(1600,1000);
		distribution.divide(3,1);
		distribution.setAxisTitleSize(20);
		distribution.setAxisFontSize(20);
		distribution.setTitleSize(20);
		// I plot z from  from 0.1 to 0.9 because first and last bin are useless

		distribution.cd(0);distribution.draw(Z_Dist_Overall);
		distribution.cd(1);distribution.draw(Pt_Dist_Overall);
		distribution.cd(2);distribution.draw(XB_Dist_Overall);
		String strg1= String.format("%s/Distributions.png",PathFolder);
		System.out.println("Saving plots in "+ PathFolder);
		distribution.save(strg1);	



		
}

}
