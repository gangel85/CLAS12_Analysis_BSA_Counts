package GoodAnalysis;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.H3F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.jnp.physics.Particle;
import org.jlab.jnp.physics.Vector3;





public class Histos {

	private int DCHIT=60; // I am taking the x,y,z on DC1 here 
	private int DC1=61; //The index used to identify the first DC Region
	private int DC2=62 ;//The index used to identify the second DC region 
	private int DC3=63; //The index used to identify the third DC region 	 

	private int PCALHIT=70; // I am taking the x,y,z,on PCAL 
	private int PCALID=71; //The index used to idenfity the PCAL
	private int ECAL1ID=74; //The index used to idenfity the PCAL
	private int ECAL2ID=77; //The index used to idenfity the PCAL



	private String s1 = " ";
	//public Electron Selection Electron;
	private static H1F H_lu_eBC,H_lv_eBC,H_lw_eBC;
	private static H2F H_e_xB_Q2BC, H_e_y_Q2BC, H_e_W_Q2BC,H_e_vx_vyBC,H_ESampl_ECalBC;
	private static H1F HQ2BC,H_e_xBC,H_e_yBC,H_e_WBC,H_e_vzBC,H_e_ThetaBC,H_e_PhiBC,H_e_MomBC;


	// Electron Histograms  AC
	private static H1F H_lu_eAC,H_lv_eAC,H_lw_eAC;
	private static H2F H_e_xB_Q2AC, H_e_y_Q2AC, H_e_W_Q2AC,H_e_vx_vyAC,H_ESampl_ECalAC;
	private static H1F HQ2AC,H_e_xAC,H_e_yAC,H_e_WAC,H_e_vzAC,H_e_ThetaAC,H_e_PhiAC,H_e_MomAC;

	// Other stuff
	private static H1F H_e_Theta_DIS,H_e_Theta_F, H_e_Mom_DIS,H_e_Mom_F,H_eP_Theta_DIS,H_eP_Mom_DIS;
	private static H2F H_Ecal_Escin,H_ECal_xy,H_ECal_xyNCut,H_DC_xy,H_DC_xyNCut;

	// Differences MC to REC
	private static  H2F H_Q2_difference, H_xB_difference, H_W_difference, H_Y_difference  ;
	private static H1F H_Diff_Q2, H_Diff_xB , H_Diff_W, H_Diff_Y;
	//mixed plots
	private static H2F H_xb_Q2_difference, H_W_Q2_difference, H_Y_Q2_difference;

	// Electron Histogram AC

	// Photons Histograms:
	private static H1F H_Phi_difference, H_PhiT_difference;
	private static H1F H_ph_pxBC,H_ph_pyBC,H_ph_pzBC,H_ph_MomBC,H_ph_ThetaBC,H_ph_PhiBC;
	private static H1F H_ph_pxAC,H_ph_pyAC,H_ph_pzAC,H_ph_MomAC,H_ph_ThetaAC,H_ph_PhiAC;



	// Charged Particle Histograms 
	private static H1F H_z_difference,H_PT_difference,H_Theta_difference, H_Energy_difference;
	private static H2F H_dz_z,H_dz_Pt,H_dz_phi,H_dz_phiT,H_dz_theta,H_dPt_Pt,H_dPt_z,H_dPt_phi,H_dPt_phiT,H_dPt_theta; 
	private static H2F H_dE_E ,H_dE_theta ,H_dE_phi,H_dE_phiT,H_PhiT_z,H_PhiT_phi,H_PhiT_theta,H_PhiT_phiT,H_PhiT_Pt;
	private static H2F H_dTheta_e, H_dTheta_theta, H_dTheta_phi, H_dPhi_e, H_dPhi_theta, H_dPhi_phi;
	

	private static H1F H_generated_p, H_reconstructed_p;
	
	
	//  Electron resolutions:
	H1F H_e_dep,H_e_dtheta,H_e_dphi,H_e_dQ2,H_e_dxB,H_e_dW;
	H2F H_e_dp_p , H_e_dp_theta, H_e_dp_phi, H_e_dtheta_p, H_e_dtheta_theta,H_e_dtheta_phi,H_e_dphi_p,H_e_dphi_theta,H_e_dphi_phi;
	H2F H_e_dQ2_Q2,H_e_dQ2_xB,H_e_dQ2_W,H_e_dxB_Q2,H_e_dxB_xB,H_e_dxB_W,H_e_dW_Q2,H_e_dW_xB,H_e_dW_W;

	H2F Gen_Q2_xB,Gen_Q2_W, Gen_xB_W,Rec_Q2_xB,Rec_Q2_W, Rec_xB_W, Gen_e_theta, Gen_e_phi, Gen_theta_phi, Rec_e_theta, Rec_e_phi, Rec_theta_phi;





	public Histos( String SavingDir )
	{
		Plots_el();
		Plots_ph();
		s1 = SavingDir;	
		// Inizialize all the histograms right here 
	}

	private void Plots_el() {
		//Electrons 
		H_e_xB_Q2BC = new H2F("H_e_xB_Q2","H_e_xB_Q2",100,0,1,100,0,12);
		H_e_xB_Q2BC.setTitle("Q^2 vs xB Before Cuts");
		H_e_xB_Q2BC.setTitleX("xB");
		H_e_xB_Q2BC.setTitleY("Q^2");

		H_e_xB_Q2AC = new H2F("H_e_xB_Q2","H_e_xB_Q2",100,0,1,100,0,12);
		H_e_xB_Q2AC.setTitle("Q^2 vs xB After Cuts");
		H_e_xB_Q2AC.setTitleX("xB");
		H_e_xB_Q2AC.setTitleY("Q^2");

		H_e_y_Q2BC = new H2F("H_e_y_Q2","H_e_xB_Q2",100,0.2,1,100,0,12);
		H_e_y_Q2BC.setTitle("Q^2 vs y Before Cuts");
		H_e_y_Q2BC.setTitleX("y");
		H_e_y_Q2BC.setTitleY("Q^2");

		H_e_y_Q2AC = new H2F("H_e_y_Q2","H_e_xB_Q2",100,0.2,1,100,0,12);
		H_e_y_Q2AC.setTitle("Q^2 vs y After Cuts");
		H_e_y_Q2AC.setTitleX("y");
		H_e_y_Q2AC.setTitleY("Q^2");

		H_e_W_Q2BC = new H2F("H_e_W_Q2","H_e_W_Q2",100,1,5,100,1,12);
		H_e_W_Q2BC.setTitle("Q^2 vs W Before Cuts");
		H_e_W_Q2BC.setTitleX("W");
		H_e_W_Q2BC.setTitleY("Q^2");

		H_e_W_Q2AC = new H2F("H_e_W_Q2","H_e_W_Q2",100,1,5,100,1,12);
		H_e_W_Q2AC.setTitle("Q^2 vs W After Cuts");
		H_e_W_Q2AC.setTitleX("W");
		H_e_W_Q2AC.setTitleY("Q^2");

		HQ2BC = new H1F("HQ2","HQ2",100,0,12);
		HQ2BC.setTitle("Q2 Distribution Before Cuts");
		HQ2BC.setTitleX("Q^2 [GeV^2]");

		HQ2AC = new H1F("HQ2","HQ2",100,0,12);
		HQ2AC.setTitle("Q2 Distribution After Cuts");
		HQ2AC.setTitleX("Q^2 [GeV^2]");

		H_e_xB_Q2AC = new H2F("H_e_xB_Q2AC","H_e_xB_Q2AC",100,0,1,100,0,12);
		H_e_xB_Q2AC.setTitle("Q^2 vs xB After Cut");
		H_e_xB_Q2AC.setTitleX("xB");
		H_e_xB_Q2AC.setTitleY("Q^2");

		H_e_y_Q2AC = new H2F("H_e_y_Q2AC","H_e_xB_Q2AC",100,0.2,1,100,0,12);
		H_e_y_Q2AC.setTitle("Q^2 vs y After cut");
		H_e_y_Q2AC.setTitleX("y");
		H_e_y_Q2AC.setTitleY("Q^2");

		H_e_W_Q2AC = new H2F("H_e_W_Q2AC","H_e_W_Q2AC",100,1,5,100,1,12);
		H_e_W_Q2AC.setTitle("Q^2 vs W After cut ");
		H_e_W_Q2AC.setTitleX("W");
		H_e_W_Q2AC.setTitleY("Q^2");

		H_e_WBC = new H1F("H_e_W","H_e_W",100,1.5,4.5);
		H_e_WBC.setTitle("Electron W Before Cuts");
		H_e_WBC.setTitleX("W");
		H_e_WBC.setFillColor(2);


		H_e_WAC = new H1F("H_e_W","H_e_W",100,1.5,4.5);
		H_e_WAC.setTitle("Electron W After Cuts");
		H_e_WAC.setTitleX("W");
		H_e_WAC.setFillColor(2);

		H_e_vzBC = new H1F("H_e_vz","H_e_vz",200,-20,20);
		H_e_vzBC.setTitle("Electron longitudinal vertex Before Cuts");
		H_e_vzBC.setTitleX("v_{z} (cm)");
		H_e_vzBC.setFillColor(2);

		H_e_vzAC = new H1F("H_e_vz","H_e_vz",200,-20,20);
		H_e_vzAC.setTitle("Electron longitudinal vertex After  Cuts");
		H_e_vzAC.setTitleX("v_{z} (cm)");
		H_e_vzAC.setFillColor(2);

		H_e_vx_vyBC = new H2F("H_e_vx_vy ", "H_e_vx_vy ", 100, -20,20, 100, -20, 20 );
		H_e_vx_vyBC.setTitle("Electron vertex x vs y Before Cuts");
		H_e_vx_vyBC.setTitleX(" El. x vertex [cm]");
		H_e_vx_vyBC.setTitleY(" El. y vertex [cm]");

		H_e_vx_vyAC = new H2F("H_e_vx_vy ", "H_e_vx_vy ", 100, -20,20, 100, -20, 20 );
		H_e_vx_vyAC.setTitle("Electron vertex x vs y After Cuts");
		H_e_vx_vyAC.setTitleX(" El. x vertex [cm]");
		H_e_vx_vyAC.setTitleY(" El. y vertex [cm]");			

		H_e_ThetaBC = new H1F("H_e_Theta","H_e_Theta",100,0,60);
		H_e_ThetaBC.setTitle(" Ange of electrons Before Cuts ");
		H_e_ThetaBC.setTitleX("Angle [theta]");
		H_e_ThetaBC.setFillColor(5);

		H_e_ThetaAC = new H1F("H_e_Theta","H_e_Theta",100,0,60);
		H_e_ThetaAC.setTitle(" Ange of electrons After Cuts ");
		H_e_ThetaAC.setTitleX("Angle [theta]");
		H_e_ThetaAC.setFillColor(5);

		H_e_PhiBC = new H1F("H_e_Theta","H_e_Theta",100,-180,180);
		H_e_PhiBC.setTitle(" Ange of electrons Before Cuts ");
		H_e_PhiBC.setTitleX("Angle [theta]");
		H_e_PhiBC.setFillColor(5);

		H_e_PhiAC = new H1F("H_e_Phi","H_e_Phi",100,-180,180);
		H_e_PhiAC.setTitle(" Az. Ange of electrons After Cuts ");
		H_e_PhiAC.setTitleX("Angle [Phi]");
		H_e_PhiAC.setFillColor(5);

		H_e_MomBC = new H1F("H_e_Mom","H_e_Mom",100,0,9);
		H_e_MomBC.setTitle(" Momentum of electrons Before Cuts ");
		H_e_MomBC.setTitleX("Momentum [GeV/c]");
		H_e_MomBC.setFillColor(5);

		H_e_MomAC = new H1F("H_e_Mom","H_e_Mom",100,0,9);
		H_e_MomAC.setTitle(" Momentum of electrons After Cuts ");
		H_e_MomAC.setTitleX("Momentum [GeV/c]");
		H_e_MomAC.setFillColor(5);

		H_e_Theta_F = new H1F("H_e_Theta_F","H_e_Theta_F",100,0,60);
		H_e_Theta_F.setTitle(" Ange of electrons before cuts");
		H_e_Theta_F.setTitleX("Angle [theta]");
		H_e_Theta_F.setFillColor(5);

		H_e_Mom_F = new H1F("H_e_Mom_F","H_e_Mom_F",100,0,9);
		H_e_Mom_F.setTitle(" Momentum of electrons before cuts");
		H_e_Mom_F.setTitleX("Momentum [GeV/c]");
		H_e_Mom_F.setFillColor(5);

		H_ESampl_ECalBC = new H2F("H_ESampl_ECal","H_ESampl_ECal",100,0,11,100,0,0.5);
		H_ESampl_ECalBC.setTitle("Electron ECAL Sampling Fraction Before Cuts ");
		H_ESampl_ECalBC.setTitleX("p (GeV/c)");
		H_ESampl_ECalBC.setTitleY("Edep/p");

		H_ESampl_ECalAC = new H2F("H_ESampl_ECal","H_ESampl_ECal",100,0,11,100,0,0.5);
		H_ESampl_ECalAC.setTitle("Electron ECAL Sampling Fraction After Cuts ");
		H_ESampl_ECalAC.setTitleX("p (GeV/c)");
		H_ESampl_ECalAC.setTitleY("Edep/p");

		H_Ecal_Escin = new H2F("H_Ecal_Escin","H_Ecal_Escin",100,0,11,100,0,110);
		H_Ecal_Escin.setTitle("Electron Ecal vs E scintillator");
		H_Ecal_Escin.setTitleX("Edep Cal ");
		H_Ecal_Escin.setTitleY("Edep Scin");

		H_ECal_xy = new H2F("H_Ecal_xy","H_ECal_xy",100,-500,500,100,-500,500 );
		H_ECal_xy.setTitle("Electron Ecal X and Y position ");
		H_ECal_xy.setTitleX("X [cm]");
		H_ECal_xy.setTitleY("Y [cm]");

		H_ECal_xyNCut = new H2F("H_Ecal_xyNCut","H_ECal_xyNCut",100,-500,500,100,-500,500 );
		H_ECal_xyNCut.setTitle("Electron Cutted Ecal X and Y position ");
		H_ECal_xyNCut.setTitleX("X [cm]");
		H_ECal_xyNCut.setTitleY("Y [cm]");

		H_DC_xy = new H2F("H_DC_xy","H_DC_xy",100,-200,200,100,-200,200 );
		H_DC_xy.setTitle("Electron DC X and Y position ");
		H_DC_xy.setTitleX("X [cm]");
		H_DC_xy.setTitleY("Y [cm]");

		H_DC_xyNCut = new H2F("H_DC_xyNCut","H_DC_xyNCut",100,-200,200,100,-200,200 );
		H_DC_xyNCut.setTitle("All electrons DC X and Y position ");
		H_DC_xyNCut.setTitleX("X [cm]");
		H_DC_xyNCut.setTitleY("Y [cm]");

		H_lu_eBC = new H1F("H_lu_eBC","H_lu_eBC", 100, 0 ,450);
		H_lu_eBC.setTitle("lu for e Before Cuts");
		H_lu_eBC.setTitleX("lu for e");
		H_lu_eBC.setFillColor(4);

		H_lu_eAC = new H1F("H_lu_eAC","H_lu_eAC", 100, 0 ,450);
		H_lu_eAC.setTitle("lu for e After Cuts");
		H_lu_eAC.setTitleX("lu for e");
		H_lu_eAC.setFillColor(4);

		H_lv_eBC = new H1F("H_lv_eBC","H_lv_eBC", 100, 0 ,450);
		H_lv_eBC.setTitle("lv for e Before Cuts");
		H_lv_eBC.setTitleX("lv for e");
		H_lv_eBC.setFillColor(4);

		H_lv_eAC = new H1F("H_lv_eAC","H_lv_eAC", 100, 0 ,450);
		H_lv_eAC.setTitle("lv for e After Cuts");
		H_lv_eAC.setTitleX("lv for e");
		H_lv_eAC.setFillColor(4);

		H_lw_eBC = new H1F("H_lw_eBC","H_lw_eBC", 100, 0 ,450);
		H_lw_eBC.setTitle("lw for e Before Cuts");
		H_lw_eBC.setTitleX("lw for e");
		H_lw_eBC.setFillColor(4);

		H_lw_eAC = new H1F("H_lw_eAC","H_lw_eAC", 100, 0 ,450);
		H_lw_eAC.setTitle("lw for e After Cuts");
		H_lw_eAC.setTitleX("lw for e");
		H_lw_eAC.setFillColor(4);

		H_e_WBC = new H1F ("H_e_WBC" , "H_e_WBC", 100, 1, 5) ;
		H_e_WBC.setTitle( " W distribution Before Cuts");
		H_e_WBC.setTitleX( " W ");

		H_e_WAC = new H1F ("H_e_WAC" , "H_e_WAC", 100, 1, 5) ;
		H_e_WAC.setTitle( " W distribution After Cuts");
		H_e_WAC.setTitleX( " W ");

		H_e_xBC = new H1F ("H_e_xBC" , "H_e_xBC", 100, 0, 1) ;
		H_e_xBC.setTitle( " xB distribution Before Cuts");
		H_e_xBC.setTitleX( " xB ");

		H_e_xAC = new H1F ("H_e_xAC" , "H_e_xAC", 100, 0, 1) ;
		H_e_xAC.setTitle( " xB distribution After Cuts");
		H_e_xAC.setTitleX( " xB ");

		H_e_yBC = new H1F ("H_e_yBC" , "H_e_yBC", 100, 0, 1) ;
		H_e_yBC.setTitle( " Y  distribution Before Cuts");
		H_e_yBC.setTitleX( " Y ");

		H_e_yAC = new H1F ("H_e_yAC" , "H_e_yAC", 100, 0, 1) ;
		H_e_yAC.setTitle( " Y  distribution After  Cuts");
		H_e_yAC.setTitleX( " Y ");

		H_Q2_difference = new H2F ("H_Q2_difference","H_Q2_difference", 100, 1,12,40, -4,4 );
		H_Q2_difference.setTitle(" Differences in Q2 Generated-Reconstructed");
		H_Q2_difference.setTitleY("Q2 Difference [GeV]");
		H_Q2_difference.setTitleX("Q2 Generated electron [GeV]");

		H_xB_difference = new H2F ("H_xB_difference","H_xB_difference",100,0,1,40, -0.4, 0.4);
		H_xB_difference.setTitle(" Differences in xB Generated-Reconstructed");
		H_xB_difference.setTitleY("xB Difference ");
		H_xB_difference.setTitleX(" xB of generated electron");

		 H_W_difference = new H2F ("H_W_difference","H_W_difference",100,0,5,40, -2, 2 );
		H_W_difference.setTitle(" Differences in W Generated-Reconstructed");
		H_W_difference.setTitleY("W Difference");
		H_W_difference.setTitleX(" W of generated electron ");

		 H_Y_difference = new H2F ("H_Y_difference","H_Y_difference",100,0,1,40, -0.4, 0.4  );
		H_Y_difference.setTitle("Differences in W Generated-Reconstructed");
		H_Y_difference.setTitleY("Y Difference ");
		H_Y_difference.setTitleX(" Y of generated electron");

		
		H_Diff_Q2 = new H1F ( "H_Diff_Q2 ", 50, -0.2, 0.2);
		H_Diff_Q2.setTitle("Q2 LUND - Q2 REC");
		H_Diff_Q2.setTitleX("Difference in Q2 [GeV]");
		H_Diff_Q2.setFillColor(3);
		
		H_Diff_xB = new H1F ("H_Diff_xB ", 50, -0.04, 0.04 );
		H_Diff_xB.setTitle("xB LUND - xB REC");
		H_Diff_xB.setTitleX("Difference in xB ");
		H_Diff_xB.setFillColor(4);
		
		H_Diff_W = new H1F ("H_Diff_W ", 50, -0.15, 0.15 );
		H_Diff_W.setTitle("W LUND - W REC");
		H_Diff_W.setTitleX("Difference in W");
		H_Diff_W.setFillColor(5);
		
		H_Diff_Y = new H1F ("H_Diff_Y ", 50, -0.03, 0.03);
		H_Diff_Y.setTitle("Y LUND - Y REC");
		H_Diff_Y.setTitleX("Difference in Y ");
		H_Diff_Y.setFillColor(6);
	
		H_z_difference = new H1F ("H_z_difference", 200, -0.06, 0.06);
		H_z_difference.setTitle("Z LUND - Z REC");
		H_z_difference.setTitleX("Difference in Z ");
		H_z_difference.setFillColor(3);
		
		H_PT_difference = new H1F ("H_PT_difference", 200, -0.15, 0.15);
		H_PT_difference.setTitle("PT LUND - PT REC");
		H_PT_difference.setTitleX("Difference in PT ");
		H_PT_difference.setFillColor(4);
		
		
		H_Theta_difference = new H1F ("H_Theta_difference", 200, -8 , 8 );
		H_Theta_difference.setTitle("Theta - Theta REC");
		H_Theta_difference.setTitleX("Difference in Theta ");
		H_Theta_difference.setFillColor(5);
		
		H_Phi_difference = new H1F ("H_Phi_difference", 200, -20 , 20 );
		H_Phi_difference.setTitle("Phi - Phi REC");
		H_Phi_difference.setTitleX("Difference in Phi ");
		H_Phi_difference.setFillColor(6);
		
		H_PhiT_difference = new H1F ("H_PhiT_difference", 200, -50 , 50 );
		H_PhiT_difference.setTitle("Phi trento - Phi REC trento");
		H_PhiT_difference.setTitleX("Difference in Phi Trento ");
		H_PhiT_difference.setFillColor(7);
		
		H_Energy_difference = new H1F ("H_Energy_difference", 100, -0.1, 0.1);
		H_Energy_difference.setTitle("Energy Lund - Energy REC");
		H_Energy_difference.setTitleX("Difference in Energy ");
		H_Energy_difference.setFillColor(8);
		
		//mixed plots
		
		H_dz_z = new H2F("H_dz_z",100,0,1,100,-0.06,0.06);
		H_dz_z.setTitle("Δz vs z ");
		H_dz_z.setTitleX("z");
		H_dz_z.setTitleY("Δz");
		
		H_dz_Pt = new H2F("H_dz_Pt",100,0,1.5,100,-0.06,0.06);
		H_dz_Pt.setTitle("Δz vs PT^2 ");
		H_dz_Pt.setTitleX("PT^2 [GeV^2]");
		H_dz_Pt.setTitleY("Δz");
		
		H_dz_phi = new H2F("H_dz_phi", 100,-180,180,100, -0.06,0.06);
		
		
		H_dz_theta = new H2F("H_dz_theta", 100, 0,40,100,-0.06,0.06);
		
		
		H_dz_phiT = new H2F("H_dz_phiT", 100, 0,360, 100, -0.06,0.06);
		H_dz_phiT.setTitle("Δz vs φ_{trento} ");
		H_dz_phiT.setTitleX(" Φ_{trento} [degrees]");
		H_dz_phiT.setTitleY(" Δz ");
		
		
		H_dPt_Pt= new H2F("H_dPt_Pt",100,0,1.6,100,-0.2,0.2);
		H_dPt_Pt.setTitle("ΔPT^2 vs PT^2 ");
		H_dPt_Pt.setTitleX("PT^2 [GeV^2]");
		H_dPt_Pt.setTitleY("ΔPT^2[GeV^2]");
		
		
		H_dPt_z = new H2F("H_dPt_z",100,0,1,100,-0.2,0.2);
		H_dPt_z.setTitle("ΔPT^2 vs z ");
		H_dPt_z.setTitleX(" z ");
		H_dPt_z.setTitleY("ΔPT^2[GeV^2]");
		
		H_dPt_phi = new H2F("H_dz_phi",100,-180,180,100,-0.2,0.2);
		H_dPt_theta = new H2F("H_dz_theta",100,0,40,100,-0.2,0.2);
		H_dPt_phiT = new H2F("H_dPt_phiT",100,0,360,100,-0.2,0.2);
		H_dPt_phiT.setTitle("ΔPT^2 vs Φ_{trento}");
		H_dPt_phiT.setTitleX(" Φ_{trento} [degrees] ");
		H_dPt_phiT.setTitleY("ΔPT^2 [GeV^2]");
	
		H_dE_E = new H2F("H_dE_E",100,0,8,100,-0.15,0.15);
		H_dE_E.setTitle("ΔE vs E");
		H_dE_E.setTitleX(" E [GeV] ");
		H_dE_E.setTitleY("ΔE [GeV]");
		
		H_dE_theta = new H2F("H_dE_theta", 100,0,45,100, -0.15,0.15);
		H_dE_theta.setTitle("ΔE vs θ");
		H_dE_theta.setTitleX(" θ [degrees] ");
		H_dE_theta.setTitleY("ΔE [GeV]");
		
		
		H_dE_phi = new H2F("H_dE_phi",100,-180,180,100,-0.15,0.15);
		H_dE_phi.setTitle("ΔE vs Φ");
		H_dE_phi.setTitleX(" Φ [degrees] ");
		H_dE_phi.setTitleY("ΔE [GeV]");
		H_dE_phiT = new H2F ("H_dE_phiT",100,0,360,100, -0.15,0.15);
		
		H_PhiT_Pt = new H2F("H_PhiT_Pt",100,0,1.6,100,-20,20);
		H_PhiT_Pt.setTitle("ΔΦ_{trento} vs PT^2 ");
		H_PhiT_Pt.setTitleX(" PT^2 [GeV^2] ");
		H_PhiT_Pt.setTitleY("ΔΦ_{trento} [degrees]");
		
		H_PhiT_z = new H2F("H_PhiT_z",100,0,1,100,-20,20);
		H_PhiT_z.setTitle("ΔΦ_{trento} vs z ");
		H_PhiT_z.setTitleX(" z ");
		H_PhiT_z.setTitleY("ΔΦ_{trento} [degrees]");
		
		H_PhiT_phi = new H2F("H_PhiT_phi",100,-180,180,100,-20,20);
		H_PhiT_theta = new H2F("H_PhiT_theta",100,0,45,100,-20,20);
		H_PhiT_phiT = new H2F("H_PhiT_phiT",100,0,360,100,-20,20);
		H_PhiT_phiT.setTitle("ΔΦ_{trento} vs Φ_{trento} ");
		H_PhiT_phiT.setTitleX(" Φ_{trento} ");
		H_PhiT_phiT.setTitleY("ΔΦ_{trento} [degrees]");
		
		H_dTheta_e = new H2F("H_dTheta_e",100,0,8,100,-5,5);
		H_dTheta_e.setTitle("Δθ vs E ");
		H_dTheta_e.setTitleX(" E [GeV] ");
		H_dTheta_e.setTitleY("Δθ [degrees]");
		
		H_dTheta_theta = new H2F("H_dTheta_theta",100,0,45,100,-5,5);
		H_dTheta_theta.setTitle("Δθ vs θ ");
		H_dTheta_theta.setTitleX(" θ [degrees] ");
		H_dTheta_theta.setTitleY("Δθ [degrees]");
		
		H_dTheta_phi = new H2F("H_dTheta_phi",100,-180,180,100,-5,5);
		H_dTheta_phi.setTitle("Δθ vs Φ ");
		H_dTheta_phi.setTitleX(" Φ [degrees] ");
		H_dTheta_phi.setTitleY("Δθ [degrees]");
		
		H_dPhi_e = new H2F("H_dPhi_e",100,0,8,100,-5,5);
		H_dPhi_e.setTitle("ΔΦ vs E ");
		H_dPhi_e.setTitleX(" E [GeV] ");
		H_dPhi_e.setTitleY("ΔΦ [degrees]");
		
		H_dPhi_theta = new H2F("H_dPhi_theta",100,0,45,100,-5,5);
		H_dPhi_theta.setTitle("ΔΦ vs Θ ");
		H_dPhi_theta.setTitleX(" θ [degrees] ");
		H_dPhi_theta.setTitleY("ΔΦ [degrees]");
		
		H_dPhi_phi = new H2F("H_dPhi_phi",100,-180,180,100,-5,5);
		H_dPhi_phi.setTitle("ΔΦ vs Φ  ");
		H_dPhi_phi.setTitleX(" Φ [degrees] ");
		H_dPhi_phi.setTitleY("ΔΦ [degrees]");
		
		
		//private static H2F H_xb_Q2_difference, H_W_Q2_difference, H_Y_Q2_difference;

		// Electrons :
		
		 H_e_dep  = new H1F("H_e_dep",100,-0.2,0.2);
		 H_e_dep.setTitle(" electron Δp/p  ");
		 H_e_dep.setTitleX("Δp/p");
		
		 
		 H_e_dp_p = new H2F("H_e_dp_p",100,0,10,100,-0.2,0.2);
		 H_e_dp_p.setTitle("electron Δp/p vs p ");
		 H_e_dp_p.setTitleX(" p [GeV] ");
		 H_e_dp_p.setTitleY("Δp/p");
			
		 H_e_dp_theta = new H2F("H_e_dp_theta",100,0,40,100,-0.2,0.2);
		 H_e_dp_theta.setTitle("electron Δp/p vs θ ");
		 H_e_dp_theta.setTitleX(" θ [degrees] ");
		 H_e_dp_theta.setTitleY("Δp/p");
		 
		 H_e_dp_phi = new H2F("H_e_dp_phi",100,-180,180,100,-0.2,0.2);
		 H_e_dp_phi.setTitle("electron Δp/p vs Φ ");
		 H_e_dp_phi.setTitleX(" Φ [degrees] ");
		 H_e_dp_phi.setTitleY("Δp/p");
		 
		 H_e_dtheta  = new H1F("H_e_dtheta",100,-1.0,1.0);
		 H_e_dtheta.setTitle(" electron Δθ  ");
		 H_e_dtheta.setTitleX(" Δθ  ");
		
		 
		 H_e_dtheta_p = new H2F("H_e_dtheta_p",100,0,10,100,-1.5,1.5);
		 H_e_dtheta_p.setTitle("electron Δθ vs p ");
		 H_e_dtheta_p.setTitleX(" p [GeV] ");
		 H_e_dtheta_p.setTitleY("Δθ [degrees]");
		 
		 H_e_dtheta_theta = new H2F("H_e_dtheta_theta",100,0,40,100,-1.5,1.5);
		 H_e_dtheta_theta.setTitle("electron Δθ vs θ ");
		 H_e_dtheta_theta.setTitleX(" θ [degrees] ");
		 H_e_dtheta_theta.setTitleY("Δθ [degrees]");
		 
		 H_e_dtheta_phi = new H2F("H_e_dtheta_phi",100,-180,180,100,-1.5,1.5);
		 H_e_dtheta_phi.setTitle("electron Δθ vs Φ ");
		 H_e_dtheta_phi.setTitleX(" Φ [degrees] ");
		 H_e_dtheta_phi.setTitleY("Δθ [degrees]");
		 
		 
	     H_e_dphi = new H1F("H_e_dphi", 100, -3,3);
	     H_e_dphi.setTitle(" electron ΔΦ  ");
	     H_e_dphi.setTitleX(" ΔΦ  ");
		 
	     H_e_dphi_p = new H2F("H_e_dphi_p",100,0,10,100,-3,3);
	     H_e_dphi_p.setTitle("electron ΔΦ  vs  p");
	     H_e_dphi_p.setTitleX(" p [GeV] ");
	     H_e_dphi_p.setTitleY("ΔΦ [degrees]");
	     
	     H_e_dphi_theta = new H2F("H_e_dphi_theta",100,0,40,100,-3,3);
	     H_e_dphi_theta.setTitle("electron ΔΦ  vs θ ");
	     H_e_dphi_theta.setTitleX(" θ [degrees] ");
	     H_e_dphi_theta.setTitleY("ΔΦ [degrees]");
	     
	     H_e_dphi_phi = new H2F("H_e_dphi_phi",100,-180,180,100,-3,3);
	     H_e_dphi_phi.setTitle("electron ΔΦ vs Φ ");
	     H_e_dphi_phi.setTitleX(" Φ [degrees] ");
	     H_e_dphi_phi.setTitleY("ΔΦ [degrees]");
	     
	     
	     H_e_dQ2 = new H1F("H_e_dQ2", 100, -0.3,0.3);
	     
	     H_e_dQ2_Q2 = new H2F("H_e_dQ2_Q2",100,0,10,100,-0.4,0.4);
	     H_e_dQ2_Q2.setTitle("electron ΔQ^2 vs Q^2 ");
	     H_e_dQ2_Q2.setTitleX(" Q^2 [GeV^2] ");
	     H_e_dQ2_Q2.setTitleY("ΔQ^2 [GeV^2]");
	     
	     H_e_dQ2_xB = new H2F("H_e_dQ2_xB",100,0,1,100,-0.4,0.4);
	     H_e_dQ2_xB.setTitle("electron ΔQ^2 vs xB ");
	     H_e_dQ2_xB.setTitleX(" xB  ");
	     H_e_dQ2_xB.setTitleY("ΔQ^2 [GeV^2]");
	     
	     H_e_dQ2_W = new H2F("H_e_dQ2_W",100,1.5,4.5,100,-0.4,0.4);
	     H_e_dQ2_W.setTitle("electron ΔQ^2 vs W ");
	     H_e_dQ2_W.setTitleX(" W  ");
	     H_e_dQ2_W.setTitleY("ΔQ^2 [GeV^2]");
	     
	     H_e_dxB = new H1F("H_e_dxB", 100, -0.3,0.3);
	     H_e_dxB.setTitle("electron ΔxB ");
	     H_e_dxB.setTitleX(" Gen. xB - Rec. xB  ");
	     
	     H_e_dxB_Q2 = new H2F("H_e_dxB_Q2",100,0,10,100,-0.4,0.4);
	     H_e_dxB_Q2.setTitle("electron ΔxB vs  Q^2 ");
	     H_e_dxB_Q2.setTitleX(" Q^2 [GeV^2]  ");
	     H_e_dxB_Q2.setTitleY("ΔxB");
	     
	     H_e_dxB_xB = new H2F("H_e_dxB_xB",100,0,1,100,-0.4,0.4);
	     H_e_dxB_xB.setTitle("electron ΔxB vs  xB ");
	     H_e_dxB_xB.setTitleX(" xB  ");
	     H_e_dxB_xB.setTitleY("ΔxB");
	     
	     H_e_dxB_W = new H2F("H_e_dxB_W",100,1.5,4.5,100,-0.4,0.4);
	     H_e_dxB_W.setTitle("electron ΔxB vs W ");
	     H_e_dxB_W.setTitleX(" W  ");
	     H_e_dxB_W.setTitleY("ΔxB");
	     
	     H_e_dW = new H1F("H_e_dW", 100, -0.3,0.3);
	     H_e_dxB.setTitle("electron ΔW ");
	     H_e_dxB.setTitleX(" Gen. W - Rec. W  ");
	     
	     H_e_dW_Q2 = new H2F("H_e_dW_Q2",100,0,10,100,-0.4,0.4);
	     H_e_dW_Q2.setTitle("electron ΔW vs  Q^2 ");
	     H_e_dW_Q2.setTitleX(" Q^2 [GeV^2]  ");
	     H_e_dW_Q2.setTitleY("ΔW");
	     
	     H_e_dW_xB = new H2F("H_e_dW_xB",100,0,1,100,-0.4,0.4);
	     H_e_dW_xB.setTitle("electron ΔW vs xB ");
	     H_e_dW_xB.setTitleX(" xB  ");
	     H_e_dW_xB.setTitleY("ΔW");
	     
	     H_e_dW_W = new H2F("H_e_dW_W",100,1.5,4.5,100,-0.4,0.4);
	     H_e_dW_W.setTitle("electron ΔW vs W ");
	     H_e_dW_W.setTitleX(" W ");
	     H_e_dW_W.setTitleY("ΔW");
	    
       //Generated and Reconstructed electron 
	     
	 	Gen_Q2_xB = new H2F("Gen_Q2_xB",100,0,1,100,0,14); Rec_Q2_xB= new H2F("Rec_Q2_xB",100,0,1,100,0,14);
	 	Gen_Q2_xB.setTitle("Generated Q^2 vs xB ");Rec_Q2_xB.setTitle("Reconstructed Q^2 vs xB ");
	 	Gen_Q2_xB.setTitleX(" xB ");Rec_Q2_xB.setTitleX(" xB ");
	 	Gen_Q2_xB.setTitleY("Q^2 [GeV^2] ");Rec_Q2_xB.setTitleY("Q^2 [GeV^2] ");
	    
	 	
	 	Gen_Q2_W = new H2F("Gen_Q2_W",100,1,5,100,0,14);Rec_Q2_W=new H2F("Rec_Q2_W",100,1,5,100,0,14);
	 	Gen_Q2_W.setTitle("Generated Q^2 vs W ");Rec_Q2_W.setTitle("Reconstructed Q^2 vs W ");
	 	Gen_Q2_W.setTitleX(" W ");Rec_Q2_W.setTitleX(" W ");
	 	Gen_Q2_W.setTitleY("Q^2 [GeV^2] ");Rec_Q2_W.setTitleY("Q^2 [GeV^2] ");
	 	
	 	Gen_xB_W = new H2F("Gen_xB_W",100,1,5,100,0,1);Rec_xB_W=new H2F("Rec_xB_W",100,1,5,100,0,1);
	 	Gen_xB_W.setTitle("Generated xB vs W ");Rec_xB_W.setTitle("Reconstructed xB vs W ");
	 	Gen_xB_W.setTitleX(" W ");Rec_xB_W.setTitleX(" W ");
	 	Gen_xB_W.setTitleY("xB ");Rec_xB_W.setTitleY("xB ");
	 	
	    Gen_e_theta= new H2F("Gen_e_theta",100,0,50,100,0,10); Rec_e_theta= new H2F("Rec_e_theta",100,0,50,100,0,10);
	    Gen_e_theta.setTitle("Generated E vs θ ");Rec_e_theta.setTitle("Reconstructed E vs θ ");
	    Gen_e_theta.setTitleX(" θ [degrees] ");Rec_e_theta.setTitleX(" θ [degrees] ");
	    Gen_e_theta.setTitleY("E [GeV] ");Rec_e_theta.setTitleY("E [GeV] ");
	 	
	    Gen_e_phi= new H2F("Gen_e_phi",100,-180,180,100,0,10); Rec_e_phi= new H2F("Rec_e_phi",100,-180,180,100,0,10);
	    Gen_e_phi.setTitle("Generated E vs Φ");Rec_e_phi.setTitle("Reconstructed E vs Φ ");
	    Gen_e_phi.setTitleX(" Φ [degrees] ");Rec_e_phi.setTitleX(" Φ [degrees] ");
	    Gen_e_phi.setTitleY("E [GeV] ");Rec_e_phi.setTitleY("E [GeV] ");
	    
	    
	    Gen_theta_phi=new H2F("Gen_theta_phi",100,-180,180,100,0,50);Rec_theta_phi=new H2F("Rec_theta_phi",100,-180,180,100,0,50);
	    Gen_theta_phi.setTitle("Generated θ vs Φ");Rec_theta_phi.setTitle("Reconstructed  θ vs Φ ");
	    Gen_theta_phi.setTitleX(" Φ [degrees] ");Rec_theta_phi.setTitleX(" Φ [degrees] ");
	    Gen_theta_phi.setTitleY("θ [degrees] ");Rec_theta_phi.setTitleY("θ [degrees] ");
		 
	    
	   
	    
	    H_generated_p  = new H1F("H_generated_p",36,0,8);
	    H_generated_p.setTitle(" Generated Hadron Momentum distribution  ");
	    H_generated_p.setTitleX(" p [GeV/c]  ");
	    
	    H_reconstructed_p  = new H1F("H_reconstructed_p",36,0,8);
	    H_reconstructed_p.setTitle(" Reconstructed Hadron Momentum distribution  ");
	    H_reconstructed_p.setTitleX(" p [GeV/c]  ");
	    
	}


	private void Plots_ph() {


		H_ph_pxBC = new H1F("H_ph_pxBC","H_ph_pxBC",100,0,5);
		H_ph_pxBC.setTitle(" Px photon Before Cuts");
		H_ph_pxBC.setTitleX("Momentum px [GeV]");
		H_ph_pxBC.setFillColor(4);

		H_ph_pxAC = new H1F("H_ph_pxAC","H_ph_pxBAC",100,0,5);
		H_ph_pxAC.setTitle(" Px photon After Cuts");
		H_ph_pxAC.setTitleX("Momentum px [GeV]");
		H_ph_pxAC.setFillColor(5);

		H_ph_pyBC = new H1F("H_ph_pyBC","H_ph_pyBC",100,0,5);
		H_ph_pyBC.setTitle(" Py photon Before Cuts");
		H_ph_pyBC.setTitleX("Momentum py [GeV]");
		H_ph_pyBC.setFillColor(4);

		H_ph_pyAC = new H1F("H_ph_pyAC","H_ph_pyAC",100,0,5);
		H_ph_pyAC.setTitle(" Py photon After Cuts");
		H_ph_pyAC.setTitleX("Momentum py [GeV]");
		H_ph_pyAC.setFillColor(5);

		H_ph_pzBC = new H1F("H_ph_pzBC","H_ph_pzBC",100,0,5);
		H_ph_pzBC.setTitle(" Pz photon Before Cuts");
		H_ph_pzBC.setTitleX("Momentum pz [GeV]");
		H_ph_pzBC.setFillColor(4);

		H_ph_pzAC = new H1F("H_ph_pzAC","H_ph_pzAC",100,0,5);
		H_ph_pzAC.setTitle(" Pz photon After Cuts");
		H_ph_pzAC.setTitleX("Momentum pz [GeV]");
		H_ph_pzAC.setFillColor(5);


		H_ph_MomBC = new H1F("H_ph_MomBC"," H_ph_MomBC",100,0,5);
		H_ph_MomBC.setTitle(" Momentum photon Before Cuts");
		H_ph_MomBC.setTitleX("Momentum  [GeV]");
		H_ph_MomBC.setFillColor(4);

		H_ph_MomAC = new H1F("H_ph_MomAC"," H_ph_MomAC",100,0,5);
		H_ph_MomAC.setTitle(" Momentum photon After Cuts");
		H_ph_MomAC.setTitleX("Momentum [GeV]");
		H_ph_MomAC.setFillColor(5);

		H_ph_ThetaBC = new H1F("H_ph_ThetaBC"," H_ph_ThetaBC",100,0,50);
		H_ph_ThetaBC.setTitle(" Polar Angle (θ) photon Before Cuts");
		H_ph_ThetaBC.setTitleX("Polar Angle (θ) [deg]");
		H_ph_ThetaBC.setFillColor(4);

		H_ph_ThetaAC = new H1F("H_ph_ThetaAC"," H_ph_ThetaAC",100,0,50);
		H_ph_ThetaAC.setTitle(" Polar Angle (θ) photon After Cuts");
		H_ph_ThetaAC.setTitleX("Polar Angle (θ) [deg]");
		H_ph_ThetaAC.setFillColor(5);

		H_ph_PhiBC = new H1F("H_ph_PhiBC"," H_ph_PhiBC",100,-180,180);
		H_ph_PhiBC.setTitle(" Azimuthal Angle (φ) photon Before Cuts");
		H_ph_PhiBC.setTitleX("Azimuthal Angle (φ) [deg]");
		H_ph_PhiBC.setFillColor(4);

		H_ph_PhiAC = new H1F("H_ph_PhiAC"," H_ph_PhiAC",100,-180,180);
		H_ph_PhiAC.setTitle(" Azimuthal Angle (φ) photon After Cuts");
		H_ph_PhiAC.setTitleX("Azimuthal Angle (φ) [deg]");
		H_ph_PhiAC.setFillColor(5);
		
		
		
	

	}
	// Rimuovere tutti i doppioni
	// Avere un plot del vertice 

	public void SetUnSkim_el(ParticleREC electron,double Q2, double xB, double W_var, double y_var) {
		//Add vertex somewhere. 
		//System.out.println("Vertex is "+ electron.vx() +" , " + electron.vx()+ " , "+ electron.vy());
		H_e_vzBC.fill(electron.vz());
		H_e_vx_vyBC.fill(electron.vx(),electron.vy());
		if (electron.getDetectorHits(DC1)!=null) {
		H_DC_xyNCut.fill(electron.getDetectorHits(DC1).x(), electron.getDetectorHits(DC1).y());}
		if ( electron.getDetectorHits(PCALID)!= null ){
		H_lu_eBC.fill(electron.getDetectorHits(PCALID).x());
		H_lv_eBC.fill(electron.getDetectorHits(PCALID).y());
		H_lw_eBC.fill(electron.getDetectorHits(PCALID).z());
		}
		
		
		//H_DC_xyNCut.fill(electron.getDetectorHits(DC1).get(0).x(), electron.getDetectorHits(DC1).get(0).y());
		//H_lu_eBC.fill(electron.getDetectorHits(PCALID).get(0).x());
		//H_lv_eBC.fill(electron.getDetectorHits(PCALID).get(0).y());
		//H_lw_eBC.fill(electron.getDetectorHits(PCALID).get(0).z());
		// Sample fraction should be redon looking at Example from GAGIK 
		// H_ESampl_ECalBC.fill(electron.p(),electron.e_samplFrac);
		// System.out.println(" Electron angle and momentum " + Math.toDegrees(electron.theta())+ " / " + electron.p());
		//System.out.println(" Variables Q " + Q2+ " y_var " + y_var+ " W_var "+W_var );
		H_e_ThetaBC.fill(Math.toDegrees(electron.theta()));
		H_e_MomBC.fill(electron.p());
		H_e_PhiBC.fill(Math.toDegrees(electron.phi()));
		H_e_y_Q2BC.fill(y_var,Q2);
		H_e_W_Q2BC.fill(W_var,Q2);
		H_e_xB_Q2BC.fill(xB, Q2);	
		H_e_WBC.fill(W_var);
		

		// if(electron.DC1_cut==false) H_DC_xyNCut.fill(electron.e_DC_X_1st_NS, electron.e_DC_Y_1st_NS);		
	}

	public void PrintUnSkim_el(String time) {
		new File(s1+"/Plots").mkdir();
		String PathFolder = s1+"/Plots/"+time;
		File dir = new File(PathFolder);
		dir.mkdir();
		PathFolder=s1+"/Plots/"+time+"/el_Before_Cut";
		new File (PathFolder).mkdir();
		//String PathFolder = s1+"/Plots/el_Before_Cut/"+time;
		System.out.println("I am Plotting Elctron Histogram");		
		EmbeddedCanvas electrons = new EmbeddedCanvas();
		electrons.setSize(1600,1000);
		electrons.divide(3,2);
		electrons.setAxisTitleSize(22);
		electrons.setAxisFontSize(22);
		electrons.setTitleSize(22);
		electrons.cd(0);electrons.draw(H_e_vzBC);
		electrons.cd(1);electrons.draw(H_e_vx_vyBC);
		electrons.cd(2);electrons.draw(H_lu_eBC);
		electrons.cd(3);electrons.draw(H_lv_eBC);
		electrons.cd(4);electrons.draw(H_lw_eBC);
		electrons.cd(5);electrons.draw(H_ESampl_ECalBC);
		electrons.cd(6);electrons.draw(H_e_WBC);
		electrons.cd(7);electrons.draw(H_DC_xyNCut);
		String strg0 = String.format("%s/Electrons.png",PathFolder);
		System.out.println("Saving plots in "+ PathFolder);
		electrons.save(strg0);		

		EmbeddedCanvas e_Kin = new EmbeddedCanvas();
		e_Kin.setSize(1600,1000);
		e_Kin.divide(3,2);
		e_Kin.setAxisTitleSize(28);
		e_Kin.setAxisFontSize(28);
		e_Kin.setTitleSize(28);
		e_Kin.cd(0);e_Kin.draw(H_e_ThetaBC);
		e_Kin.cd(1);e_Kin.draw(H_e_MomBC);	
		e_Kin.cd(2);e_Kin.draw(H_e_PhiBC);
		e_Kin.cd(3);e_Kin.draw(H_e_xB_Q2BC);
		e_Kin.cd(4);e_Kin.draw(H_e_y_Q2BC);
		e_Kin.cd(5);e_Kin.draw(H_e_W_Q2BC);
		String strg = String.format("%s/ele_Kienematic.png",PathFolder);
		System.out.println("Saving plots in "+ strg);
		e_Kin.save(strg);

		EmbeddedCanvas fiducial = new EmbeddedCanvas();
		fiducial.setSize(1200,1200);
		fiducial.divide(2,2);
		fiducial.setAxisTitleSize(28);
		fiducial.setAxisFontSize(28);
		fiducial.setTitleSize(28);	
		fiducial.cd(0);fiducial.draw(H_ECal_xy);
		fiducial.cd(1);fiducial.draw(H_ECal_xyNCut);
		fiducial.cd(2);fiducial.draw(H_DC_xy);
		fiducial.cd(3);fiducial.draw(H_DC_xyNCut);
		String strgF = String.format("%s/Electron_fiducial.png",PathFolder);
		System.out.println("Saving plots in "+ strgF);
		fiducial.save(strgF);		

	}

	public void SetSkim_el(ParticleREC electron,double Q2AC,double xBAC,double W_varAC,double y_varAC, double MC_Q2, double MC_xB, double MC_W, double MC_y) {
		H_e_vzAC.fill(electron.vz());
		H_e_vx_vyAC.fill(electron.vx(),electron.vy());
		H_e_ThetaAC.fill(Math.toDegrees(electron.theta()));
		H_e_MomAC.fill(electron.p());
		H_e_PhiAC.fill(Math.toDegrees(electron.phi()));
		H_e_y_Q2AC.fill(y_varAC, Q2AC);
		H_e_W_Q2AC.fill(W_varAC,Q2AC);
		HQ2AC.fill(Q2AC);
		H_e_xB_Q2AC.fill(xBAC,Q2AC);
		H_e_WAC.fill(W_varAC);
		H_e_xAC.fill(xBAC);
		H_e_yAC.fill(y_varAC);
		if (electron.getDetectorHits(PCALID)!=null && electron.getDetectorHits(DC1)!=null) {
		H_DC_xy.fill(electron.getDetectorHits(DC1).x(), electron.getDetectorHits(DC1).y());
		//System.out.println(electron.getDetectorHits(PCALID).x()+" "+electron.getDetectorHits(PCALID).y()+" "+ electron.getDetectorHits(PCALID).z());
		H_lu_eAC.fill(electron.getDetectorHits(PCALID).x());
		H_lv_eAC.fill(electron.getDetectorHits(PCALID).y());
		H_lw_eAC.fill(electron.getDetectorHits(PCALID).z());
}
		//H_DC_xy.fill(electron.getDetectorHits(DC1).get(0).x(), electron.getDetectorHits(DC1).get(0).y());
		//H_lu_eAC.fill(electron.getDetectorHits(PCALID).get(0).x());
		//H_lv_eAC.fill(electron.getDetectorHits(PCALID).get(0).y());
		//H_lw_eAC.fill(electron.getDetectorHits(PCALID).get(0).z());
		double distanceQ2=MC_Q2-Q2AC;
		H_Q2_difference.fill(MC_Q2,distanceQ2);
		H_Y_difference.fill(MC_y,(MC_y-y_varAC));
		H_xB_difference.fill(MC_xB,(MC_xB-xBAC));
		H_W_difference.fill(W_varAC,(MC_W-W_varAC));
		H_Diff_Q2.fill(distanceQ2);
		H_Diff_xB.fill((MC_xB-xBAC));
		H_Diff_Y.fill((MC_y-y_varAC));
		H_Diff_W.fill((MC_y-y_varAC));
		// }
	}

	public void PrintSkim_el(String time) {
		new File(s1+"/Plots").mkdir();
		String PathFolder = s1+"/Plots/"+time;
		File dir = new File(PathFolder);
		dir.mkdir();
		PathFolder=s1+"/Plots/"+time+"/el_After_Cut";
		new File (PathFolder).mkdir();
		//String PathFolder = s1+"/Plots/el_Before_Cut/"+time;
		TDirectory directory = new TDirectory();
		directory.mkdir("/electron/"); directory.cd("/electron/");			
		directory.addDataSet(H_e_vzAC);
		directory.addDataSet(H_e_ThetaAC);
		directory.addDataSet(H_e_MomAC);	
		directory.addDataSet(H_e_PhiAC);
		directory.addDataSet(HQ2AC);
		directory.addDataSet(H_e_yAC);
		directory.addDataSet(H_e_xAC);
		directory.addDataSet(H_e_WAC);
		directory.addDataSet(H_lu_eAC);
		directory.addDataSet(H_lv_eAC);
		directory.addDataSet(H_lw_eAC);
		directory.addDataSet(H_e_xB_Q2AC);
		directory.addDataSet(H_DC_xy);
		directory.writeFile(PathFolder+"/Histograms_Electrons.hipo");
		
		 	
		System.out.println("I am Plotting Elctron Histogram");		
		/*
			new File(s1+"/Plots").mkdir();
			new File (s1+"/Plots/el_After_Cut").mkdir();
			Date nowAC = new Date();
			SimpleDateFormat dateFormatAC = new SimpleDateFormat("MMMdd_hh_mm");

			String timeAC = dateFormatAC.format(nowAC);
			String PathFolderAC = s1+"/Plots/el_After_Cut/"+timeAC;
			File dirAC = new File(PathFolderAC);
			dirAC.mkdir();

			  
		 */
		System.out.println("I am Plotting Elctron Histogram After Cuts");

		EmbeddedCanvas electrons = new EmbeddedCanvas();
		electrons.setSize(1600,1000);
		electrons.divide(3,2);
		electrons.setAxisTitleSize(22);
		electrons.setAxisFontSize(22);
		electrons.setTitleSize(22);
		electrons.cd(0);electrons.draw(H_e_vzAC);
		electrons.cd(1);electrons.draw(H_e_vx_vyAC);
		electrons.cd(2);electrons.draw(H_e_WAC);
		electrons.cd(3);electrons.draw(H_DC_xy);
		electrons.cd(4);electrons.draw(H_lu_eAC);
		electrons.cd(5);electrons.draw(H_lv_eAC);
		String strg0 = String.format("%s/Electrons.png",PathFolder);
		System.out.println("Saving plots in "+ strg0);
		electrons.save(strg0);		

		EmbeddedCanvas e_Kin = new EmbeddedCanvas();
		e_Kin.setSize(1600,1000);
		e_Kin.divide(3,2);
		e_Kin.setAxisTitleSize(28);
		e_Kin.setAxisFontSize(28);
		e_Kin.setTitleSize(28);
		e_Kin.cd(0);e_Kin.draw(H_e_ThetaAC);
		e_Kin.cd(1);e_Kin.draw(H_e_MomAC);	
		e_Kin.cd(2);e_Kin.draw(H_e_PhiAC);
		e_Kin.cd(3);e_Kin.draw(H_e_xB_Q2AC);
		e_Kin.cd(4);e_Kin.draw(H_e_y_Q2AC);
		e_Kin.cd(5);e_Kin.draw(H_e_W_Q2AC);
		String strg = String.format("%s/ele_Kienematic.png",PathFolder);
		System.out.println("Saving plots in "+ strg);
		e_Kin.save(strg);

		EmbeddedCanvas Differences = new EmbeddedCanvas();
		Differences.setSize(1600,1000);
		Differences.setAxisTitleSize(28);
		Differences.setAxisFontSize(28);
		Differences.setTitleSize(28);
		Differences.divide(4, 2);
		Differences.cd(0);Differences.draw(H_Q2_difference);
		Differences.cd(1);Differences.draw(H_xB_difference);
		Differences.cd(2);Differences.draw(H_Y_difference);
		Differences.cd(3); Differences.draw(H_W_difference);
		Differences.cd(4);Differences.draw(H_Diff_Q2);
		Differences.cd(5);Differences.draw(H_Diff_xB);
		Differences.cd(6);Differences.draw(H_Diff_Y);
		Differences.cd(7); Differences.draw(H_Diff_W);
		String strg2 = String.format("%s/MC_Differences.png",PathFolder);
		System.out.println("Saving plots in "+ strg2);
		Differences.save(strg2);

		
	}

	public void SetSkim_ChargedHadron(double Mc_z, double Rec_z, double Mc_Pt, double Rec_pt, double Mc_theta, double Rec_Theta, double Mc_E, double Rec_E, double Mc_PhiClas, double Rec_PhiClas,double Mc_phiTrento, double Rec_phiTrento, double Mc_p, double Rec_p) {

		 H_generated_p.fill(Mc_p);
		 H_reconstructed_p.fill(Rec_p);
		//System.out.println(" Mc phi trento "+ Mc_phiTrento + " rec trento " + Rec_phiTrento);
		//phi trento etc. 
		H_z_difference.fill(Mc_z-Rec_z);
		H_dz_z.fill(Mc_z,Mc_z-Rec_z);
		H_dz_Pt.fill(Mc_Pt,Mc_z-Rec_z);
		H_dz_phi.fill(Mc_PhiClas,Mc_z-Rec_z);
		H_dz_theta.fill(Mc_theta,Mc_z-Rec_z);
		H_dz_phiT.fill(Mc_phiTrento,Mc_z-Rec_z);
		
		H_PT_difference.fill(Mc_Pt-Rec_pt);
		H_dPt_Pt.fill(Mc_Pt,Mc_Pt-Rec_pt);
		H_dPt_z.fill(Mc_z,Mc_Pt-Rec_pt);
		H_dPt_phi.fill(Mc_PhiClas,Mc_Pt-Rec_pt);
		H_dPt_theta.fill(Mc_theta,Mc_Pt-Rec_pt);
		H_dPt_phiT.fill(Mc_phiTrento,Mc_Pt-Rec_pt);
		
		H_PhiT_difference.fill(Mc_phiTrento-Rec_phiTrento);
		H_PhiT_Pt.fill(Mc_Pt,Mc_phiTrento-Rec_phiTrento);
		H_PhiT_z.fill(Mc_z,Mc_phiTrento-Rec_phiTrento);
		H_PhiT_phi.fill(Mc_PhiClas,Mc_phiTrento-Rec_phiTrento);
		H_PhiT_theta.fill(Mc_theta,Mc_phiTrento-Rec_phiTrento);
		H_PhiT_phiT.fill(Mc_phiTrento,Mc_phiTrento-Rec_phiTrento);
		
		H_Energy_difference.fill(Mc_E-Rec_E);
		H_dE_E.fill(Mc_E,Mc_E-Rec_E);
		H_dE_theta.fill(Mc_theta,Mc_E-Rec_E);
		H_dE_phi.fill(Mc_PhiClas,Mc_E-Rec_E);
		H_dE_phiT.fill(Mc_phiTrento, Mc_E-Rec_E);
			
		H_Theta_difference.fill(Mc_theta-Rec_Theta);
		H_dTheta_e.fill(Mc_E,Mc_theta-Rec_Theta);
		H_dTheta_theta.fill(Mc_theta,Mc_theta-Rec_Theta);
		H_dTheta_phi.fill(Mc_PhiClas,Mc_theta-Rec_Theta);
		
		
		H_Phi_difference.fill(Mc_PhiClas-Rec_PhiClas);
		H_dPhi_e.fill(Mc_E,Mc_PhiClas-Rec_PhiClas);
		H_dPhi_theta.fill(Mc_theta,Mc_PhiClas-Rec_PhiClas);
		H_dPhi_phi.fill(Mc_PhiClas,Mc_PhiClas-Rec_PhiClas);
		
		
	}

	public void ElectronGenerated(double mc_e_p,double mc_e_theta,double mc_e_phi,double mc_e_Q2,double mc_e_xB,double mc_e_W) {
		Gen_Q2_xB.fill(mc_e_xB, mc_e_Q2);
		Gen_Q2_W.fill(mc_e_W, mc_e_Q2);
		Gen_xB_W.fill(mc_e_W,mc_e_xB);
		Gen_xB_W.fill(mc_e_W,mc_e_xB);
		Gen_e_theta.fill(mc_e_theta,mc_e_p);
		Gen_e_phi.fill(mc_e_phi, mc_e_p);
		 Gen_theta_phi.fill(mc_e_phi, mc_e_theta);
	}
	public void ElectronResolutions(double mc_e_p,double mc_e_theta,double mc_e_phi,double mc_e_Q2,double mc_e_xB,double mc_e_W,ParticleREC electronRec,double rec_e_Q2,double rec_e_xB,double rec_e_W){
		
		//Gen_Q2_xB.fill(mc_e_xB, mc_e_Q2);  
	 	//Gen_Q2_W.fill(mc_e_W, mc_e_Q2); 
	 	//Gen_xB_W.fill(mc_e_W,mc_e_xB);
	   // Gen_e_theta.fill(mc_e_theta,mc_e_p);
	   // Gen_e_phi.fill(mc_e_phi, mc_e_p); 
	   // Gen_theta_phi.fill(mc_e_phi, mc_e_theta);
		 
	    
	    Rec_Q2_xB.fill(rec_e_xB, rec_e_Q2);
	    Rec_Q2_W.fill(rec_e_W, rec_e_Q2);
	    Rec_xB_W.fill(rec_e_W,rec_e_xB);
	    Rec_e_theta.fill(Math.toDegrees(electronRec.theta()), electronRec.p());
	    Rec_e_phi.fill(Math.toDegrees(electronRec.phi()),electronRec.p());
	    Rec_theta_phi.fill(Math.toDegrees(electronRec.phi()), Math.toDegrees(electronRec.theta()));
		 
	    
		H_e_dep.fill(mc_e_p-electronRec.p());
		H_e_dp_p.fill(mc_e_p,(mc_e_p-electronRec.p())/mc_e_p);
		H_e_dp_theta.fill(mc_e_theta,(mc_e_p-electronRec.p())/mc_e_p);
		H_e_dp_phi.fill(mc_e_phi,(mc_e_p-electronRec.p())/mc_e_p);
		
		//System.out.println(" valori elettrone angolo " + mc_e_theta+ " ricostruito " + Math.toDegrees(electronRec.theta()));
		
		H_e_dtheta.fill(mc_e_theta-Math.toDegrees(electronRec.theta()));
		H_e_dtheta_p.fill(mc_e_p,mc_e_theta-Math.toDegrees(electronRec.theta()));
		H_e_dtheta_theta.fill(mc_e_theta,mc_e_theta-Math.toDegrees(electronRec.theta()));
		H_e_dtheta_phi.fill(mc_e_phi,mc_e_theta-Math.toDegrees(electronRec.theta()));
		
		//System.out.println(" Valori phi mc "+ mc_e_phi + "  valori phi rec "+ Math.toDegrees(electronRec.phi()));
     	 H_e_dphi.fill(mc_e_phi-Math.toDegrees(electronRec.phi()));
     	 H_e_dphi_p.fill(mc_e_p,mc_e_phi-Math.toDegrees(electronRec.phi()));
     	 H_e_dphi_theta.fill(mc_e_theta,mc_e_phi-Math.toDegrees(electronRec.phi()));
     	 H_e_dphi_phi.fill(mc_e_phi,mc_e_phi-Math.toDegrees(electronRec.phi()));
     	
     	 
     	 //Q2
		H_e_dQ2.fill(mc_e_Q2-rec_e_Q2);
		H_e_dQ2_Q2.fill(mc_e_Q2,mc_e_Q2-rec_e_Q2);
		H_e_dQ2_xB.fill(mc_e_xB,mc_e_Q2-rec_e_Q2);
		H_e_dQ2_W.fill(mc_e_W,mc_e_Q2-rec_e_Q2);
		
		H_e_dxB.fill(mc_e_xB-rec_e_xB);
		H_e_dxB_Q2.fill(mc_e_Q2,mc_e_xB-rec_e_xB);
		H_e_dxB_xB.fill(mc_e_xB,mc_e_xB-rec_e_xB);
		H_e_dxB_W.fill(mc_e_W,mc_e_xB-rec_e_xB);
		
		H_e_dW.fill(mc_e_W-rec_e_W);
		H_e_dW_Q2.fill(mc_e_Q2,mc_e_W-rec_e_W);
		H_e_dW_xB.fill(mc_e_xB,mc_e_W-rec_e_W);
		H_e_dW_W.fill(mc_e_W,mc_e_W-rec_e_W);
		
		
	}
	
	public void Print_ChargedHadron (String time) {
		// Printing the histos into a file with histograms. 
		
		
		
			new File(s1+"/Plots").mkdir();
			String PathFolder = s1+"/Plots/"+time;
			File dir = new File(PathFolder);
			dir.mkdir();
			PathFolder=s1+"/Plots/"+time+"/Hadron";
			new File (PathFolder).mkdir();
			TDirectory directory = new TDirectory();
			TDirectory directoryel = new TDirectory();
			directory.mkdir("/hadron/"); directory.cd("/hadron/");	
			directory.addDataSet(H_dz_z);
			directory.addDataSet(H_dz_Pt);
			directory.addDataSet(H_dz_phi);
			directory.addDataSet(H_dz_phiT);
			directory.addDataSet(H_dz_theta);
			directory.addDataSet(H_dPt_Pt);
			directory.addDataSet(H_dPt_z);
			directory.addDataSet(H_dPt_phi);
			directory.addDataSet(H_dPt_theta);
			directory.addDataSet(H_dE_E);
			directory.addDataSet(H_dE_theta);
			directory.addDataSet(H_dE_phi);
			directory.addDataSet(H_z_difference);
			directory.addDataSet(H_PT_difference);
			directory.addDataSet(H_Theta_difference);
			directory.addDataSet(H_Energy_difference);
			directory.addDataSet(H_dz_phiT);
			directory.addDataSet(H_dPt_phiT);
			directory.addDataSet(H_dE_phiT);
			directory.addDataSet(H_Phi_difference);
			directory.addDataSet(H_PhiT_difference);
			directory.addDataSet(H_PhiT_Pt);
			directory.addDataSet(H_PhiT_z);
			directory.addDataSet(H_PhiT_phi);
			directory.addDataSet(H_PhiT_theta);
			directory.addDataSet(H_PhiT_phiT);
			directory.addDataSet(H_dTheta_e);
			directory.addDataSet(H_dTheta_theta);
			directory.addDataSet(H_dTheta_phi);
			directory.addDataSet(H_dPhi_e);
			directory.addDataSet(H_dPhi_theta);
			directory.addDataSet(H_dPhi_phi);
		   directory.addDataSet(H_generated_p);
		   directory.addDataSet(H_reconstructed_p);
			
			//directory.cd("..");
			directoryel.mkdir("/electron/"); directoryel.cd("/electron/");
			directoryel.addDataSet(H_e_dep);
			directoryel.addDataSet(H_e_dp_p);
			directoryel.addDataSet(H_e_dp_theta);
			directoryel.addDataSet(H_e_dp_phi);
			directoryel.addDataSet(H_e_dtheta);
			directoryel.addDataSet(H_e_dtheta_p);
			directoryel.addDataSet(H_e_dtheta_theta);
			directoryel.addDataSet(H_e_dtheta_phi);
			directoryel.addDataSet(H_e_dphi);
			directoryel.addDataSet(H_e_dphi_p);
			directoryel.addDataSet(H_e_dphi_theta);
			directoryel.addDataSet(H_e_dphi_phi);
			directoryel.addDataSet(H_e_dQ2);
			directoryel.addDataSet(H_e_dQ2_Q2);
			directoryel.addDataSet(H_e_dQ2_xB);
			directoryel.addDataSet(H_e_dQ2_W);
			directoryel.addDataSet(H_e_dxB);
			directoryel.addDataSet(H_e_dxB_Q2);
			directoryel.addDataSet(H_e_dxB_xB);
			directoryel.addDataSet(H_e_dxB_W);
			directoryel.addDataSet(H_e_dW);
			directoryel.addDataSet(H_e_dW_Q2);
			directoryel.addDataSet(H_e_dW_xB);
			directoryel.addDataSet(H_e_dW_W);
		    directoryel.addDataSet(Gen_Q2_xB);
		    directoryel.addDataSet(Rec_Q2_xB);
		    directoryel.addDataSet(Gen_Q2_W);
		    directoryel.addDataSet(Rec_Q2_W);
		    directoryel.addDataSet(Gen_xB_W);
		    directoryel.addDataSet(Rec_xB_W);
		    directoryel.addDataSet(Gen_e_theta);
		    directoryel.addDataSet(Rec_e_theta);
		    directoryel.addDataSet(Gen_e_phi);
		    directoryel.addDataSet(Rec_e_phi);
		    directoryel.addDataSet(Gen_theta_phi);
		    directoryel.addDataSet(Rec_theta_phi);
		    
		    			 
			
			// Add electron plots in a differentfolder "electron" 
			
			directory.writeFile(PathFolder+"/Histograms_Hadron.hipo");
			directoryel.writeFile(PathFolder+"/Histograms_Electron.hipo");
			System.out.println("I am Plotting Hadron  and Electron comparisons");
/*
			F1D fun = new F1D("f1d","[a]*gaus(x,[m],[s])",-0.05,0.05);
			fun.setParameter(0, 1000);
			fun.setParameter(1, 0);
			fun.setParameter(2, 0.01);
			
			fun.setLineColor(2);
			fun.setLineWidth(5);
			DataFitter.fit(fun, H_PT_difference, "Q");
			
			System.out.println(" ----FITTING RESULTS----- ");
			System.out.println(" Par 0 " +fun.getParameter(0));
			System.out.println(" Par 1 " +fun.getParameter(1));
			System.out.println(" Par 2 " +fun.getParameter(2));
			System.out.println(fun.getParameterEstimate());
			*/
			File File_txt = new File(PathFolder+"/Histos.txt"); 
			
				try {
					File_txt.createNewFile();
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
				if(File_txt.exists()) {
					System.out.println(" Output File Found" );
					try {
						PrintWriter outP = new PrintWriter(File_txt);			
						outP.println("Z  mean "+ H_z_difference.getMean()+ " RMS " + H_z_difference.getRMS());
						outP.println("PT mean "+ H_PT_difference.getMean()+ " RMS " + H_PT_difference.getRMS());
						outP.println("Theta mean "+ H_Theta_difference.getMean()+ " RMS " + H_Theta_difference.getRMS());
						outP.println("Energy mean "+ H_Energy_difference.getMean()+ " RMS " + H_Energy_difference.getRMS());
						outP.println("Phi mean "+ H_Phi_difference.getMean()+ " RMS " + H_Phi_difference.getRMS());
						outP.println("Phi Trento mean "+ H_PhiT_difference.getMean()+ " RMS " + H_PhiT_difference.getRMS());
					    outP.println(" ========= Electron : =====");
					    outP.println("El p mean "+ H_e_dep.getMean() + " RMS "+ H_e_dep.getRMS());
						outP.close();
					} catch (FileNotFoundException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			
			EmbeddedCanvas hadron = new EmbeddedCanvas();
			hadron.setSize(1600,1000);
			hadron.divide(3,2);
			hadron.setAxisTitleSize(22);
			hadron.setAxisFontSize(22);
			hadron.setTitleSize(22);
			hadron.cd(0);hadron.draw(H_z_difference);
			hadron.cd(1);hadron.draw(H_PT_difference);
			hadron.cd(2);hadron.draw(H_Theta_difference);
			hadron.cd(3);hadron.draw(H_Energy_difference);
			hadron.cd(4);hadron.draw(H_Phi_difference);
			hadron.cd(5);hadron.draw(H_PhiT_difference);
			String strg0 = String.format("%s/Hadrons.png",PathFolder);
			System.out.println("Saving plots in "+ strg0);
			hadron.save(strg0);		
			
	}
	
	public void SetUnSkim_ph(Photons ph) {
		//Add vertex somewhere. 

		for (int i=0; i<ph.Photons_px.size(); i++) {
			//check if the size is correct 
			//System.out.println("Index "+ i );
			H_ph_pxBC.fill(ph.Photons_px.get(i));
			H_ph_pyBC.fill(ph.Photons_py.get(i));
			H_ph_pzBC.fill(ph.Photons_pz.get(i));
			H_ph_MomBC.fill(ph.Photons_mom.get(i));
			H_ph_ThetaBC.fill(ph.Photons_theta.get(i));
			H_ph_PhiBC.fill(ph.Photons_phi.get(i));
		}

	}	

	public void SetSkim_ph(Photons ph) {
		//Add vertex somewhere. 
		for (int i=0; i<ph.Photons_pxAC.size(); i++) {
			//check if the size is correct 
			H_ph_pxAC.fill(ph.Photons_pxAC.get(i));
			H_ph_pyAC.fill(ph.Photons_pyAC.get(i));
			H_ph_pzAC.fill(ph.Photons_pzAC.get(i));
			H_ph_MomAC.fill(ph.Photons_momAC.get(i));
			H_ph_ThetaAC.fill(ph.Photons_thetaAC.get(i));
			H_ph_PhiAC.fill(ph.Photons_phiAC.get(i));
		}

	}	

	public void PrintUnSkim_ph(String time) {

		new File(s1+"/Plots").mkdir();
		String PathFolder = s1+"/Plots/"+time;
		File dir = new File(PathFolder);
		dir.mkdir();
		PathFolder=s1+"/Plots/"+time+"/ph_Before_Cut";
		new File (PathFolder).mkdir();

		 

		System.out.println("I am Plotting Photons Histogram");

		EmbeddedCanvas photons = new EmbeddedCanvas();
		photons.setSize(1600,1000);
		photons.divide(3,2);
		photons.setAxisTitleSize(22);
		photons.setAxisFontSize(22);
		photons.setTitleSize(22);
		photons.cd(0);photons.draw(H_ph_pxBC);
		photons.cd(1);photons.draw(H_ph_pyBC);
		photons.cd(2);photons.draw(H_ph_pzBC);
		photons.cd(3);photons.draw(H_ph_MomBC);
		photons.cd(4);photons.draw(H_ph_ThetaBC);
		photons.cd(5);photons.draw(H_ph_PhiBC);
		String strg0 = String.format("%s/Photons.png",PathFolder);
		System.out.println("Saving plots in "+ PathFolder);
		photons.save(strg0);		

	}

	public void PrintSkim_ph(String time) {

		new File(s1+"/Plots").mkdir();
		String PathFolder = s1+"/Plots/"+time;
		File dir = new File(PathFolder);
		dir.mkdir();
		PathFolder=s1+"/Plots/"+time+"/ph_After_Cut";
		new File (PathFolder).mkdir();
		 

		System.out.println("I am Plotting Photons Histogram After Cuts");
		EmbeddedCanvas photons = new EmbeddedCanvas();
		photons.setSize(1600,1000);
		photons.divide(3,2);
		photons.setAxisTitleSize(22);
		photons.setAxisFontSize(22);
		photons.setTitleSize(22);
		photons.cd(0);photons.draw(H_ph_pxAC);
		photons.cd(1);photons.draw(H_ph_pyAC);
		photons.cd(2);photons.draw(H_ph_pzAC);
		photons.cd(3);photons.draw(H_ph_MomAC);
		photons.cd(4);photons.draw(H_ph_ThetaAC);
		photons.cd(5);photons.draw(H_ph_PhiAC);
		String strg0 = String.format("%s/Photons.png",PathFolder);
		System.out.println("Saving plots in "+ PathFolder);
		photons.save(strg0);


		TDirectory directory = new TDirectory();
		directory.mkdir("/photon/"); directory.cd("/photon/");			
		directory.addDataSet(H_ph_pxAC);
		directory.addDataSet(H_ph_pyAC);
		directory.addDataSet(H_ph_pzAC);	
		directory.addDataSet(H_ph_MomAC);
		directory.addDataSet(H_ph_ThetaAC);
		directory.addDataSet(H_ph_PhiAC);
		directory.writeFile(PathFolder+"/Histograms_Photons.hipo");


	}

}
