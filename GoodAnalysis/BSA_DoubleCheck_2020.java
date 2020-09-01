package GoodAnalysis;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Scanner;

import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.math.F1D;
import org.jlab.groot.math.Func1D;
import org.jlab.groot.ui.TGCanvas;

import TheoreticalCurves.Grids;

public class BSA_DoubleCheck_2020 {

	private static final String Outputdir = null;
	public BSA_DoubleCheck_2020() {
		// TODO Auto-generated constructor stub
	}

	  static double[][][][] TotalPions ;
      static double[][][][] TPion_helicity ;
      static GraphErrors[][][][] BSA_Plots ;
      static H1F[][][][] BSA_Histos ;
      static  F1D[][][][] Fit_Functions ;
            
      
      static double[][][][][] HelicityP ;
      static double[][][][][] HelicityM ;
      
      static double[][][][] BSA_Values ;
      static double[][][][] BSA_Errors ;
      static double[][][][] Str_F_Values ;
      static double[][][][] Str_F_errors ;
      static double[][][][] RChiSquared;
      
       static double[][][][] Q2_Values ;
      static double[][][][] xB_Values ;
      static double[][][][] W_Values ;
      static double[][][][] y_Values ;
      static double[][][][] eps_Values ;
      static double[][][][] Fattore_Values ;
      
      static double[][][][] z_Values ;
	  static double[][][][] Pt_Values;
	  
	  static H1F[][][][] Dis_Q2 ;
      static H1F[][][][] Dis_xB ;
      static H1F[][][][] Dis_W ;
      static H1F[][][][] Dis_y ;
      static H1F[][][][] Dis_z;
      static H1F[][][][] Dis_Pt;
	  

	    
	   // **** BINS ****
	   static int Phi_Bins=12;
		static int Z_Bins = 1; // They get updated by the histograms themself
		//int Z_Bins=1;
		static int PT_Bins=1; // They get updated by the histograms themself
		static int Q2_Bins = 4;
		static int XB_Bins =7;
		
		
	public static void main(String[] args) throws IOException {
		
		
		// TODO Auto-generated method stub
		
		
		int index=5; // 0 : 1D Q distribution | 1: 1D x distribution | 2: 1D z dist | 3 : 1D Pt dit | 4 :2D Q,x dist | 5: 4D |
		
		String FileData = "ID_"+index+"/";
		String NewPath = "/Users/gangelini/work/Analysis_Results/Latest/"+FileData;
		
		
		
	       String Outputdir = NewPath+"/Results_BSA/";
			          
		 // **** BINS ****
		    Phi_Bins=12;
		
			   if(index==0) { Q2_Bins=13 ; XB_Bins=1; Z_Bins =1 ; PT_Bins=1;} // FileData = "New_ID0";}
			   if(index==1) { Q2_Bins=1 ; XB_Bins=11; Z_Bins =1 ; PT_Bins=1;} // FileData = "New_ID1";}
			   if(index==2) { Q2_Bins=1 ; XB_Bins=1; Z_Bins =12 ; PT_Bins=1;} // FileData = "New_ID2";}
			   if(index==3) { Q2_Bins=1 ; XB_Bins=1; Z_Bins =1 ; PT_Bins=14;} // FileData = "New_ID3";}
			   if(index==4) { Q2_Bins=4 ; XB_Bins=7; Z_Bins =1 ; PT_Bins=1;}  //FileData = "New_ID4";}
			   if(index==5) { Q2_Bins=4 ; XB_Bins=7; Z_Bins =9 ; PT_Bins=8;}  //FileData = "New_ID5";}
			   
			   		   

		// **** CANVAS ****
	 
		    TDirectory SidisData_PiP = new TDirectory();
		    TDirectory countsData_PiP = new TDirectory();

		        
		        double Radconv = 0.01745;
		     TotalPions = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		       TPion_helicity = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		        BSA_Plots = new GraphErrors[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		            BSA_Histos = new  H1F[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		            Fit_Functions = new  F1D[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		              
		            HelicityP = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins][Phi_Bins];
				     HelicityM = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins][Phi_Bins];
				     
		      BSA_Values = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		       BSA_Errors = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		        Str_F_Values = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		        Str_F_errors = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		        RChiSquared= new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		        
		        
		        Q2_Values = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		        xB_Values = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		        W_Values = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		        y_Values = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		        eps_Values = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		       Fattore_Values = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
		        
		       z_Values = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
			      Pt_Values = new double[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
			      
			      
			      Dis_Q2 = new H1F[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
			      Dis_xB = new H1F[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
			       Dis_W = new H1F[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
			       Dis_y = new H1F[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
			       Dis_z= new H1F[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];
			       Dis_Pt= new H1F[Q2_Bins][XB_Bins][Z_Bins][PT_Bins];

for(int a=0 ; a<Q2_Bins ; a++) {
		      for (int b=0; b<XB_Bins; b++) {    
		       // Reading the file in the Q,x folder:               
		      countsData_PiP.readFile(NewPath+"/MultiBins/Q2_" + a+"/Xb_" + b +"/HistogramsCounts.hipo");
		      SidisData_PiP.readFile(NewPath+"/MultiBins/Q2_" + a+"/Xb_" + b +"/Histograms_SIDIS.hipo");
		      
				H2F PionCounts_PiP =  (H2F) countsData_PiP.getObject("/main/PionCounts");
				H1F Dist_Q2 = (H1F) countsData_PiP.getObject("/main/Dist_Q2"); Dist_Q2.setTitleX("Q2"); Dist_Q2.setFillColor(2);Dist_Q2.setOptStat("1111");
				H1F Dist_xB = (H1F) countsData_PiP.getObject("/main/Dist_xB"); Dist_xB.setTitleX("xB"); Dist_xB.setFillColor(3);Dist_xB.setOptStat("1111");
				H1F Dist_W = (H1F) countsData_PiP.getObject("/main/Dist_W"); Dist_W.setTitleX(" W "); Dist_W.setFillColor(4);Dist_W.setOptStat("1111");
				H1F Dist_y = (H1F) countsData_PiP.getObject("/main/Dist_y");Dist_y.setTitleX(" y "); Dist_y.setFillColor(5);Dist_y.setOptStat("1111");
				H1F Dist_eps = (H1F) countsData_PiP.getObject("/main/Dist_eps");Dist_eps.setTitleX(" epsilon"); Dist_eps.setFillColor(6);Dist_eps.setOptStat("1111");
				double fattore = Math.sqrt(2*Dist_eps.getMean()*(1-Dist_eps.getMean()));
				 System.out.println(" EPS MIN "+ Dist_eps.getMean());
				 System.out.println(" FATTORE MIN "+ fattore);
				System.out.println(" Q2 mean " + Dist_Q2.getMean() + " bin Q2 "+a + " bin XB "+b);
				
				//Definitions of Z and PT bins using Histograms
				 Z_Bins=PionCounts_PiP.getXAxis().getNBins() ; PT_Bins=PionCounts_PiP.getYAxis().getNBins();
				System.out.println(" Z bins are " + Z_Bins + " PT bins are "+ PT_Bins);
		
				TGCanvas Test = new TGCanvas("canvas_BSA","canvas_BSA",1600,1400); // Create Canvas with (1,1) divisions
				Test.divide(2, 2); Test.cd(0);Test.draw(PionCounts_PiP.getSlicesX().get(0));
					int valoritotali=0;
				for (int k=0 ; k<PionCounts_PiP.getSlicesX().get(0).getxAxis().getNBins();k++) {
					valoritotali+=PionCounts_PiP.getSlicesX().get(0).getBinContent(k);
				}
				System.out.println(" ###### Valori totali sono " + valoritotali);
				// Running over Z_Bins and PT_BINS;
				        
				
				
				
						    for (int i=0; i<Z_Bins; i++) {
						    	 for ( int j=0; j< PT_Bins; j++) { 
						    		System.out.println(" Z bins " + i + " PT bins" + j );
						            H1F BSA_PiP  = (H1F) SidisData_PiP.getObject("/Z"+i+"PT"+j+"/PhiBSA");		
						            H1F HelicityPositive = (H1F)SidisData_PiP.getObject("/Z"+i+"PT"+j+"/HelicityPositive");
						            H1F HelicityNegative = (H1F)SidisData_PiP.getObject("/Z"+i+"PT"+j+"/HelicityNegative");
						            H1F BSAPi= new H1F ("BSAPi",Phi_Bins,0,6.283); 
						        	H1F Dist_z = (H1F) SidisData_PiP.getObject("/Z"+i+"PT"+j+"/Dist_z"); Dist_z.setTitleX(" z "); Dist_z.setFillColor(7);Dist_z.setOptStat("1111");
						        	H1F Dist_Pt = (H1F) SidisData_PiP.getObject("/Z"+i+"PT"+j+"/Dist_Pt");  Dist_Pt.setTitleX(" Pt "); Dist_Pt.setFillColor(8);Dist_Pt.setOptStat("1111");
						        	
						        	// = Filling the Distributions 
						        	Fattore_Values[a][b][i][j]=fattore;
									Q2_Values[a][b][i][j]=Dist_Q2.getMean();
									xB_Values[a][b][i][j]=Dist_xB.getMean();
									W_Values[a][b][i][j]=Dist_W.getMean();
									y_Values[a][b][i][j]=Dist_y.getMean();
									eps_Values[a][b][i][j]=Dist_eps.getMean();
						        	z_Values[a][b][i][j]= Dist_z.getMean();
						        	Pt_Values[a][b][i][j]=Dist_Pt.getMean();
						        	
						      	    Dis_Q2[a][b][i][j] =Dist_Q2 ;
						            Dis_xB[a][b][i][j] =Dist_xB ;
						            Dis_W[a][b][i][j] =Dist_W ;
						            Dis_y[a][b][i][j] =Dist_y ;
						            Dis_z[a][b][i][j] =Dist_z ;
						            Dis_Pt[a][b][i][j] =Dist_Pt ;
						            
      
						        	System.out.println( " ----> z mean value " + Dist_z.getMean() + " z variable loop " + i + " pt variable loop " + j);
						        	System.out.println( " ----> Pt mean value " + Dist_Pt.getMean());
						        	
						        	double DataBCount_PiP= PionCounts_PiP.getBinContent(i, j);
						        	
						            F1D fun = new F1D("[a]*sin(x)","[p0]*sin(x)", 0,6.283);
						            fun.setParameter(0, 0.02);	
						            fun.setLineColor(2);
						            fun.setLineWidth(4);

     						     GraphErrors BSAPlot = new GraphErrors();
     						      BSAPlot.setTitleX("Angle Φ [rad]");
     						      BSAPlot.setTitleY(" BSA ");
     						     	int totalpion =0;
						   // Loop over Angle phi  
     						     /*	
     						     for(int k=0; k< BSA_PiP.getxAxis().getNBins();k++) {
     						    	 if(k==0) {
     						    	 HelicityM[a][b][i][j][11]=HelicityNegative.getBinContent(k);
     							    HelicityP[a][b][i][j][11]=HelicityPositive.getBinContent(k);
     						    	 }
     						    	 if(k==1) {
         						    	 HelicityM[a][b][i][j][10]=HelicityNegative.getBinContent(k);
         							    HelicityP[a][b][i][j][10]=HelicityPositive.getBinContent(k);
         						    	 }
     						    	 
     						    	if(k==2) {
        						    	 HelicityM[a][b][i][j][9]=HelicityNegative.getBinContent(k);
        							    HelicityP[a][b][i][j][9]=HelicityPositive.getBinContent(k);
        						    	 }
     						    	if(k==3) {
       						    	 HelicityM[a][b][i][j][8]=HelicityNegative.getBinContent(k);
       							    HelicityP[a][b][i][j][8]=HelicityPositive.getBinContent(k);
       						    	 }
     						    	if(k==4) {
          						    	 HelicityM[a][b][i][j][7]=HelicityNegative.getBinContent(k);
          							    HelicityP[a][b][i][j][7]=HelicityPositive.getBinContent(k);
          						    	 }
     						    	if(k==5) {
         						    	 HelicityM[a][b][i][j][6]=HelicityNegative.getBinContent(k);
         							    HelicityP[a][b][i][j][6]=HelicityPositive.getBinContent(k);
         						    	 }
     						    	if(k==6) {
        						    	 HelicityM[a][b][i][j][5]=HelicityNegative.getBinContent(k);
        							    HelicityP[a][b][i][j][5]=HelicityPositive.getBinContent(k);
        						    	 }
     						    	if(k==7) {
       						    	 HelicityM[a][b][i][j][4]=HelicityNegative.getBinContent(k);
       							    HelicityP[a][b][i][j][4]=HelicityPositive.getBinContent(k);
       						    	 }
     						    	if(k==8) {
       						    	 HelicityM[a][b][i][j][3]=HelicityNegative.getBinContent(k);
       							    HelicityP[a][b][i][j][3]=HelicityPositive.getBinContent(k);
       						    	 }
     						    	if(k==9) {
          						    	 HelicityM[a][b][i][j][2]=HelicityNegative.getBinContent(k);
          							    HelicityP[a][b][i][j][2]=HelicityPositive.getBinContent(k);
          						    	 }
     						    	if(k==10) {
         						    	 HelicityM[a][b][i][j][1]=HelicityNegative.getBinContent(k);
         							    HelicityP[a][b][i][j][1]=HelicityPositive.getBinContent(k);
         						    	 }
     						    	if(k==11) {
        						    	 HelicityM[a][b][i][j][0]=HelicityNegative.getBinContent(k);
        							    HelicityP[a][b][i][j][0]=HelicityPositive.getBinContent(k);
        						    	 }
     						     }
     						   */  	
     						 if(Phi_Bins != BSA_PiP.getxAxis().getNBins()) { System.out.println(" Attention Bins Defined in code are different from bins in data (xAxis of BSA "); return;}
						   for(int k =0 ; k< BSA_PiP.getxAxis().getNBins();k++)
						     {
							   /*
							   double BSAnum =  HelicityP[a][b][i][j][k] - HelicityM[a][b][i][j][k];
						    	double BSAden =  HelicityP[a][b][i][j][k] + HelicityM[a][b][i][j][k];
							   double BSACalculation =(1.1587486)* (BSAnum/BSAden);
							   */
							   HelicityP[a][b][i][j][k]=HelicityPositive.getBinContent(k);
							   HelicityM[a][b][i][j][k]=HelicityNegative.getBinContent(k);
						    double BSAnum =  HelicityPositive.getBinContent(k) - HelicityNegative.getBinContent(k);
						    double BSAden =  HelicityPositive.getBinContent(k) + HelicityNegative.getBinContent(k);
						    totalpion+=BSAden;
						    double BSACalculation =-(1.1587486)* (BSAnum/BSAden);  	
						    	BSAPi.setBinContent(k, BSACalculation); 
						       if(HelicityNegative.getBinContent(k)==0&&HelicityPositive.getBinContent(k)==0) {
						    	System.out.println(" This entry is 0 , k: " + k);
						    	BSACalculation=0;
						    }
						      // System.out.println("DataBCount for bin "+ j +" is "+ DataBCount_PiP);
						       //System.out.println("Total Pion "+totalpion);
						    	 double Denominator = (HelicityPositive.getBinContent(k)+HelicityNegative.getBinContent(k))*(HelicityPositive.getBinContent(k)+HelicityNegative.getBinContent(k));
						    	 double part1 = (1.15874*HelicityNegative.getBinContent(k)*2*Math.sqrt(HelicityPositive.getBinContent(k)))/Denominator;
						    	 double part2 = (1.15874*HelicityPositive.getBinContent(k)*2*Math.sqrt(HelicityPositive.getBinContent(k)))/Denominator;
						    	 double errore = Math.sqrt(part1*part1+part2*part2);
						    	 if(HelicityNegative.getBinContent(k)==0 && HelicityPositive.getBinContent(k)==0) errore=0;
						    BSAPi.setBinError(k, errore); // set the error on the BSA plot
						    BSAPlot.addPoint(BSA_PiP.getxAxis().getBinCenter(k)*0.01745,  BSACalculation,(8.66)*0.01745, errore); 						    	
						     }
						   //Filling Pion coutns
						   TPion_helicity[a][b][i][j]=totalpion ;
						   TotalPions[a][b][i][j]=DataBCount_PiP ;
						 
						    
				           //Fitting
						   DataFitter.fit(fun, BSAPlot, "QR"); 
						   BSA_Plots[a][b][i][j]= BSAPlot; 
						  System.out.println(" ----====== >>> Valore massimo" + BSA_Plots[a][b][i][j].getMax() + " valore a " + a );
						   BSA_Histos[a][b][i][j]= BSAPi; 
						   RChiSquared[a][b][i][j] = fun.getChiSquare()/11;
						   double Computedp0= fun.parameter(0).value();
						   double err = fun.parameter(0).error();
						   System.out.println(" ===== RISULTATI ======="  );

						   System.out.println(" Valore Parametero " +Computedp0 + " valore errore "+err  );
						   System.out.println(" ==============================="  );
						   Fit_Functions[a][b][i][j]=fun;
						   BSA_Values[a][b][i][j]=Computedp0;
						   BSA_Errors[a][b][i][j]=err;
						   Str_F_Values[a][b][i][j]=Computedp0 /fattore;
						   Str_F_errors[a][b][i][j]=err /fattore;
						   System.out.println("FATTORE " + fattore);
						    } // END PT
		        	} // END Z 
						    
						    
						    
						 
						    
						    /*
						    TGCanvas canvas_BSA_z = new TGCanvas("canvas_BSA","canvas_BSA",1600,1400); // Create Canvas with (1,1) divisions
						    int XCanvas=((int) Math.floor(Math.sqrt(BSAPlots.size())));
						    int YCanvas=(int) Math.floor(Math.sqrt(BSAPlots.size()));
						
						    canvas_BSA_z.divide(XCanvas+1, YCanvas+1);
							if(a==Print_Bin_Q && b==Print_Bin_X) {
								for(int i=0;  i<BSAPlots.size(); i++) {		
								canvas_BSA_z.cd(i);
								canvas_BSA_z.draw(BSAPlots.get(i));
								c_Bin_Distri.cd(5);c_Bin_Distri.draw(BSAPlots.get(i));
								//canvas_BSA_z.draw(BSAHistos.get(i));					
								}
							}
						*/
						
		       }//for Xb
		}// for Q2
if(index==0) Print_outoput_ID0(Outputdir);
if(index==1) Print_outoput_ID1(Outputdir);   
if(index==2) Print_outoput_ID2(Outputdir); 
if(index==3) Print_outoput_ID3(Outputdir); 
if(index==4) Print_outoput_ID4(Outputdir);
if(index==5) Print_outoput_ID5(Outputdir);
	}


	private static void Print_outoput_ID0(String Outputdir) throws FileNotFoundException {
		TGCanvas Canvas_BSA = new TGCanvas("canvas_BSA","canvas_BSA",1600,1400); // Create Canvas with (1,1) divisions
	       TGCanvas Canvas_Dis_Q2 = new TGCanvas("canvas_Q2","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_xB = new TGCanvas("canvas_XB","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_W = new TGCanvas("canvas_W","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_y = new TGCanvas("canvas_y","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_z = new TGCanvas("canvas_z","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_Pt = new TGCanvas("canvas_Pt","canvas_BSA",1600,1400); 
	  
	       
	    int XCanvas=((int) Math.floor(Math.sqrt(Q2_Bins)));
	    int YCanvas=(int) Math.floor(Math.sqrt(Q2_Bins));
	    Canvas_BSA.divide(XCanvas+1, YCanvas+1);     Canvas_Dis_Q2.divide(XCanvas+1, XCanvas+1);      Canvas_Dis_xB.divide(XCanvas+1, XCanvas+1);
	    Canvas_Dis_W.divide(XCanvas+1, YCanvas+1);     Canvas_Dis_y.divide(XCanvas+1, XCanvas+1);      Canvas_Dis_z.divide(XCanvas+1, XCanvas+1);
	    Canvas_Dis_Pt.divide(XCanvas+1, XCanvas+1);
	    
	    for (int i=0; i<Q2_Bins; i++ ) {
		   System.out.println( "valore i " + i + " valore Q2 bins "+ Q2_Bins);
		   GraphErrors BSAQ2 = BSA_Plots[i][0][0][0];
		   BSAQ2.setTitleX(" Azimutal Angle φ_{trento} [rad]");
		   BSAQ2.setTitleY(" Beam Single Spin Asymmetry ");
		   BSAQ2.setTitle(" Bin "+i );
		   BSAQ2.getFunction().setOptStat("1111");
		   Canvas_BSA.cd(i);Canvas_BSA.draw(BSAQ2);
		   
		   Dis_Q2[i][0][0][0].setTitle("Q2 Bin " + i); Dis_xB[i][0][0][0].setTitle("Q2 Bin " + i);
		   Dis_W[i][0][0][0].setTitle("Q2 Bin " + i); Dis_z[i][0][0][0].setTitle("Q2 Bin " + i);
		   Dis_Pt[i][0][0][0].setTitle("Q2 Bin " + i); 
		   Canvas_Dis_Q2.cd(i);Canvas_Dis_Q2.draw(Dis_Q2[i][0][0][0]);
		   Canvas_Dis_xB.cd(i);Canvas_Dis_xB.draw(Dis_xB[i][0][0][0]);
		   Canvas_Dis_W.cd(i);Canvas_Dis_W.draw(Dis_W[i][0][0][0]);
		   Canvas_Dis_y.cd(i);Canvas_Dis_y.draw(Dis_y[i][0][0][0]);
		   Canvas_Dis_z.cd(i);Canvas_Dis_z.draw(Dis_z[i][0][0][0]);
		   Canvas_Dis_Pt.cd(i);Canvas_Dis_Pt.draw(Dis_Pt[i][0][0][0]);

		   
	   }
	    Canvas_Dis_Q2.save(Outputdir+"Q2_dist.png");
	    Canvas_Dis_xB.save(Outputdir+"xB_dist.png");
	    Canvas_Dis_W.save(Outputdir+"W_dist.png");
	    Canvas_Dis_y.save(Outputdir+"y_dist.png");
	    Canvas_Dis_z.save(Outputdir+"z_dist.png");
	    Canvas_Dis_Pt.save(Outputdir+"Pt_dist.png");
	    
	   //Compute ALU and F/F over the Q2 range
	   GraphErrors ALUvsQ2 =  new GraphErrors();
	   GraphErrors FFvsQ2 =  new GraphErrors();
       ALUvsQ2.setMarkerColor(4); FFvsQ2.setMarkerColor(5);
       ALUvsQ2.setTitle("BSA Q^2 Dependence");  FFvsQ2.setTitle("F/F Q^2 Dependence");
       ALUvsQ2.setTitleX(" Q^2 [GeV^2]");FFvsQ2.setTitleX(" Q^2 [GeV^2]");

	   for (int i=0; i<Q2_Bins; i++ ) {
		   ALUvsQ2.addPoint(Q2_Values[i][0][0][0], BSA_Values[i][0][0][0], 0, BSA_Errors[i][0][0][0]);
		   FFvsQ2.addPoint(Q2_Values[i][0][0][0], Str_F_Values[i][0][0][0], 0, Str_F_errors[i][0][0][0]);
		
		   
	   }
	   Canvas_BSA.cd(Q2_Bins);Canvas_BSA.draw(ALUvsQ2);
	   Canvas_BSA.cd(Q2_Bins+1);Canvas_BSA.draw(FFvsQ2);
	   Canvas_BSA.save(Outputdir+"BSA.png");
	   File folder = new File(Outputdir+"");
		folder.mkdir();
		File File_txt = new File(Outputdir+"/1D_Q2_table.txt");
		try {
			File_txt.createNewFile();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if(File_txt.exists()) {
			System.out.println(" Output File Found" );
			PrintWriter outP = new PrintWriter(File_txt);
			NumberFormat nf = NumberFormat.getInstance();
	        nf.setMinimumFractionDigits(4);
			outP.println("======================================================================================================================================");
		outP.println("====================================     1D Q^2 BSA Analysis from Giovanni Angelini code       =======================================");
	    outP.println(" ====== All the Chisquared/NDF are to be considered having an error of +/- 0.408 since 12 data points are used for the fits =======");
		outP.println("===================================================================================================================================");
		outP.println("|<Q2> | <xB>| <W> | <y> | <z> | <Pt> | <epsilon> | Red.ChiSqr | ALU | Err_ALU | FLU/FUU | Err_FLU/UU | Pions Nr |");
		outP.println("===================================================================================================================================");
		   for (int i=0; i<Q2_Bins; i++ ) {
		outP.println(nf.format(Q2_Values[i][0][0][0])+" "+nf.format(xB_Values[i][0][0][0])+" "+nf.format(W_Values[i][0][0][0])+" "+nf.format(y_Values[i][0][0][0])+" "+nf.format(z_Values[i][0][0][0])+" "+nf.format(Pt_Values[i][0][0][0])+" "+nf.format(eps_Values[i][0][0][0])+" "+
				nf.format(RChiSquared[i][0][0][0])+" "+nf.format(BSA_Values[i][0][0][0]) + " "+nf.format(BSA_Errors[i][0][0][0])+" "+nf.format(Str_F_Values[i][0][0][0])+" "+nf.format(Str_F_errors[i][0][0][0])+" "+nf.format(TPion_helicity[i][0][0][0]) );
		   }
		   outP.close();
		}
		  
	}
	
	private static void Print_outoput_ID1(String Outputdir) throws FileNotFoundException {
	       TGCanvas Canvas_Dis_Q2 = new TGCanvas("canvas_Q2","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_xB = new TGCanvas("canvas_XB","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_W = new TGCanvas("canvas_W","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_y = new TGCanvas("canvas_y","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_z = new TGCanvas("canvas_z","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_Pt = new TGCanvas("canvas_Pt","canvas_BSA",1600,1400); 
		TGCanvas Canvas_BSA = new TGCanvas("canvas_BSA","canvas_BSA",1600,1400); // Create Canvas with (1,1) divisions
		
		//test 
		TGCanvas Helicity = new TGCanvas("Helicity","Helicity",1600,1400); // Create Canvas with (1,1) divisions
		Helicity.divide(3, 3);

		 for (int j=0; j<XB_Bins; j++ ) {
			 if(j<9) {
					H1F HellP = new H1F("HellP",12,0,360);
					H1F HellM = new H1F("HellM",12,0,360); HellM.setLineColor(2);
				 System.out.println(" J value "+ j);
				 Helicity.cd(j);
				 for (int k=0; k<Phi_Bins; k++) {
					 
		   HellM.setBinContent(k, HelicityM[0][j][0][0][k]);  
		   HellP.setBinContent(k, HelicityP[0][j][0][0][k]);  
		  
				 }
				 HellM.setTitle(" Bin "+ j);HellP.setTitle(" Bin "+ j);
				 Helicity.draw(HellM);Helicity.draw(HellP,"same");
			 }
		 }
		
	    int XCanvas=((int) Math.floor(Math.sqrt(XB_Bins)));
	    int YCanvas=(int) Math.floor(Math.sqrt(XB_Bins));
	    Canvas_BSA.divide(XCanvas+1, YCanvas+1);     Canvas_Dis_Q2.divide(XCanvas+1, XCanvas+1);      Canvas_Dis_xB.divide(XCanvas+1, XCanvas+1);
	    Canvas_Dis_W.divide(XCanvas+1, YCanvas+1);     Canvas_Dis_y.divide(XCanvas+1, XCanvas+1);      Canvas_Dis_z.divide(XCanvas+1, XCanvas+1);
	    Canvas_Dis_Pt.divide(XCanvas+1, XCanvas+1);
	    for (int j=0; j<XB_Bins; j++ ) {
		   GraphErrors BSAxB = BSA_Plots[0][j][0][0];
		   BSAxB.setTitleX(" Azimutal Angle φ_{trento} [rad]");
		   BSAxB.setTitleY(" Beam Single Spin Asymmetry ");
		   BSAxB.setTitle(" Bin "+j );
		   BSAxB.getFunction().setOptStat("1111");
		   Canvas_BSA.cd(j);Canvas_BSA.draw(BSAxB);
		   
		   Dis_Q2[0][j][0][0].setTitle("xB Bin " + j); Dis_xB[0][j][0][0].setTitle("xB Bin " + j);
		   Dis_W[0][j][0][0].setTitle("xB Bin " + j); Dis_z[0][j][0][0].setTitle("xB Bin " + j);
		   Dis_Pt[0][j][0][0].setTitle("xB Bin " + j); 
		   Canvas_Dis_Q2.cd(j);Canvas_Dis_Q2.draw(Dis_Q2[0][j][0][0]);
		   Canvas_Dis_xB.cd(j);Canvas_Dis_xB.draw(Dis_xB[0][j][0][0]);
		   Canvas_Dis_W.cd(j);Canvas_Dis_W.draw(Dis_W[0][j][0][0]);
		   Canvas_Dis_y.cd(j);Canvas_Dis_y.draw(Dis_y[0][j][0][0]);
		   Canvas_Dis_z.cd(j);Canvas_Dis_z.draw(Dis_z[0][j][0][0]);
		   Canvas_Dis_Pt.cd(j);Canvas_Dis_Pt.draw(Dis_Pt[0][j][0][0]);
	   }
	    
	    Canvas_Dis_Q2.save(Outputdir+"Q2_dist.png");
	    Canvas_Dis_xB.save(Outputdir+"xB_dist.png");
	    Canvas_Dis_W.save(Outputdir+"W_dist.png");
	    Canvas_Dis_y.save(Outputdir+"y_dist.png");
	    Canvas_Dis_z.save(Outputdir+"z_dist.png");
	    Canvas_Dis_Pt.save(Outputdir+"Pt_dist.png");
	    
	   //Compute ALU and F/F over the Q2 range
	   GraphErrors ALUvsxB =  new GraphErrors();
	   GraphErrors FFvsxB =  new GraphErrors();
       ALUvsxB.setMarkerColor(4); FFvsxB.setMarkerColor(5);
       ALUvsxB.setTitle("BSA xB Dependence");  FFvsxB.setTitle("F/F xB Dependence");
       ALUvsxB.setTitleX(" xB ");FFvsxB.setTitleX("xB ");
	   for (int j=0; j<XB_Bins; j++ ) {
		   ALUvsxB.addPoint(xB_Values[0][j][0][0], BSA_Values[0][j][0][0], 0, BSA_Errors[0][j][0][0]);
		   FFvsxB.addPoint(xB_Values[0][j][0][0], Str_F_Values[0][j][0][0], 0, Str_F_errors[0][j][0][0]);
	   }
	   Canvas_BSA.cd(XB_Bins);Canvas_BSA.draw(ALUvsxB);
	   Canvas_BSA.cd(XB_Bins+1);Canvas_BSA.draw(FFvsxB);
	   Canvas_BSA.save(Outputdir+"BSA.png");
	   
	   File folder = new File(Outputdir+"");
		folder.mkdir();
		File File_txt = new File(Outputdir+"/1D_xB_table.txt");
		try {
			File_txt.createNewFile();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if(File_txt.exists()) {
			System.out.println(" Output File Found" );
			PrintWriter outP = new PrintWriter(File_txt);
			NumberFormat nf = NumberFormat.getInstance();
	        nf.setMinimumFractionDigits(4);
			outP.println("======================================================================================================================================");
		outP.println("====================================     1D xB BSA Analysis from Giovanni Angelini code       =======================================");
	    outP.println(" ====== All the Chisquared/NDF are to be considered having an error of +/- 0.408 since 12 data points are used for the fits =======");
		outP.println("===================================================================================================================================");
		outP.println("|<Q2> | <xB>| <W> | <y> | <z> | <Pt> | <epsilon> | Red.ChiSqr | ALU | Err_ALU | FLU/FUU | Err_FLU/UU | Pions Nr |");
		outP.println("===================================================================================================================================");
		   for (int j=0; j<XB_Bins; j++ ) {
		outP.println(nf.format(Q2_Values[0][j][0][0])+" "+nf.format(xB_Values[0][j][0][0])+" "+nf.format(W_Values[0][j][0][0])+" "+nf.format(y_Values[0][j][0][0])+" "+nf.format(z_Values[0][j][0][0])+" "+nf.format(Pt_Values[0][j][0][0])+" "+nf.format(eps_Values[0][j][0][0])+" "+
				nf.format(RChiSquared[0][j][0][0])+" "+nf.format(BSA_Values[0][j][0][0]) + " "+nf.format(BSA_Errors[0][j][0][0])+" "+nf.format(Str_F_Values[0][j][0][0])+" "+nf.format(Str_F_errors[0][j][0][0])+" "+nf.format(TPion_helicity[0][j][0][0]) );
		   }
		   outP.close();
		}
		  
	}
	
	
	private static void Print_outoput_ID2(String Outputdir) throws FileNotFoundException {
	       TGCanvas Canvas_Dis_Q2 = new TGCanvas("canvas_Q2","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_xB = new TGCanvas("canvas_XB","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_W = new TGCanvas("canvas_W","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_y = new TGCanvas("canvas_y","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_z = new TGCanvas("canvas_z","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_Pt = new TGCanvas("canvas_Pt","canvas_BSA",1600,1400); 
		TGCanvas Canvas_BSA = new TGCanvas("canvas_BSA","canvas_BSA",1600,1400); // Create Canvas with (1,1) divisions
	    int XCanvas=((int) Math.floor(Math.sqrt(Z_Bins)));
	    int YCanvas=(int) Math.floor(Math.sqrt(Z_Bins));
	    Canvas_BSA.divide(XCanvas+1, YCanvas+1);     Canvas_Dis_Q2.divide(XCanvas+1, XCanvas+1);      Canvas_Dis_xB.divide(XCanvas+1, XCanvas+1);
	    Canvas_Dis_W.divide(XCanvas+1, YCanvas+1);     Canvas_Dis_y.divide(XCanvas+1, XCanvas+1);      Canvas_Dis_z.divide(XCanvas+1, XCanvas+1);
	    Canvas_Dis_Pt.divide(XCanvas+1, XCanvas+1);
	    for (int k=0; k<Z_Bins; k++ ) {
		   GraphErrors BSAz = BSA_Plots[0][0][k][0];
		   BSAz.setTitleX(" Azimutal Angle φ_{trento} [rad]");
		   BSAz.setTitleY(" Beam Single Spin Asymmetry ");
		   BSAz.setTitle(" Bin "+k );
		   BSAz.getFunction().setOptStat("1111");
		   Canvas_BSA.cd(k);Canvas_BSA.draw(BSAz);
		   
		   Dis_Q2[0][0][k][0].setTitle("xB Bin " + k); Dis_xB[0][0][k][0].setTitle("xB Bin " + k);
		   Dis_W[0][0][k][0].setTitle("xB Bin " + k); Dis_z[0][0][k][0].setTitle("xB Bin " + k);
		   Dis_Pt[0][0][k][0].setTitle("xB Bin " + k); 
		   Canvas_Dis_Q2.cd(k);Canvas_Dis_Q2.draw(Dis_Q2[0][0][k][0]);
		   Canvas_Dis_xB.cd(k);Canvas_Dis_xB.draw(Dis_xB[0][0][k][0]);
		   Canvas_Dis_W.cd(k);Canvas_Dis_W.draw(Dis_W[0][0][k][0]);
		   Canvas_Dis_y.cd(k);Canvas_Dis_y.draw(Dis_y[0][0][k][0]);
		   Canvas_Dis_z.cd(k);Canvas_Dis_z.draw(Dis_z[0][0][k][0]);
		   Canvas_Dis_Pt.cd(k);Canvas_Dis_Pt.draw(Dis_Pt[0][0][k][0]);
	   }
	    
	    Canvas_Dis_Q2.save(Outputdir+"Q2_dist.png");
	    Canvas_Dis_xB.save(Outputdir+"xB_dist.png");
	    Canvas_Dis_W.save(Outputdir+"W_dist.png");
	    Canvas_Dis_y.save(Outputdir+"y_dist.png");
	    Canvas_Dis_z.save(Outputdir+"z_dist.png");
	    Canvas_Dis_Pt.save(Outputdir+"Pt_dist.png");
	    
	   //Compute ALU and F/F over the Q2 range
	   GraphErrors ALUvsxB =  new GraphErrors();
	   GraphErrors FFvsxB =  new GraphErrors();
    ALUvsxB.setMarkerColor(4); FFvsxB.setMarkerColor(5);
    ALUvsxB.setTitle("BSA Z Dependence");  FFvsxB.setTitle("F/F Z Dependence");
    ALUvsxB.setTitleX(" z ");FFvsxB.setTitleX(" z");
	   for (int k=0; k<Z_Bins; k++ ) {
		   ALUvsxB.addPoint(z_Values[0][0][k][0], BSA_Values[0][0][k][0], 0, BSA_Errors[0][0][k][0]);
		   FFvsxB.addPoint(z_Values[0][0][k][0], Str_F_Values[0][0][k][0], 0, Str_F_errors[0][0][k][0]);
	   }
	   Canvas_BSA.cd(Z_Bins);Canvas_BSA.draw(ALUvsxB);
	   Canvas_BSA.cd(Z_Bins+1);Canvas_BSA.draw(FFvsxB);
	   Canvas_BSA.save(Outputdir+"BSA.png");
	   
	   File folder = new File(Outputdir+"");
		folder.mkdir();
		File File_txt = new File(Outputdir+"/1D_z_table.txt");
		try {
			File_txt.createNewFile();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if(File_txt.exists()) {
			System.out.println(" Output File Found" );
			PrintWriter outP = new PrintWriter(File_txt);
			NumberFormat nf = NumberFormat.getInstance();
	        nf.setMinimumFractionDigits(4);
			outP.println("======================================================================================================================================");
		outP.println("====================================     1D Z BSA Analysis from Giovanni Angelini code       =======================================");
	    outP.println(" ====== All the Chisquared/NDF are to be considered having an error of +/- 0.408 since 12 data points are used for the fits =======");
		outP.println("===================================================================================================================================");
		outP.println("|<Q2> | <xB>| <W> | <y> | <z> | <Pt> | <epsilon> | Red.ChiSqr | ALU | Err_ALU | FLU/FUU | Err_FLU/UU | Pions Nr |");
		outP.println("===================================================================================================================================");
		   for (int k=0; k<Z_Bins; k++ ) {
		outP.println(nf.format(Q2_Values[0][0][k][0])+" "+nf.format(xB_Values[0][0][k][0])+" "+nf.format(W_Values[0][0][k][0])+" "+nf.format(y_Values[0][0][k][0])+" "+nf.format(z_Values[0][0][k][0])+" "+nf.format(Pt_Values[0][0][k][0])+" "+nf.format(eps_Values[0][0][k][0])+" "+
				nf.format(RChiSquared[0][0][k][0])+" "+nf.format(BSA_Values[0][0][k][0]) + " "+nf.format(BSA_Errors[0][0][k][0])+" "+nf.format(Str_F_Values[0][0][k][0])+" "+nf.format(Str_F_errors[0][0][k][0])+" "+nf.format(TPion_helicity[0][0][k][0]) );
		   }
		   outP.close();
		}
		  
	}
	
	
	private static void Print_outoput_ID3(String Outputdir) throws FileNotFoundException {
	       TGCanvas Canvas_Dis_Q2 = new TGCanvas("canvas_Q2","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_xB = new TGCanvas("canvas_XB","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_W = new TGCanvas("canvas_W","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_y = new TGCanvas("canvas_y","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_z = new TGCanvas("canvas_z","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_Pt = new TGCanvas("canvas_Pt","canvas_BSA",1600,1400); 
		TGCanvas Canvas_BSA = new TGCanvas("canvas_BSA","canvas_BSA",1600,1400); // Create Canvas with (1,1) divisions
	    int XCanvas=((int) Math.floor(Math.sqrt(PT_Bins)));
	    int YCanvas=(int) Math.floor(Math.sqrt(PT_Bins));
	    Canvas_BSA.divide(XCanvas+1, YCanvas+2);     Canvas_Dis_Q2.divide(XCanvas+1, XCanvas+2);      Canvas_Dis_xB.divide(XCanvas+1, XCanvas+2);
	    Canvas_Dis_W.divide(XCanvas+1, YCanvas+2);     Canvas_Dis_y.divide(XCanvas+1, XCanvas+2);      Canvas_Dis_z.divide(XCanvas+1, XCanvas+2);
	    Canvas_Dis_Pt.divide(XCanvas+1, XCanvas+2);
	    for (int l=0; l<PT_Bins; l++ ) {
		   GraphErrors BSAPT = BSA_Plots[0][0][0][l];
		   BSAPT.setTitleX(" Azimutal Angle φ_{trento} [rad]");
		   BSAPT.setTitleY(" Beam Single Spin Asymmetry ");
		   BSAPT.setTitle(" Bin "+l );
		   BSAPT.getFunction().setOptStat("1111");
		   Canvas_BSA.cd(l);Canvas_BSA.draw(BSAPT);
		   
		   Dis_Q2[0][0][0][l].setTitle("Pt Bin " + l); Dis_xB[0][0][0][l].setTitle("Pt Bin " + l);
		   Dis_W[0][0][0][l].setTitle("Pt Bin " + l); Dis_z[0][0][0][l].setTitle("Pt Bin " + l);
		   Dis_Pt[0][0][0][l].setTitle("Pt Bin " + l); 
		   Canvas_Dis_Q2.cd(l);Canvas_Dis_Q2.draw(Dis_Q2[0][0][0][l]);
		   Canvas_Dis_xB.cd(l);Canvas_Dis_xB.draw(Dis_xB[0][0][0][l]);
		   Canvas_Dis_W.cd(l);Canvas_Dis_W.draw(Dis_W[0][0][0][l]);
		   Canvas_Dis_y.cd(l);Canvas_Dis_y.draw(Dis_y[0][0][0][l]);
		   Canvas_Dis_z.cd(l);Canvas_Dis_z.draw(Dis_z[0][0][0][l]);
		   Canvas_Dis_Pt.cd(l);Canvas_Dis_Pt.draw(Dis_Pt[0][0][0][l]);
	   }
	    
	    Canvas_Dis_Q2.save(Outputdir+"Q2_dist.png");
	    Canvas_Dis_xB.save(Outputdir+"xB_dist.png");
	    Canvas_Dis_W.save(Outputdir+"W_dist.png");
	    Canvas_Dis_y.save(Outputdir+"y_dist.png");
	    Canvas_Dis_z.save(Outputdir+"z_dist.png");
	    Canvas_Dis_Pt.save(Outputdir+"Pt_dist.png");
	    
	   //Compute ALU and F/F over the Q2 range
	   GraphErrors ALUvsxB =  new GraphErrors();
	   GraphErrors FFvsxB =  new GraphErrors();
 ALUvsxB.setMarkerColor(4); FFvsxB.setMarkerColor(5);
 ALUvsxB.setTitle("BSA Pt Dependence");  FFvsxB.setTitle("F/F Pt Dependence");
 ALUvsxB.setTitleX(" Pt ");FFvsxB.setTitleX(" Pt ");
	   for (int l=0; l<(PT_Bins-1); l++ ) { //skipping the lowest bin
		   ALUvsxB.addPoint(Pt_Values[0][0][0][l], BSA_Values[0][0][0][l], 0, BSA_Errors[0][0][0][l]);
		   FFvsxB.addPoint(Pt_Values[0][0][0][l], Str_F_Values[0][0][0][l], 0, Str_F_errors[0][0][0][l]);
		   
		   }
	   Canvas_BSA.cd(PT_Bins);Canvas_BSA.draw(ALUvsxB);
	   Canvas_BSA.cd(PT_Bins+1);Canvas_BSA.draw(FFvsxB);
	   Canvas_BSA.save(Outputdir+"BSA.png");
	   
	   File folder = new File(Outputdir+"");
		folder.mkdir();
		File File_txt = new File(Outputdir+"/1D_PT_table.txt");
		try {
			File_txt.createNewFile();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if(File_txt.exists()) {
			System.out.println(" Output File Found" );
			PrintWriter outP = new PrintWriter(File_txt);
			NumberFormat nf = NumberFormat.getInstance();
	        nf.setMinimumFractionDigits(4);
			outP.println("======================================================================================================================================");
		outP.println("====================================     1D PT BSA Analysis from Giovanni Angelini code       =======================================");
	    outP.println(" ====== All the Chisquared/NDF are to be considered having an error of +/- 0.408 since 12 data points are used for the fits =======");
		outP.println("===================================================================================================================================");
		outP.println("|<Q2> | <xB>| <W> | <y> | <z> | <Pt> | <epsilon> | Red.ChiSqr | ALU | Err_ALU | FLU/FUU | Err_FLU/UU | Pions Nr |");
		outP.println("===================================================================================================================================");
		   for (int l=0; l<PT_Bins; l++ ) {
		outP.println(nf.format(Q2_Values[0][0][0][l])+" "+nf.format(xB_Values[0][0][0][l])+" "+nf.format(W_Values[0][0][0][l])+" "+nf.format(y_Values[0][0][0][l])+" "+nf.format(z_Values[0][0][0][l])+" "+nf.format(Pt_Values[0][0][0][l])+" "+nf.format(eps_Values[0][0][0][l])+" "+
				nf.format(RChiSquared[0][0][0][l])+" "+nf.format(BSA_Values[0][0][0][l]) + " "+nf.format(BSA_Errors[0][0][0][l])+" "+nf.format(Str_F_Values[0][0][0][l])+" "+nf.format(Str_F_errors[0][0][0][l])+" "+nf.format(TPion_helicity[0][0][0][l]) );
		   }
		   outP.close();
		}
		  
	}
	
	private static void Print_outoput_ID4(String Outputdir) throws FileNotFoundException {	
		// nine bins analysis
		   TGCanvas Canvas_Dis_Q2 = new TGCanvas("canvas_Q2","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_xB = new TGCanvas("canvas_XB","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_W = new TGCanvas("canvas_W","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_y = new TGCanvas("canvas_y","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_z = new TGCanvas("canvas_z","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_Dis_Pt = new TGCanvas("canvas_Pt","canvas_BSA",1600,1400); 
	       TGCanvas Canvas_BSA = new TGCanvas("canvas_BSA","canvas_BSA",1600,1400); // Create Canvas with (1,1) divisions
		
		Canvas_Dis_Q2.divide(3, 3);Canvas_Dis_xB.divide(3, 3);Canvas_Dis_W.divide(3, 3);Canvas_Dis_y.divide(3, 3);
		Canvas_Dis_z.divide(3, 3);Canvas_Dis_Pt.divide(3, 3);Canvas_BSA.divide(3, 3);
		 File folder = new File(Outputdir+"");
			folder.mkdir();
			File File_txt = new File(Outputdir+"/1D_PT_table.txt");
			try {
				File_txt.createNewFile();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if(File_txt.exists()) {
				System.out.println(" Output File Found" );
				PrintWriter outP = new PrintWriter(File_txt);
				NumberFormat nf = NumberFormat.getInstance();
		        nf.setMinimumFractionDigits(4);
				outP.println("======================================================================================================================================");
			outP.println("====================================     2D Q and X BSA Analysis from Giovanni Angelini code       =======================================");
		    outP.println(" ====== All the Chisquared/NDF are to be considered having an error of +/- 0.408 since 12 data points are used for the fits =======");
			outP.println("===================================================================================================================================");
			outP.println("|<Q2> | <xB>| <W> | <y> | <z> | <Pt> | <epsilon> | Red.ChiSqr | ALU | Err_ALU | FLU/FUU | Err_FLU/UU | Pions Nr |");
			outP.println("===================================================================================================================================");
		for(int b=1;b<=9;b++) {
		Canvas_Dis_Q2.cd(b-1);Canvas_Dis_xB.cd(b-1);Canvas_Dis_W.cd(b-1);Canvas_Dis_y.cd(b-1);
		Canvas_Dis_z.cd(b-1);Canvas_Dis_Pt.cd(b-1);Canvas_BSA.cd(b-1);
		
		if(b==1) {
			  Dis_Q2[0][0][0][0].setTitle("Bin " + b); Dis_xB[0][0][0][0].setTitle("Bin " + b);Dis_y[0][0][0][0].setTitle("Bin " + b);
			   Dis_W[0][0][0][0].setTitle("Bin " + b); Dis_z[0][0][0][0].setTitle("Bin " + b);
			   Dis_Pt[0][0][0][0].setTitle("Bin " + b); 
			   GraphErrors BSAPT = BSA_Plots[0][0][0][0];
			   BSAPT.setTitleX(" Azimutal Angle φ_{trento} [rad]");
			   BSAPT.setTitleY(" Beam Single Spin Asymmetry ");
			   BSAPT.setTitle(" Bin "+ b );
			   BSAPT.getFunction().setOptStat("1111");
				Canvas_Dis_Q2.draw(Dis_Q2[0][0][0][0]);Canvas_Dis_xB.draw(Dis_xB[0][0][0][0]);Canvas_Dis_W.draw(Dis_W[0][0][0][0]);Canvas_Dis_y.draw(Dis_y[0][0][0][0]);
				Canvas_Dis_z.draw(Dis_z[0][0][0][0]);Canvas_Dis_Pt.draw(Dis_Pt[0][0][0][0]);Canvas_BSA.draw(BSAPT);
		
				outP.println(nf.format(Q2_Values[0][0][0][0])+" "+nf.format(xB_Values[0][0][0][0])+" "+nf.format(W_Values[0][0][0][0])+" "+nf.format(y_Values[0][0][0][0])+" "+nf.format(z_Values[0][0][0][0])+" "+nf.format(Pt_Values[0][0][0][0])+" "+nf.format(eps_Values[0][0][0][0])+" "+
						nf.format(RChiSquared[0][0][0][0])+" "+nf.format(BSA_Values[0][0][0][0]) + " "+nf.format(BSA_Errors[0][0][0][0])+" "+nf.format(Str_F_Values[0][0][0][0])+" "+nf.format(Str_F_errors[0][0][0][0])+" "+nf.format(TPion_helicity[0][0][0][0]) );
		}
		
		if(b==2) {
			  Dis_Q2[0][1][0][0].setTitle("Bin " + b); Dis_xB[0][1][0][0].setTitle("Bin " + b);Dis_y[0][1][0][0].setTitle("Bin " + b);
			   Dis_W[0][1][0][0].setTitle("Bin " + b); Dis_z[0][1][0][0].setTitle("Bin " + b);
			   Dis_Pt[0][1][0][0].setTitle("Bin " + b); 
			   GraphErrors BSAPT = BSA_Plots[0][1][0][0];
			   BSAPT.setTitleX(" Azimutal Angle φ_{trento} [rad]");
			   BSAPT.setTitleY(" Beam Single Spin Asymmetry ");
			   BSAPT.setTitle(" Bin "+ b );
			   BSAPT.getFunction().setOptStat("1111");
				Canvas_Dis_Q2.draw(Dis_Q2[0][1][0][0]);Canvas_Dis_xB.draw(Dis_xB[0][1][0][0]);Canvas_Dis_W.draw(Dis_W[0][1][0][0]);Canvas_Dis_y.draw(Dis_y[0][1][0][0]);
				Canvas_Dis_z.draw(Dis_z[0][1][0][0]);Canvas_Dis_Pt.draw(Dis_Pt[0][1][0][0]);Canvas_BSA.draw(BSAPT);
				
				outP.println(nf.format(Q2_Values[0][1][0][0])+" "+nf.format(xB_Values[0][1][0][0])+" "+nf.format(W_Values[0][1][0][0])+" "+nf.format(y_Values[0][1][0][0])+" "+nf.format(z_Values[0][1][0][0])+" "+nf.format(Pt_Values[0][1][0][0])+" "+nf.format(eps_Values[0][0][0][0])+" "+
						nf.format(RChiSquared[0][1][0][0])+" "+nf.format(BSA_Values[0][1][0][0]) + " "+nf.format(BSA_Errors[0][1][0][0])+" "+nf.format(Str_F_Values[0][1][0][0])+" "+nf.format(Str_F_errors[0][1][0][0])+" "+nf.format(TPion_helicity[0][1][0][0]) );
		}
		if(b==3) {
			  Dis_Q2[1][1][0][0].setTitle("Bin " + b); Dis_xB[1][1][0][0].setTitle("Bin " + b);Dis_y[1][1][0][0].setTitle("Bin " + b);
			   Dis_W[1][1][0][0].setTitle("Bin " + b); Dis_z[1][1][0][0].setTitle("Bin " + b);
			   Dis_Pt[1][1][0][0].setTitle("Bin " + b); 
			   GraphErrors BSAPT = BSA_Plots[1][1][0][0];
			   BSAPT.setTitleX(" Azimutal Angle φ_{trento} [rad]");
			   BSAPT.setTitleY(" Beam Single Spin Asymmetry ");
			   BSAPT.setTitle(" Bin "+ b );
			   BSAPT.getFunction().setOptStat("1111");
				Canvas_Dis_Q2.draw(Dis_Q2[1][1][0][0]);Canvas_Dis_xB.draw(Dis_xB[1][1][0][0]);Canvas_Dis_W.draw(Dis_W[1][1][0][0]);Canvas_Dis_y.draw(Dis_y[1][1][0][0]);
				Canvas_Dis_z.draw(Dis_z[1][1][0][0]);Canvas_Dis_Pt.draw(Dis_Pt[1][1][0][0]);Canvas_BSA.draw(BSAPT);
	
				outP.println(nf.format(Q2_Values[1][1][0][0])+" "+nf.format(xB_Values[1][1][0][0])+" "+nf.format(W_Values[1][1][0][0])+" "+nf.format(y_Values[1][1][0][0])+" "+nf.format(z_Values[1][1][0][0])+" "+nf.format(Pt_Values[1][1][0][0])+" "+nf.format(eps_Values[0][0][0][0])+" "+
						nf.format(RChiSquared[1][1][0][0])+" "+nf.format(BSA_Values[1][1][0][0]) + " "+nf.format(BSA_Errors[1][1][0][0])+" "+nf.format(Str_F_Values[1][1][0][0])+" "+nf.format(Str_F_errors[1][1][0][0])+" "+nf.format(TPion_helicity[1][1][0][0]) );
		
		}
		if(b==4) {
			  Dis_Q2[0][2][0][0].setTitle("Bin " + b); Dis_xB[0][2][0][0].setTitle("Bin " + b);Dis_y[0][2][0][0].setTitle("Bin " + b);
			   Dis_W[0][2][0][0].setTitle("Bin " + b); Dis_z[0][2][0][0].setTitle("Bin " + b);
			   Dis_Pt[0][2][0][0].setTitle("Bin " + b); 
			   GraphErrors BSAPT = BSA_Plots[0][2][0][0];
			   BSAPT.setTitleX(" Azimutal Angle φ_{trento} [rad]");
			   BSAPT.setTitleY(" Beam Single Spin Asymmetry ");
			   BSAPT.setTitle(" Bin "+ b );
			   BSAPT.getFunction().setOptStat("1111");
				Canvas_Dis_Q2.draw(Dis_Q2[0][2][0][0]);Canvas_Dis_xB.draw(Dis_xB[0][2][0][0]);Canvas_Dis_W.draw(Dis_W[0][2][0][0]);Canvas_Dis_y.draw(Dis_y[0][2][0][0]);
				Canvas_Dis_z.draw(Dis_z[0][2][0][0]);Canvas_Dis_Pt.draw(Dis_Pt[0][2][0][0]);Canvas_BSA.draw(BSAPT);
	
				outP.println(nf.format(Q2_Values[0][2][0][0])+" "+nf.format(xB_Values[0][2][0][0])+" "+nf.format(W_Values[0][2][0][0])+" "+nf.format(y_Values[0][2][0][0])+" "+nf.format(z_Values[0][2][0][0])+" "+nf.format(Pt_Values[0][2][0][0])+" "+nf.format(eps_Values[0][0][0][0])+" "+
						nf.format(RChiSquared[0][2][0][0])+" "+nf.format(BSA_Values[0][2][0][0]) + " "+nf.format(BSA_Errors[0][2][0][0])+" "+nf.format(Str_F_Values[0][2][0][0])+" "+nf.format(Str_F_errors[0][2][0][0])+" "+nf.format(TPion_helicity[0][2][0][0]) );
		
		}
		
		if(b==5) {
			  Dis_Q2[1][2][0][0].setTitle("Bin " + b); Dis_xB[1][2][0][0].setTitle("Bin " + b);Dis_y[1][2][0][0].setTitle("Bin " + b);
			   Dis_W[1][2][0][0].setTitle("Bin " + b); Dis_z[1][2][0][0].setTitle("Bin " + b);
			   Dis_Pt[1][2][0][0].setTitle("Bin " + b); 
			   GraphErrors BSAPT = BSA_Plots[1][2][0][0];
			   BSAPT.setTitleX(" Azimutal Angle φ_{trento} [rad]");
			   BSAPT.setTitleY(" Beam Single Spin Asymmetry ");
			   BSAPT.setTitle(" Bin "+ b );
			   BSAPT.getFunction().setOptStat("1111");
				Canvas_Dis_Q2.draw(Dis_Q2[1][2][0][0]);Canvas_Dis_xB.draw(Dis_xB[1][2][0][0]);Canvas_Dis_W.draw(Dis_W[1][2][0][0]);Canvas_Dis_y.draw(Dis_y[1][2][0][0]);
				Canvas_Dis_z.draw(Dis_z[1][2][0][0]);Canvas_Dis_Pt.draw(Dis_Pt[1][2][0][0]);Canvas_BSA.draw(BSAPT);
				
				outP.println(nf.format(Q2_Values[1][2][0][0])+" "+nf.format(xB_Values[1][2][0][0])+" "+nf.format(W_Values[1][2][0][0])+" "+nf.format(y_Values[1][2][0][0])+" "+nf.format(z_Values[1][2][0][0])+" "+nf.format(Pt_Values[1][2][0][0])+" "+nf.format(eps_Values[0][0][0][0])+" "+
						nf.format(RChiSquared[1][2][0][0])+" "+nf.format(BSA_Values[1][2][0][0]) + " "+nf.format(BSA_Errors[1][2][0][0])+" "+nf.format(Str_F_Values[1][2][0][0])+" "+nf.format(Str_F_errors[1][2][0][0])+" "+nf.format(TPion_helicity[1][2][0][0]) );
		
		}
		if(b==6) {
			  Dis_Q2[0][3][0][0].setTitle("Bin " + b); Dis_xB[0][3][0][0].setTitle("Bin " + b);Dis_y[0][3][0][0].setTitle("Bin " + b);
			   Dis_W[0][3][0][0].setTitle("Bin " + b); Dis_z[0][3][0][0].setTitle("Bin " + b);
			   Dis_Pt[0][3][0][0].setTitle("Bin " + b); 
			   GraphErrors BSAPT = BSA_Plots[0][3][0][0];
			   BSAPT.setTitleX(" Azimutal Angle φ_{trento} [rad]");
			   BSAPT.setTitleY(" Beam Single Spin Asymmetry ");
			   BSAPT.setTitle(" Bin "+ b );
			   BSAPT.getFunction().setOptStat("1111");
				Canvas_Dis_Q2.draw(Dis_Q2[0][3][0][0]);Canvas_Dis_xB.draw(Dis_xB[0][3][0][0]);Canvas_Dis_W.draw(Dis_W[0][3][0][0]);Canvas_Dis_y.draw(Dis_y[0][3][0][0]);
				Canvas_Dis_z.draw(Dis_z[0][3][0][0]);Canvas_Dis_Pt.draw(Dis_Pt[0][3][0][0]);Canvas_BSA.draw(BSAPT);
		
				outP.println(nf.format(Q2_Values[0][3][0][0])+" "+nf.format(xB_Values[0][3][0][0])+" "+nf.format(W_Values[0][3][0][0])+" "+nf.format(y_Values[0][3][0][0])+" "+nf.format(z_Values[0][3][0][0])+" "+nf.format(Pt_Values[0][3][0][0])+" "+nf.format(eps_Values[0][0][0][0])+" "+
						nf.format(RChiSquared[0][3][0][0])+" "+nf.format(BSA_Values[0][3][0][0]) + " "+nf.format(BSA_Errors[0][3][0][0])+" "+nf.format(Str_F_Values[0][3][0][0])+" "+nf.format(Str_F_errors[0][3][0][0])+" "+nf.format(TPion_helicity[0][3][0][0]) );
				
		}
		if(b==7) {
			  Dis_Q2[1][3][0][0].setTitle("Bin " + b); Dis_xB[1][3][0][0].setTitle("Bin " + b);Dis_y[1][3][0][0].setTitle("Bin " + b);
			   Dis_W[1][3][0][0].setTitle("Bin " + b); Dis_z[1][3][0][0].setTitle("Bin " + b);
			   Dis_Pt[1][3][0][0].setTitle("Bin " + b); 
			   GraphErrors BSAPT = BSA_Plots[1][3][0][0];
			   BSAPT.setTitleX(" Azimutal Angle φ_{trento} [rad]");
			   BSAPT.setTitleY(" Beam Single Spin Asymmetry ");
			   BSAPT.setTitle(" Bin "+ b );
			   BSAPT.getFunction().setOptStat("1111");
				Canvas_Dis_Q2.draw(Dis_Q2[1][3][0][0]);Canvas_Dis_xB.draw(Dis_xB[1][3][0][0]);Canvas_Dis_W.draw(Dis_W[1][3][0][0]);Canvas_Dis_y.draw(Dis_y[1][3][0][0]);
				Canvas_Dis_z.draw(Dis_z[1][3][0][0]);Canvas_Dis_Pt.draw(Dis_Pt[1][3][0][0]);Canvas_BSA.draw(BSAPT);
		
				outP.println(nf.format(Q2_Values[1][3][0][0])+" "+nf.format(xB_Values[1][3][0][0])+" "+nf.format(W_Values[1][3][0][0])+" "+nf.format(y_Values[1][3][0][0])+" "+nf.format(z_Values[1][3][0][0])+" "+nf.format(Pt_Values[1][3][0][0])+" "+nf.format(eps_Values[0][0][0][0])+" "+
						nf.format(RChiSquared[1][3][0][0])+" "+nf.format(BSA_Values[1][3][0][0]) + " "+nf.format(BSA_Errors[1][3][0][0])+" "+nf.format(Str_F_Values[1][3][0][0])+" "+nf.format(Str_F_errors[1][3][0][0])+" "+nf.format(TPion_helicity[1][3][0][0]) );
		}
		
		if(b==8) {
			  Dis_Q2[0][4][0][0].setTitle("Bin " + b); Dis_xB[0][4][0][0].setTitle("Bin " + b);Dis_y[0][4][0][0].setTitle("Bin " + b);
			   Dis_W[0][4][0][0].setTitle("Bin " + b); Dis_z[0][4][0][0].setTitle("Bin " + b);
			   Dis_Pt[0][4][0][0].setTitle("Bin " + b); 
			   GraphErrors BSAPT = BSA_Plots[0][4][0][0];
			   BSAPT.setTitleX(" Azimutal Angle φ_{trento} [rad]");
			   BSAPT.setTitleY(" Beam Single Spin Asymmetry ");
			   BSAPT.setTitle(" Bin "+ b );
			   BSAPT.getFunction().setOptStat("1111");
				Canvas_Dis_Q2.draw(Dis_Q2[0][4][0][0]);Canvas_Dis_xB.draw(Dis_xB[0][4][0][0]);Canvas_Dis_W.draw(Dis_W[0][4][0][0]);Canvas_Dis_y.draw(Dis_y[0][4][0][0]);
				Canvas_Dis_z.draw(Dis_z[0][4][0][0]);Canvas_Dis_Pt.draw(Dis_Pt[0][4][0][0]);Canvas_BSA.draw(BSAPT);
		
				outP.println(nf.format(Q2_Values[0][4][0][0])+" "+nf.format(xB_Values[0][4][0][0])+" "+nf.format(W_Values[0][4][0][0])+" "+nf.format(y_Values[0][4][0][0])+" "+nf.format(z_Values[0][4][0][0])+" "+nf.format(Pt_Values[0][4][0][0])+" "+nf.format(eps_Values[0][0][0][0])+" "+
						nf.format(RChiSquared[0][4][0][0])+" "+nf.format(BSA_Values[0][4][0][0]) + " "+nf.format(BSA_Errors[0][4][0][0])+" "+nf.format(Str_F_Values[0][4][0][0])+" "+nf.format(Str_F_errors[0][4][0][0])+" "+nf.format(TPion_helicity[0][4][0][0]) );
		
		}
		if(b==9) {
			  Dis_Q2[1][4][0][0].setTitle("Bin " + b); Dis_xB[1][4][0][0].setTitle("Bin " + b);Dis_y[1][4][0][0].setTitle("Bin " + b);
			   Dis_W[1][4][0][0].setTitle("Bin " + b); Dis_z[1][4][0][0].setTitle("Bin " + b);
			   Dis_Pt[1][4][0][0].setTitle("Bin " + b); 
			   GraphErrors BSAPT = BSA_Plots[1][4][0][0];
			   BSAPT.setTitleX(" Azimutal Angle φ_{trento} [rad]");
			   BSAPT.setTitleY(" Beam Single Spin Asymmetry ");
			   BSAPT.setTitle(" Bin "+ b );
			   BSAPT.getFunction().setOptStat("1111");
			   
				Canvas_Dis_Q2.draw(Dis_Q2[1][4][0][0]);Canvas_Dis_xB.draw(Dis_xB[1][4][0][0]);Canvas_Dis_W.draw(Dis_W[1][4][0][0]);Canvas_Dis_y.draw(Dis_y[1][4][0][0]);
				Canvas_Dis_z.draw(Dis_z[1][4][0][0]);Canvas_Dis_Pt.draw(Dis_Pt[1][4][0][0]);Canvas_BSA.draw(BSAPT);
		}
	}
		outP.close();
			}
    Canvas_Dis_Q2.save(Outputdir+"Q2_dist.png");
    Canvas_Dis_xB.save(Outputdir+"xB_dist.png");
    Canvas_Dis_W.save(Outputdir+"W_dist.png");
    Canvas_Dis_y.save(Outputdir+"y_dist.png");
    Canvas_Dis_z.save(Outputdir+"z_dist.png");
    Canvas_Dis_Pt.save(Outputdir+"Pt_dist.png");
    Canvas_BSA.save(Outputdir+"BSA.png");
   
	}
	
	private static void Print_outoput_ID5(String Outputdir) throws IOException {	
		
		for(int b=1;b<=9;b++) {
			
			if(b==1) {
			String Newfolder= Outputdir+"/ElBin_"+b+"/";
			 File folder = new File(Newfolder);
				folder.mkdir();
				File File_txt = new File(Newfolder+"/Bin_"+b+"_4D_z_PT_table.txt");
				
			
				if (!File_txt.getParentFile().exists())
					File_txt.getParentFile().mkdirs();
				if (!File_txt.exists())File_txt.createNewFile();
				
				
				
				
				try {
					File_txt.createNewFile();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(File_txt.exists()) {
					System.out.println(" Output File Found" );
					PrintWriter outP = new PrintWriter(File_txt);
					NumberFormat nf = NumberFormat.getInstance();
			        nf.setMinimumFractionDigits(4);
					outP.println("======================================================================================================================================");
				outP.println("====================================     2D Q and X BSA Analysis from Giovanni Angelini code       =======================================");
			    outP.println(" ====== All the Chisquared/NDF are to be considered having an error of +/- 0.408 since 12 data points are used for the fits =======");
				outP.println("===================================================================================================================================");
				outP.println("|<Q2> | <xB>| <W> | <y> | <z> | <Pt> | <epsilon> | Red.ChiSqr | ALU | Err_ALU | FLU/FUU | Err_FLU/UU | Pions Nr |");
				outP.println("===================================================================================================================================");
				
				for(int l=0; l<PT_Bins;l++) {
					
				    for(int k=0;k<Z_Bins;k++) {
				    	
					outP.println(nf.format(Q2_Values[0][0][k][l])+" "+nf.format(xB_Values[0][0][k][l])+" "+nf.format(W_Values[0][0][k][l])+" "+nf.format(y_Values[0][0][k][l])+" "+nf.format(z_Values[0][0][k][l])+" "+nf.format(Pt_Values[0][0][k][l])+" "+nf.format(eps_Values[0][0][k][l])+" "+
							nf.format(RChiSquared[0][0][k][l])+" "+nf.format(BSA_Values[0][0][k][l]) + " "+nf.format(BSA_Values[0][0][k][l])+" "+nf.format(Str_F_Values[0][0][k][l])+" "+nf.format(Str_F_errors[0][0][k][l])+" "+nf.format(TPion_helicity[0][0][k][l]) );
				}
				   
			}
				
			outP.close();
			
			TGCanvas Canvas_As_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
			TGCanvas Canvas_FF_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
			int XCanvas=((int) Math.floor(Math.sqrt(PT_Bins)));
			    int YCanvas=(int) Math.floor(Math.sqrt(PT_Bins));
			    Canvas_As_PT.divide(XCanvas+1, YCanvas+1);Canvas_FF_PT.divide(XCanvas+1, YCanvas+1);
			    for(int l=0; l<PT_Bins;l++) {
			    	Canvas_As_PT.cd(l);Canvas_FF_PT.cd(l);
					GraphErrors BSAvsPt = new GraphErrors();BSAvsPt.setMarkerColor(1);BSAvsPt.setTitle("PT Bin:"+l);;BSAvsPt.setTitleX(" z ");BSAvsPt.setTitleY("ALU");
					GraphErrors FFvsPt = new GraphErrors();FFvsPt.setMarkerColor(2);FFvsPt.setTitle("PT Bin:"+l );FFvsPt.setTitleX(" z ");FFvsPt.setTitleY("FLU/FUU");
					for(int k=2;k<Z_Bins;k++) {
					
						BSAvsPt.addPoint(z_Values[0][0][k][l], BSA_Values[0][0][k][l], 0, BSA_Values[0][0][k][l]);
				    	FFvsPt.addPoint(z_Values[0][0][k][l], Str_F_Values[0][0][k][l], 0, Str_F_errors[0][0][k][l]);
					}
					 Canvas_As_PT.draw(BSAvsPt);Canvas_FF_PT.draw(FFvsPt);
				}
			    Canvas_As_PT.save(Newfolder+"/BSA_fixedPt_vs_Z.png");
				Canvas_FF_PT.save(Newfolder+"/FF_fixedPt_vs_Z.png");
				
				TGCanvas Canvas_As_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				TGCanvas Canvas_FF_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				int XCanvasZ=((int) Math.floor(Math.sqrt(Z_Bins)));
				    int YCanvasZ=(int) Math.floor(Math.sqrt(Z_Bins));
				    Canvas_As_z.divide(XCanvasZ, YCanvasZ);Canvas_FF_z.divide(XCanvasZ, YCanvasZ);
				    for(int k=0; k<Z_Bins;k++) {
				    	Canvas_As_z.cd(k);Canvas_FF_z.cd(k);
						GraphErrors BSAvsz = new GraphErrors();BSAvsz.setMarkerColor(3);BSAvsz.setTitle("z Bin:"+k);;BSAvsz.setTitleX("Pt [GeV]");BSAvsz.setTitleY("ALU");
						GraphErrors FFvsz = new GraphErrors();FFvsz.setMarkerColor(4);FFvsz.setTitle("z Bin:"+k);FFvsz.setTitleX("Pt [GeV]");FFvsz.setTitleY("FLU/FUU");
						for(int l=1;l<PT_Bins-1;l++) {
							BSAvsz.addPoint(Pt_Values[0][0][k][l], BSA_Values[0][0][k][l], 0, BSA_Values[0][0][k][l]);
					    	FFvsz.addPoint(Pt_Values[0][0][k][l], Str_F_Values[0][0][k][l], 0, Str_F_errors[0][0][k][l]);
						}
						 Canvas_As_z.draw(BSAvsz);Canvas_FF_z.draw(FFvsz);
					}
				    Canvas_As_z.save(Newfolder+"/BSA_fixedZ_vs_Pt.png");
					Canvas_FF_z.save(Newfolder+"/FF_fixedZ_vs_Pt.png");
				
				}
			}
			
			if(b==2) {
			String Newfolder= Outputdir+"/ElBin_"+b+"/";
			 File folder = new File(Newfolder);
				folder.mkdir();
				File File_txt = new File(Newfolder+"/Bin_"+b+"_4D_z_PT_table.txt");
				

				if (!File_txt.getParentFile().exists())
					File_txt.getParentFile().mkdirs();
				if (!File_txt.exists())File_txt.createNewFile();
				try {
					File_txt.createNewFile();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(File_txt.exists()) {
					System.out.println(" Output File Found" );
					PrintWriter outP = new PrintWriter(File_txt);
					NumberFormat nf = NumberFormat.getInstance();
			        nf.setMinimumFractionDigits(4);
					outP.println("======================================================================================================================================");
				outP.println("====================================     2D Q and X BSA Analysis from Giovanni Angelini code       =======================================");
			    outP.println(" ====== All the Chisquared/NDF are to be considered having an error of +/- 0.408 since 12 data points are used for the fits =======");
				outP.println("===================================================================================================================================");
				outP.println("|<Q2> | <xB>| <W> | <y> | <z> | <Pt> | <epsilon> | Red.ChiSqr | ALU | Err_ALU | FLU/FUU | Err_FLU/UU | Pions Nr |");
				outP.println("===================================================================================================================================");
			for(int l=0; l<PT_Bins;l++) {
				for(int k=0;k<Z_Bins;k++) {
					outP.println(nf.format(Q2_Values[0][1][k][l])+" "+nf.format(xB_Values[0][1][k][l])+" "+nf.format(W_Values[0][1][k][l])+" "+nf.format(y_Values[0][1][k][l])+" "+nf.format(z_Values[0][1][k][l])+" "+nf.format(Pt_Values[0][1][k][l])+" "+nf.format(eps_Values[0][1][k][l])+" "+
							nf.format(RChiSquared[0][1][k][l])+" "+nf.format(BSA_Values[0][1][k][l]) + " "+nf.format(BSA_Errors[0][1][k][l])+" "+nf.format(Str_F_Values[0][1][k][l])+" "+nf.format(Str_F_errors[0][1][k][l])+" "+nf.format(TPion_helicity[0][1][k][l]) );
				}
			}
			outP.close();
			TGCanvas Canvas_As_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
			TGCanvas Canvas_FF_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
			int XCanvas=((int) Math.floor(Math.sqrt(PT_Bins)));
			    int YCanvas=(int) Math.floor(Math.sqrt(PT_Bins));
			    Canvas_As_PT.divide(XCanvas+1, YCanvas+1);Canvas_FF_PT.divide(XCanvas+1, YCanvas+1);
			    for(int l=0; l<PT_Bins;l++) {
			    	Canvas_As_PT.cd(l);Canvas_FF_PT.cd(l);
					GraphErrors BSAvsPt = new GraphErrors();BSAvsPt.setMarkerColor(1);BSAvsPt.setTitle("PT Bin:"+l);;BSAvsPt.setTitleX(" z ");BSAvsPt.setTitleY("ALU");
					GraphErrors FFvsPt = new GraphErrors();FFvsPt.setMarkerColor(2);FFvsPt.setTitle("PT Bin:"+l );FFvsPt.setTitleX(" z ");FFvsPt.setTitleY("FLU/FUU");
					for(int k=2;k<Z_Bins;k++) {
					
						BSAvsPt.addPoint(z_Values[0][1][k][l], BSA_Values[0][1][k][l], 0, BSA_Values[0][1][k][l]);
				    	FFvsPt.addPoint(z_Values[0][1][k][l], Str_F_Values[0][1][k][l], 0, Str_F_errors[0][1][k][l]);
					}
					 Canvas_As_PT.draw(BSAvsPt);Canvas_FF_PT.draw(FFvsPt);
				}
			    Canvas_As_PT.save(Newfolder+"/BSA_fixedPt_vs_Z.png");
				Canvas_FF_PT.save(Newfolder+"/FF_fixedPt_vs_Z.png");
				
				TGCanvas Canvas_As_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				TGCanvas Canvas_FF_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				int XCanvasZ=((int) Math.floor(Math.sqrt(Z_Bins)));
				    int YCanvasZ=(int) Math.floor(Math.sqrt(Z_Bins));
				    Canvas_As_z.divide(XCanvasZ, YCanvasZ);Canvas_FF_z.divide(XCanvasZ, YCanvasZ);
				    for(int k=0; k<Z_Bins;k++) {
				    	Canvas_As_z.cd(k);Canvas_FF_z.cd(k);
						GraphErrors BSAvsz = new GraphErrors();BSAvsz.setMarkerColor(3);BSAvsz.setTitle("z Bin:"+k);;BSAvsz.setTitleX("Pt [GeV]");BSAvsz.setTitleY("ALU");
						GraphErrors FFvsz = new GraphErrors();FFvsz.setMarkerColor(4);FFvsz.setTitle("z Bin:"+k);FFvsz.setTitleX("Pt [GeV]");FFvsz.setTitleY("FLU/FUU");
						for(int l=1;l<PT_Bins-1;l++) {
							BSAvsz.addPoint(Pt_Values[0][1][k][l], BSA_Values[0][1][k][l], 0, BSA_Values[0][1][k][l]);
					    	FFvsz.addPoint(Pt_Values[0][1][k][l], Str_F_Values[0][1][k][l], 0, Str_F_errors[0][1][k][l]);
						}
						 Canvas_As_z.draw(BSAvsz);Canvas_FF_z.draw(FFvsz);
					}
				    Canvas_As_z.save(Newfolder+"/BSA_fixedZ_vs_Pt.png");
					Canvas_FF_z.save(Newfolder+"/FF_fixedZ_vs_Pt.png");
				}
				}
			
			if(b==3) {
				String Newfolder= Outputdir+"/ElBin_"+b+"/";
				 File folder = new File(Newfolder);
					folder.mkdir();
					File File_txt = new File(Newfolder+"/Bin_"+b+"_4D_z_PT_table.txt");
					try {
						File_txt.createNewFile();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					if(File_txt.exists()) {
						System.out.println(" Output File Found" );
						PrintWriter outP = new PrintWriter(File_txt);
						NumberFormat nf = NumberFormat.getInstance();
				        nf.setMinimumFractionDigits(4);
						outP.println("======================================================================================================================================");
					outP.println("====================================     2D Q and X BSA Analysis from Giovanni Angelini code       =======================================");
				    outP.println(" ====== All the Chisquared/NDF are to be considered having an error of +/- 0.408 since 12 data points are used for the fits =======");
					outP.println("===================================================================================================================================");
					outP.println("|<Q2> | <xB>| <W> | <y> | <z> | <Pt> | <epsilon> | Red.ChiSqr | ALU | Err_ALU | FLU/FUU | Err_FLU/UU | Pions Nr |");
					outP.println("===================================================================================================================================");
				for(int l=0; l<PT_Bins;l++) {
					for(int k=0;k<Z_Bins;k++) {
						
						outP.println(nf.format(Q2_Values[1][1][k][l])+" "+nf.format(xB_Values[1][1][k][l])+" "+nf.format(W_Values[1][1][k][l])+" "+nf.format(y_Values[1][1][k][l])+" "+nf.format(z_Values[1][1][k][l])+" "+nf.format(Pt_Values[1][1][k][l])+" "+nf.format(eps_Values[1][1][k][l])+" "+
								nf.format(RChiSquared[1][1][k][l])+" "+nf.format(BSA_Values[1][1][k][l]) + " "+nf.format(BSA_Errors[1][1][k][l])+" "+nf.format(Str_F_Values[1][1][k][l])+" "+nf.format(Str_F_errors[1][1][k][l])+" "+nf.format(TPion_helicity[1][1][k][l]) );
					}
				}
				outP.close();
				TGCanvas Canvas_As_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				TGCanvas Canvas_FF_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				int XCanvas=((int) Math.floor(Math.sqrt(PT_Bins)));
				    int YCanvas=(int) Math.floor(Math.sqrt(PT_Bins));
				    Canvas_As_PT.divide(XCanvas+1, YCanvas+1);Canvas_FF_PT.divide(XCanvas+1, YCanvas+1);
				    for(int l=0; l<PT_Bins;l++) {
				    	Canvas_As_PT.cd(l);Canvas_FF_PT.cd(l);
						GraphErrors BSAvsPt = new GraphErrors();BSAvsPt.setMarkerColor(1);BSAvsPt.setTitle("PT Bin:"+l);;BSAvsPt.setTitleX(" z ");BSAvsPt.setTitleY("ALU");
						GraphErrors FFvsPt = new GraphErrors();FFvsPt.setMarkerColor(2);FFvsPt.setTitle("PT Bin:"+l );FFvsPt.setTitleX(" z ");FFvsPt.setTitleY("FLU/FUU");
						for(int k=2;k<Z_Bins;k++) {
						
							BSAvsPt.addPoint(z_Values[1][1][k][l], BSA_Values[1][1][k][l], 0, BSA_Values[1][1][k][l]);
					    	FFvsPt.addPoint(z_Values[1][1][k][l], Str_F_Values[1][1][k][l], 0, Str_F_errors[1][1][k][l]);
						}
						
						 Canvas_As_PT.draw(BSAvsPt);Canvas_FF_PT.draw(FFvsPt);
					}
				    Canvas_As_PT.save(Newfolder+"/BSA_fixedPt_vs_Z.png");
					Canvas_FF_PT.save(Newfolder+"/FF_fixedPt_vs_Z.png");
					
					TGCanvas Canvas_As_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
					TGCanvas Canvas_FF_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
					int XCanvasZ=((int) Math.floor(Math.sqrt(Z_Bins)));
					    int YCanvasZ=(int) Math.floor(Math.sqrt(Z_Bins));
					    Canvas_As_z.divide(XCanvasZ, YCanvasZ);Canvas_FF_z.divide(XCanvasZ, YCanvasZ);
					    for(int k=0; k<Z_Bins;k++) {
					    	Canvas_As_z.cd(k);Canvas_FF_z.cd(k);
							GraphErrors BSAvsz = new GraphErrors();BSAvsz.setMarkerColor(3);BSAvsz.setTitle("z Bin:"+k);;BSAvsz.setTitleX("Pt [GeV]");BSAvsz.setTitleY("ALU");
							GraphErrors FFvsz = new GraphErrors();FFvsz.setMarkerColor(4);FFvsz.setTitle("z Bin:"+k);FFvsz.setTitleX("Pt [GeV]");FFvsz.setTitleY("FLU/FUU");
							for(int l=1;l<PT_Bins-1;l++) {
								
								BSAvsz.addPoint(Pt_Values[1][1][k][l], BSA_Values[1][1][k][l], 0, BSA_Values[1][1][k][l]);
						    	FFvsz.addPoint(Pt_Values[1][1][k][l], Str_F_Values[1][1][k][l], 0, Str_F_errors[1][1][k][l]);
							
							}
							 Canvas_As_z.draw(BSAvsz);Canvas_FF_z.draw(FFvsz);
						}
					    Canvas_As_z.save(Newfolder+"/BSA_fixedZ_vs_Pt.png");
						Canvas_FF_z.save(Newfolder+"/FF_fixedZ_vs_Pt.png");
					}
					}
			if(b==4) {
				String Newfolder= Outputdir+"/ElBin_"+b+"/";
				 File folder = new File(Newfolder);
					folder.mkdir();
					File File_txt = new File(Newfolder+"/Bin_"+b+"_4D_z_PT_table.txt");
					try {
						File_txt.createNewFile();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					if(File_txt.exists()) {
						System.out.println(" Output File Found" );
						PrintWriter outP = new PrintWriter(File_txt);
						NumberFormat nf = NumberFormat.getInstance();
				        nf.setMinimumFractionDigits(4);
						outP.println("======================================================================================================================================");
					outP.println("====================================     2D Q and X BSA Analysis from Giovanni Angelini code       =======================================");
				    outP.println(" ====== All the Chisquared/NDF are to be considered having an error of +/- 0.408 since 12 data points are used for the fits =======");
					outP.println("===================================================================================================================================");
					outP.println("|<Q2> | <xB>| <W> | <y> | <z> | <Pt> | <epsilon> | Red.ChiSqr | ALU | Err_ALU | FLU/FUU | Err_FLU/UU | Pions Nr |");
					outP.println("===================================================================================================================================");
				for(int l=0; l<PT_Bins;l++) {
					for(int k=0;k<Z_Bins;k++) {
						outP.println(nf.format(Q2_Values[0][2][k][l])+" "+nf.format(xB_Values[0][2][k][l])+" "+nf.format(W_Values[0][2][k][l])+" "+nf.format(y_Values[0][2][k][l])+" "+nf.format(z_Values[0][2][k][l])+" "+nf.format(Pt_Values[0][2][k][l])+" "+nf.format(eps_Values[0][2][k][l])+" "+
								nf.format(RChiSquared[0][2][k][l])+" "+nf.format(BSA_Values[0][2][k][l]) + " "+nf.format(BSA_Errors[0][2][k][l])+" "+nf.format(Str_F_Values[0][2][k][l])+" "+nf.format(Str_F_errors[0][2][k][l])+" "+nf.format(TPion_helicity[0][2][k][l]) );
					}
				}
				outP.close();
				TGCanvas Canvas_As_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				TGCanvas Canvas_FF_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				int XCanvas=((int) Math.floor(Math.sqrt(PT_Bins)));
				    int YCanvas=(int) Math.floor(Math.sqrt(PT_Bins));
				    Canvas_As_PT.divide(XCanvas+1, YCanvas+1);Canvas_FF_PT.divide(XCanvas+1, YCanvas+1);
				    for(int l=0; l<PT_Bins;l++) {
				    	Canvas_As_PT.cd(l);Canvas_FF_PT.cd(l);
						GraphErrors BSAvsPt = new GraphErrors();BSAvsPt.setMarkerColor(1);BSAvsPt.setTitle("PT Bin:"+l);;BSAvsPt.setTitleX(" z ");BSAvsPt.setTitleY("ALU");
						GraphErrors FFvsPt = new GraphErrors();FFvsPt.setMarkerColor(2);FFvsPt.setTitle("PT Bin:"+l );FFvsPt.setTitleX(" z ");FFvsPt.setTitleY("FLU/FUU");
						for(int k=2;k<Z_Bins;k++) {
						
							BSAvsPt.addPoint(z_Values[0][2][k][l], BSA_Values[0][2][k][l], 0, BSA_Values[0][2][k][l]);
					    	FFvsPt.addPoint(z_Values[0][2][k][l], Str_F_Values[0][2][k][l], 0, Str_F_errors[0][2][k][l]);
						}
						 Canvas_As_PT.draw(BSAvsPt);Canvas_FF_PT.draw(FFvsPt);
					}
				    Canvas_As_PT.save(Newfolder+"/BSA_fixedPt_vs_Z.png");
					Canvas_FF_PT.save(Newfolder+"/FF_fixedPt_vs_Z.png");
					
					TGCanvas Canvas_As_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
					TGCanvas Canvas_FF_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
					int XCanvasZ=((int) Math.floor(Math.sqrt(Z_Bins)));
					    int YCanvasZ=(int) Math.floor(Math.sqrt(Z_Bins));
					    Canvas_As_z.divide(XCanvasZ, YCanvasZ);Canvas_FF_z.divide(XCanvasZ, YCanvasZ);
					    for(int k=0; k<Z_Bins;k++) {
					    	Canvas_As_z.cd(k);Canvas_FF_z.cd(k);
							GraphErrors BSAvsz = new GraphErrors();BSAvsz.setMarkerColor(3);BSAvsz.setTitle("z Bin:"+k);;BSAvsz.setTitleX("Pt [GeV]");BSAvsz.setTitleY("ALU");
							GraphErrors FFvsz = new GraphErrors();FFvsz.setMarkerColor(4);FFvsz.setTitle("z Bin:"+k);FFvsz.setTitleX("Pt [GeV]");FFvsz.setTitleY("FLU/FUU");
							for(int l=1;l<PT_Bins-1;l++) {
								BSAvsz.addPoint(Pt_Values[0][2][k][l], BSA_Values[0][2][k][l], 0, BSA_Values[0][2][k][l]);
						    	FFvsz.addPoint(Pt_Values[0][2][k][l], Str_F_Values[0][2][k][l], 0, Str_F_errors[0][2][k][l]);
							}
							 Canvas_As_z.draw(BSAvsz);Canvas_FF_z.draw(FFvsz);
						}
					    Canvas_As_z.save(Newfolder+"/BSA_fixedZ_vs_Pt.png");
						Canvas_FF_z.save(Newfolder+"/FF_fixedZ_vs_Pt.png");
					}
					}
			
			if(b==5) {
				String Newfolder= Outputdir+"/ElBin_"+b+"/";
				 File folder = new File(Newfolder);
					folder.mkdir();
					File File_txt = new File(Newfolder+"/Bin_"+b+"_4D_z_PT_table.txt");
					try {
						File_txt.createNewFile();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					if(File_txt.exists()) {
						System.out.println(" Output File Found" );
						PrintWriter outP = new PrintWriter(File_txt);
						NumberFormat nf = NumberFormat.getInstance();
				        nf.setMinimumFractionDigits(4);
						outP.println("======================================================================================================================================");
					outP.println("====================================     2D Q and X BSA Analysis from Giovanni Angelini code       =======================================");
				    outP.println(" ====== All the Chisquared/NDF are to be considered having an error of +/- 0.408 since 12 data points are used for the fits =======");
					outP.println("===================================================================================================================================");
					outP.println("|<Q2> | <xB>| <W> | <y> | <z> | <Pt> | <epsilon> | Red.ChiSqr | ALU | Err_ALU | FLU/FUU | Err_FLU/UU | Pions Nr |");
					outP.println("===================================================================================================================================");
				for(int l=0; l<PT_Bins;l++) {
					for(int k=0;k<Z_Bins;k++) {
						outP.println(nf.format(Q2_Values[1][2][k][l])+" "+nf.format(xB_Values[1][2][k][l])+" "+nf.format(W_Values[1][2][k][l])+" "+nf.format(y_Values[1][2][k][l])+" "+nf.format(z_Values[1][2][k][l])+" "+nf.format(Pt_Values[1][2][k][l])+" "+nf.format(eps_Values[1][2][k][l])+" "+
								nf.format(RChiSquared[1][2][k][l])+" "+nf.format(BSA_Values[1][2][k][l]) + " "+nf.format(BSA_Errors[1][2][k][l])+" "+nf.format(Str_F_Values[1][2][k][l])+" "+nf.format(Str_F_errors[1][2][k][l])+" "+nf.format(TPion_helicity[1][2][k][l]) );
					}
				}
				outP.close();
				TGCanvas Canvas_As_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				TGCanvas Canvas_FF_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				int XCanvas=((int) Math.floor(Math.sqrt(PT_Bins)));
				    int YCanvas=(int) Math.floor(Math.sqrt(PT_Bins));
				    Canvas_As_PT.divide(XCanvas+1, YCanvas+1);Canvas_FF_PT.divide(XCanvas+1, YCanvas+1);
				    for(int l=0; l<PT_Bins;l++) {
				    	Canvas_As_PT.cd(l);Canvas_FF_PT.cd(l);
						GraphErrors BSAvsPt = new GraphErrors();BSAvsPt.setMarkerColor(1);BSAvsPt.setTitle("PT Bin:"+l);;BSAvsPt.setTitleX(" z ");BSAvsPt.setTitleY("ALU");
						GraphErrors FFvsPt = new GraphErrors();FFvsPt.setMarkerColor(2);FFvsPt.setTitle("PT Bin:"+l );FFvsPt.setTitleX(" z ");FFvsPt.setTitleY("FLU/FUU");
						for(int k=2;k<Z_Bins;k++) {
						
							BSAvsPt.addPoint(z_Values[1][2][k][l], BSA_Values[1][2][k][l], 0, BSA_Values[1][2][k][l]);
					    	FFvsPt.addPoint(z_Values[1][2][k][l], Str_F_Values[1][2][k][l], 0, Str_F_errors[1][2][k][l]);
						}
						 Canvas_As_PT.draw(BSAvsPt);Canvas_FF_PT.draw(FFvsPt);
					}
				    Canvas_As_PT.save(Newfolder+"/BSA_fixedPt_vs_Z.png");
					Canvas_FF_PT.save(Newfolder+"/FF_fixedPt_vs_Z.png");
					
					TGCanvas Canvas_As_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
					TGCanvas Canvas_FF_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
					int XCanvasZ=((int) Math.floor(Math.sqrt(Z_Bins)));
					    int YCanvasZ=(int) Math.floor(Math.sqrt(Z_Bins));
					    Canvas_As_z.divide(XCanvasZ, YCanvasZ);Canvas_FF_z.divide(XCanvasZ, YCanvasZ);
					    for(int k=0; k<Z_Bins;k++) {
					    	Canvas_As_z.cd(k);Canvas_FF_z.cd(k);
							GraphErrors BSAvsz = new GraphErrors();BSAvsz.setMarkerColor(3);BSAvsz.setTitle("z Bin:"+k);;BSAvsz.setTitleX("Pt [GeV]");BSAvsz.setTitleY("ALU");
							GraphErrors FFvsz = new GraphErrors();FFvsz.setMarkerColor(4);FFvsz.setTitle("z Bin:"+k);FFvsz.setTitleX("Pt [GeV]");FFvsz.setTitleY("FLU/FUU");
							for(int l=1;l<PT_Bins-1;l++) {
								BSAvsz.addPoint(Pt_Values[1][2][k][l], BSA_Values[1][2][k][l], 0, BSA_Values[1][2][k][l]);
						    	FFvsz.addPoint(Pt_Values[1][2][k][l], Str_F_Values[1][2][k][l], 0, Str_F_errors[1][2][k][l]);
							}
							 Canvas_As_z.draw(BSAvsz);Canvas_FF_z.draw(FFvsz);
						}
					    Canvas_As_z.save(Newfolder+"/BSA_fixedZ_vs_Pt.png");
						Canvas_FF_z.save(Newfolder+"/FF_fixedZ_vs_Pt.png");
					}
					}
			if(b==6) {
				String Newfolder= Outputdir+"/ElBin_"+b+"/";
				 File folder = new File(Newfolder);
					folder.mkdir();
					File File_txt = new File(Newfolder+"/Bin_"+b+"_4D_z_PT_table.txt");
					try {
						File_txt.createNewFile();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					if(File_txt.exists()) {
						System.out.println(" Output File Found" );
						PrintWriter outP = new PrintWriter(File_txt);
						NumberFormat nf = NumberFormat.getInstance();
				        nf.setMinimumFractionDigits(4);
						outP.println("======================================================================================================================================");
					outP.println("====================================     2D Q and X BSA Analysis from Giovanni Angelini code       =======================================");
				    outP.println(" ====== All the Chisquared/NDF are to be considered having an error of +/- 0.408 since 12 data points are used for the fits =======");
					outP.println("===================================================================================================================================");
					outP.println("|<Q2> | <xB>| <W> | <y> | <z> | <Pt> | <epsilon> | Red.ChiSqr | ALU | Err_ALU | FLU/FUU | Err_FLU/UU | Pions Nr |");
					outP.println("===================================================================================================================================");
				for(int l=0; l<PT_Bins;l++) {
					for(int k=0;k<Z_Bins;k++) {
						outP.println(nf.format(Q2_Values[0][3][k][l])+" "+nf.format(xB_Values[0][3][k][l])+" "+nf.format(W_Values[0][3][k][l])+" "+nf.format(y_Values[0][3][k][l])+" "+nf.format(z_Values[0][3][k][l])+" "+nf.format(Pt_Values[0][3][k][l])+" "+nf.format(eps_Values[0][3][k][l])+" "+
								nf.format(RChiSquared[0][3][k][l])+" "+nf.format(BSA_Values[0][3][k][l]) + " "+nf.format(BSA_Errors[0][3][k][l])+" "+nf.format(Str_F_Values[0][3][k][l])+" "+nf.format(Str_F_errors[0][3][k][l])+" "+nf.format(TPion_helicity[0][3][k][l]) );
					}
				}
				outP.close();
				TGCanvas Canvas_As_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				TGCanvas Canvas_FF_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				int XCanvas=((int) Math.floor(Math.sqrt(PT_Bins)));
				    int YCanvas=(int) Math.floor(Math.sqrt(PT_Bins));
				    Canvas_As_PT.divide(XCanvas+1, YCanvas+1);Canvas_FF_PT.divide(XCanvas+1, YCanvas+1);
				    for(int l=0; l<PT_Bins;l++) {
				    	Canvas_As_PT.cd(l);Canvas_FF_PT.cd(l);
						GraphErrors BSAvsPt = new GraphErrors();BSAvsPt.setMarkerColor(1);BSAvsPt.setTitle("PT Bin:"+l);;BSAvsPt.setTitleX(" z ");BSAvsPt.setTitleY("ALU");
						GraphErrors FFvsPt = new GraphErrors();FFvsPt.setMarkerColor(2);FFvsPt.setTitle("PT Bin:"+l );FFvsPt.setTitleX(" z ");FFvsPt.setTitleY("FLU/FUU");
						for(int k=2;k<Z_Bins;k++) {
						
							BSAvsPt.addPoint(z_Values[0][3][k][l], BSA_Values[0][3][k][l], 0, BSA_Values[0][3][k][l]);
					    	FFvsPt.addPoint(z_Values[0][3][k][l], Str_F_Values[0][3][k][l], 0, Str_F_errors[0][3][k][l]);
						}
						 Canvas_As_PT.draw(BSAvsPt);Canvas_FF_PT.draw(FFvsPt);
					}
				    Canvas_As_PT.save(Newfolder+"/BSA_fixedPt_vs_Z.png");
					Canvas_FF_PT.save(Newfolder+"/FF_fixedPt_vs_Z.png");
					
					TGCanvas Canvas_As_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
					TGCanvas Canvas_FF_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
					int XCanvasZ=((int) Math.floor(Math.sqrt(Z_Bins)));
					    int YCanvasZ=(int) Math.floor(Math.sqrt(Z_Bins));
					    Canvas_As_z.divide(XCanvasZ, YCanvasZ);Canvas_FF_z.divide(XCanvasZ, YCanvasZ);
					    for(int k=0; k<Z_Bins;k++) {
					    	Canvas_As_z.cd(k);Canvas_FF_z.cd(k);
							GraphErrors BSAvsz = new GraphErrors();BSAvsz.setMarkerColor(3);BSAvsz.setTitle("z Bin:"+k);;BSAvsz.setTitleX("Pt [GeV]");BSAvsz.setTitleY("ALU");
							GraphErrors FFvsz = new GraphErrors();FFvsz.setMarkerColor(4);FFvsz.setTitle("z Bin:"+k);FFvsz.setTitleX("Pt [GeV]");FFvsz.setTitleY("FLU/FUU");
							for(int l=1;l<PT_Bins-1;l++) {
								BSAvsz.addPoint(Pt_Values[0][3][k][l], BSA_Values[0][3][k][l], 0, BSA_Values[0][3][k][l]);
						    	FFvsz.addPoint(Pt_Values[0][3][k][l], Str_F_Values[0][3][k][l], 0, Str_F_errors[0][3][k][l]);
							}
							 Canvas_As_z.draw(BSAvsz);Canvas_FF_z.draw(FFvsz);
						}
					    Canvas_As_z.save(Newfolder+"/BSA_fixedZ_vs_Pt.png");
						Canvas_FF_z.save(Newfolder+"/FF_fixedZ_vs_Pt.png");
					}
					}
			
			if(b==7) {
				String Newfolder= Outputdir+"/ElBin_"+b+"/";
				 File folder = new File(Newfolder);
					folder.mkdir();
					File File_txt = new File(Newfolder+"/Bin_"+b+"_4D_z_PT_table.txt");
					try {
						File_txt.createNewFile();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					if(File_txt.exists()) {
						System.out.println(" Output File Found" );
						PrintWriter outP = new PrintWriter(File_txt);
						NumberFormat nf = NumberFormat.getInstance();
				        nf.setMinimumFractionDigits(4);
						outP.println("======================================================================================================================================");
					outP.println("====================================     2D Q and X BSA Analysis from Giovanni Angelini code       =======================================");
				    outP.println(" ====== All the Chisquared/NDF are to be considered having an error of +/- 0.408 since 12 data points are used for the fits =======");
					outP.println("===================================================================================================================================");
					outP.println("|<Q2> | <xB>| <W> | <y> | <z> | <Pt> | <epsilon> | Red.ChiSqr | ALU | Err_ALU | FLU/FUU | Err_FLU/UU | Pions Nr |");
					outP.println("===================================================================================================================================");
				for(int l=0; l<PT_Bins;l++) {
					for(int k=0;k<Z_Bins;k++) {
						outP.println(nf.format(Q2_Values[1][3][k][l])+" "+nf.format(xB_Values[1][3][k][l])+" "+nf.format(W_Values[1][3][k][l])+" "+nf.format(y_Values[1][3][k][l])+" "+nf.format(z_Values[1][3][k][l])+" "+nf.format(Pt_Values[1][3][k][l])+" "+nf.format(eps_Values[1][3][k][l])+" "+
								nf.format(RChiSquared[1][3][k][l])+" "+nf.format(BSA_Values[1][3][k][l]) + " "+nf.format(BSA_Errors[1][3][k][l])+" "+nf.format(Str_F_Values[1][3][k][l])+" "+nf.format(Str_F_errors[1][3][k][l])+" "+nf.format(TPion_helicity[1][3][k][l]) );
					}
				}
				outP.close();
				TGCanvas Canvas_As_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				TGCanvas Canvas_FF_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				int XCanvas=((int) Math.floor(Math.sqrt(PT_Bins)));
				    int YCanvas=(int) Math.floor(Math.sqrt(PT_Bins));
				    Canvas_As_PT.divide(XCanvas+1, YCanvas+1);Canvas_FF_PT.divide(XCanvas+1, YCanvas+1);
				    for(int l=0; l<PT_Bins;l++) {
				    	Canvas_As_PT.cd(l);Canvas_FF_PT.cd(l);
						GraphErrors BSAvsPt = new GraphErrors();BSAvsPt.setMarkerColor(1);BSAvsPt.setTitle("PT Bin:"+l);;BSAvsPt.setTitleX(" z ");BSAvsPt.setTitleY("ALU");
						GraphErrors FFvsPt = new GraphErrors();FFvsPt.setMarkerColor(2);FFvsPt.setTitle("PT Bin:"+l );FFvsPt.setTitleX(" z ");FFvsPt.setTitleY("FLU/FUU");
						for(int k=2;k<Z_Bins;k++) {
						
							BSAvsPt.addPoint(z_Values[1][3][k][l], BSA_Values[1][3][k][l], 0, BSA_Values[1][3][k][l]);
					    	FFvsPt.addPoint(z_Values[1][3][k][l], Str_F_Values[1][3][k][l], 0, Str_F_errors[1][3][k][l]);
						}
						 Canvas_As_PT.draw(BSAvsPt);Canvas_FF_PT.draw(FFvsPt);
					}
				    Canvas_As_PT.save(Newfolder+"/BSA_fixedPt_vs_Z.png");
					Canvas_FF_PT.save(Newfolder+"/FF_fixedPt_vs_Z.png");
					
					TGCanvas Canvas_As_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
					TGCanvas Canvas_FF_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
					int XCanvasZ=((int) Math.floor(Math.sqrt(Z_Bins)));
					    int YCanvasZ=(int) Math.floor(Math.sqrt(Z_Bins));
					    Canvas_As_z.divide(XCanvasZ, YCanvasZ);Canvas_FF_z.divide(XCanvasZ, YCanvasZ);
					    for(int k=0; k<Z_Bins;k++) {
					    	Canvas_As_z.cd(k);Canvas_FF_z.cd(k);
							GraphErrors BSAvsz = new GraphErrors();BSAvsz.setMarkerColor(3);BSAvsz.setTitle("z Bin:"+k);;BSAvsz.setTitleX("Pt [GeV]");BSAvsz.setTitleY("ALU");
							GraphErrors FFvsz = new GraphErrors();FFvsz.setMarkerColor(4);FFvsz.setTitle("z Bin:"+k);FFvsz.setTitleX("Pt [GeV]");FFvsz.setTitleY("FLU/FUU");
							for(int l=1;l<PT_Bins-1;l++) {
								BSAvsz.addPoint(Pt_Values[1][3][k][l], BSA_Values[1][3][k][l], 0, BSA_Values[1][3][k][l]);
						    	FFvsz.addPoint(Pt_Values[1][3][k][l], Str_F_Values[1][3][k][l], 0, Str_F_errors[1][3][k][l]);
							}
							 Canvas_As_z.draw(BSAvsz);Canvas_FF_z.draw(FFvsz);
						}
					    Canvas_As_z.save(Newfolder+"/BSA_fixedZ_vs_Pt.png");
						Canvas_FF_z.save(Newfolder+"/FF_fixedZ_vs_Pt.png");
					}
					}
		
			
			if(b==8) {
				String Newfolder= Outputdir+"/ElBin_"+b+"/";
				 File folder = new File(Newfolder);
					folder.mkdir();
					File File_txt = new File(Newfolder+"/Bin_"+b+"_4D_z_PT_table.txt");
					try {
						File_txt.createNewFile();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					if(File_txt.exists()) {
						System.out.println(" Output File Found" );
						PrintWriter outP = new PrintWriter(File_txt);
						NumberFormat nf = NumberFormat.getInstance();
				        nf.setMinimumFractionDigits(4);
						outP.println("======================================================================================================================================");
					outP.println("====================================     2D Q and X BSA Analysis from Giovanni Angelini code       =======================================");
				    outP.println(" ====== All the Chisquared/NDF are to be considered having an error of +/- 0.408 since 12 data points are used for the fits =======");
					outP.println("===================================================================================================================================");
					outP.println("|<Q2> | <xB>| <W> | <y> | <z> | <Pt> | <epsilon> | Red.ChiSqr | ALU | Err_ALU | FLU/FUU | Err_FLU/UU | Pions Nr |");
					outP.println("===================================================================================================================================");
				for(int l=0; l<PT_Bins;l++) {
					
					for(int k=0;k<Z_Bins;k++) {
						outP.println(nf.format(Q2_Values[0][4][k][l])+" "+nf.format(xB_Values[0][4][k][l])+" "+nf.format(W_Values[0][4][k][l])+" "+nf.format(y_Values[0][4][k][l])+" "+nf.format(z_Values[0][4][k][l])+" "+nf.format(Pt_Values[0][4][k][l])+" "+nf.format(eps_Values[0][4][k][l])+" "+
								nf.format(RChiSquared[0][4][k][l])+" "+nf.format(BSA_Values[0][4][k][l]) + " "+nf.format(BSA_Errors[0][4][k][l])+" "+nf.format(Str_F_Values[0][4][k][l])+" "+nf.format(Str_F_errors[0][4][k][l])+" "+nf.format(TPion_helicity[0][4][k][l]) );
					}
				}
				outP.close();
				TGCanvas Canvas_As_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				TGCanvas Canvas_FF_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				int XCanvas=((int) Math.floor(Math.sqrt(PT_Bins)));
				    int YCanvas=(int) Math.floor(Math.sqrt(PT_Bins));
				    Canvas_As_PT.divide(XCanvas+1, YCanvas+1);Canvas_FF_PT.divide(XCanvas+1, YCanvas+1);
				    for(int l=0; l<PT_Bins;l++) {
				    	Canvas_As_PT.cd(l);Canvas_FF_PT.cd(l);
						GraphErrors BSAvsPt = new GraphErrors();BSAvsPt.setMarkerColor(1);BSAvsPt.setTitle("PT Bin:"+l);;BSAvsPt.setTitleX(" z ");BSAvsPt.setTitleY("ALU");
						GraphErrors FFvsPt = new GraphErrors();FFvsPt.setMarkerColor(2);FFvsPt.setTitle("PT Bin:"+l );FFvsPt.setTitleX(" z ");FFvsPt.setTitleY("FLU/FUU");
						for(int k=2;k<Z_Bins;k++) {
						
							BSAvsPt.addPoint(z_Values[0][4][k][l], BSA_Values[0][4][k][l], 0, BSA_Values[0][4][k][l]);
					    	FFvsPt.addPoint(z_Values[0][4][k][l], Str_F_Values[0][4][k][l], 0, Str_F_errors[0][4][k][l]);
						}
						 Canvas_As_PT.draw(BSAvsPt);Canvas_FF_PT.draw(FFvsPt);
					}
				    Canvas_As_PT.save(Newfolder+"/BSA_fixedPt_vs_Z.png");
					Canvas_FF_PT.save(Newfolder+"/FF_fixedPt_vs_Z.png");
					
					TGCanvas Canvas_As_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
					TGCanvas Canvas_FF_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
					int XCanvasZ=((int) Math.floor(Math.sqrt(Z_Bins)));
					    int YCanvasZ=(int) Math.floor(Math.sqrt(Z_Bins));
					    Canvas_As_z.divide(XCanvasZ, YCanvasZ);Canvas_FF_z.divide(XCanvasZ, YCanvasZ);
					    for(int k=0; k<Z_Bins;k++) {
					    	Canvas_As_z.cd(k);Canvas_FF_z.cd(k);
							GraphErrors BSAvsz = new GraphErrors();BSAvsz.setMarkerColor(3);BSAvsz.setTitle("z Bin:"+k);;BSAvsz.setTitleX("Pt [GeV]");BSAvsz.setTitleY("ALU");
							GraphErrors FFvsz = new GraphErrors();FFvsz.setMarkerColor(4);FFvsz.setTitle("z Bin:"+k);FFvsz.setTitleX("Pt [GeV]");FFvsz.setTitleY("FLU/FUU");
							for(int l=1;l<PT_Bins-1;l++) {
								BSAvsz.addPoint(Pt_Values[0][4][k][l], BSA_Values[0][4][k][l], 0, BSA_Values[0][4][k][l]);
						    	FFvsz.addPoint(Pt_Values[0][4][k][l], Str_F_Values[0][4][k][l], 0, Str_F_errors[0][4][k][l]);
							}
							 Canvas_As_z.draw(BSAvsz);Canvas_FF_z.draw(FFvsz);
						}
					    Canvas_As_z.save(Newfolder+"/BSA_fixedZ_vs_Pt.png");
						Canvas_FF_z.save(Newfolder+"/FF_fixedZ_vs_Pt.png");
					}
					}
				
			
			if(b==9) {
				String Newfolder= Outputdir+"/ElBin_"+b+"/";
				 File folder = new File(Newfolder);
					folder.mkdir();
					File File_txt = new File(Newfolder+"/Bin_"+b+"_4D_z_PT_table.txt");
					try {
						File_txt.createNewFile();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					if(File_txt.exists()) {
						System.out.println(" Output File Found" );
						PrintWriter outP = new PrintWriter(File_txt);
						NumberFormat nf = NumberFormat.getInstance();
				        nf.setMinimumFractionDigits(4);
						outP.println("======================================================================================================================================");
					outP.println("====================================     2D Q and X BSA Analysis from Giovanni Angelini code       =======================================");
				    outP.println(" ====== All the Chisquared/NDF are to be considered having an error of +/- 0.408 since 12 data points are used for the fits =======");
					outP.println("===================================================================================================================================");
					outP.println("|<Q2> | <xB>| <W> | <y> | <z> | <Pt> | <epsilon> | Red.ChiSqr | ALU | Err_ALU | FLU/FUU | Err_FLU/UU | Pions Nr |");
					outP.println("===================================================================================================================================");
				for(int l=0; l<PT_Bins;l++) {
					for(int k=0;k<Z_Bins;k++) {
						outP.println(nf.format(Q2_Values[1][4][k][l])+" "+nf.format(xB_Values[1][4][k][l])+" "+nf.format(W_Values[1][4][k][l])+" "+nf.format(y_Values[1][4][k][l])+" "+nf.format(z_Values[1][4][k][l])+" "+nf.format(Pt_Values[1][4][k][l])+" "+nf.format(eps_Values[1][4][k][l])+" "+
								nf.format(RChiSquared[1][4][k][l])+" "+nf.format(BSA_Values[1][4][k][l]) + " "+nf.format(BSA_Errors[1][4][k][l])+" "+nf.format(Str_F_Values[1][4][k][l])+" "+nf.format(Str_F_errors[1][4][k][l])+" "+nf.format(TPion_helicity[1][4][k][l]) );
					}
				}
				outP.close();
				TGCanvas Canvas_As_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				TGCanvas Canvas_FF_PT = new TGCanvas("canvas_AS_PT","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
				int XCanvas=((int) Math.floor(Math.sqrt(PT_Bins)));
				    int YCanvas=(int) Math.floor(Math.sqrt(PT_Bins));
				    Canvas_As_PT.divide(XCanvas+1, YCanvas+1);Canvas_FF_PT.divide(XCanvas+1, YCanvas+1);
				    for(int l=0; l<PT_Bins;l++) {
				    	Canvas_As_PT.cd(l);Canvas_FF_PT.cd(l);
						GraphErrors BSAvsPt = new GraphErrors();BSAvsPt.setMarkerColor(1);BSAvsPt.setTitle("PT Bin:"+l);;BSAvsPt.setTitleX(" z ");BSAvsPt.setTitleY("ALU");
						GraphErrors FFvsPt = new GraphErrors();FFvsPt.setMarkerColor(2);FFvsPt.setTitle("PT Bin:"+l );FFvsPt.setTitleX(" z ");FFvsPt.setTitleY("FLU/FUU");
						for(int k=2;k<Z_Bins;k++) {
						
							BSAvsPt.addPoint(z_Values[1][4][k][l], BSA_Values[1][4][k][l], 0, BSA_Values[1][4][k][l]);
					    	FFvsPt.addPoint(z_Values[1][4][k][l], Str_F_Values[1][4][k][l], 0, Str_F_errors[1][4][k][l]);
						}
						 Canvas_As_PT.draw(BSAvsPt);Canvas_FF_PT.draw(FFvsPt);
					}
				    Canvas_As_PT.save(Newfolder+"/BSA_fixedPt_vs_Z.png");
					Canvas_FF_PT.save(Newfolder+"/FF_fixedPt_vs_Z.png");
					
					TGCanvas Canvas_As_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
					TGCanvas Canvas_FF_z = new TGCanvas("canvas_AS_z","canvas_AS_PT",1600,1400) ; // Create Canvas with (1,1) divisions
					int XCanvasZ=((int) Math.floor(Math.sqrt(Z_Bins)));
					    int YCanvasZ=(int) Math.floor(Math.sqrt(Z_Bins));
					    Canvas_As_z.divide(XCanvasZ, YCanvasZ);Canvas_FF_z.divide(XCanvasZ, YCanvasZ);
					    for(int k=0; k<Z_Bins;k++) {
					    	Canvas_As_z.cd(k);Canvas_FF_z.cd(k);
							GraphErrors BSAvsz = new GraphErrors();BSAvsz.setMarkerColor(3);BSAvsz.setTitle("z Bin:"+k);;BSAvsz.setTitleX("Pt [GeV]");BSAvsz.setTitleY("ALU");
							GraphErrors FFvsz = new GraphErrors();FFvsz.setMarkerColor(4);FFvsz.setTitle("z Bin:"+k);FFvsz.setTitleX("Pt [GeV]");FFvsz.setTitleY("FLU/FUU");
							for(int l=1;l<PT_Bins-1;l++) {
								BSAvsz.addPoint(Pt_Values[1][4][k][l], BSA_Values[1][4][k][l], 0, BSA_Values[1][4][k][l]);
						    	FFvsz.addPoint(Pt_Values[1][4][k][l], Str_F_Values[1][4][k][l], 0, Str_F_errors[1][4][k][l]);
							}
							 Canvas_As_z.draw(BSAvsz);Canvas_FF_z.draw(FFvsz);
						}
					    Canvas_As_z.save(Newfolder+"/BSA_fixedZ_vs_Pt.png");
						Canvas_FF_z.save(Newfolder+"/FF_fixedZ_vs_Pt.png");
					}
					}
			
		}
	}
	
	
}
