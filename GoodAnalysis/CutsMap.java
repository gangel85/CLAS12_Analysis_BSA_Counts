package GoodAnalysis;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jlab.groot.data.H2F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.graphics.EmbeddedCanvas;



public class CutsMap {
	private Map<Integer, List<DetectorCut>> Map = new HashMap<Integer, List<DetectorCut>>();
	/* --- Constructors
	   ---
	 */
	CutsMap(){
		
		// Generating the  standard cuts 
		EC_Cut FC_PCAL_el = new EC_Cut();
		FC_PCAL_el.setName("PCAL_el");
		EC_Cut FC_EC1_el = new EC_Cut();
		
		DC_Cut FC_DC1_el = new DC_Cut();
		FC_DC1_el.setName("DC1_el");
	
		Sampling_Fractions SF_el = new Sampling_Fractions(); 
		SF_el.setName("SamplingFraction_el");
		
		DC_Cut FC_DC3_el = new DC_Cut();
		FC_DC3_el.setName("DC3_el");
		EC_Cut FC_PCAL_ph = new EC_Cut();
		FC_PCAL_ph.setName("PCAL_ph");
		EC_Cut FC_EC1_ph = new EC_Cut();
		FC_EC1_ph.setName("EC1_ph");
		EC_Cut FC_EC2_ph = new EC_Cut();
		FC_EC2_ph.setName("EC2_ph");
		
		
		// Defining the list of cuts  for electron
		List<DetectorCut> ElectronCuts = new ArrayList<DetectorCut>();
		
		
		
		ElectronCuts.add(FC_PCAL_el);
		ElectronCuts.add(FC_DC1_el);
		ElectronCuts.add(SF_el);
		//ElectronCuts.add(FC_DC2_el);
	//	ElectronCuts.add(FC_DC3_el);
	   //ElectronCuts.add(FC_EC1_el);
	//	ElectronCuts.add(FC_EC2_el);
		Map.put(11, ElectronCuts); 
		// Defining the list of cuts for photon
	//	List<DetectorCut> PhotonCuts = new ArrayList<DetectorCut>();
		//PhotonCuts.add(FC_PCAL_ph);
		//PhotonCuts.add(FC_EC1_ph);
		//ElectronCuts.PhotonCuts.add(FC_EC2_ph);
		//Map.put(22,PhotonCuts);
		
		// Cuts for PiP
		DC_Cut FC_DC_pip = new DC_Cut();
		FC_DC_pip.setName("DC_pip");
		List<DetectorCut> PionPlusCuts = new ArrayList<DetectorCut>();
		PionPlusCuts.add(FC_DC_pip);
		Map.put(211,PionPlusCuts);
		
		// Cuts for PiM
		//DC_Cut FC_DC_pim = new DC_Cut();
		//FC_DC_pim.setName("DC_pim");
		//List<DetectorCut> PionMinusCuts = new ArrayList<DetectorCut>();
		//PionMinusCuts.add(FC_DC_pim);
		//Map.put(-211,PionMinusCuts);
		
		
	}
	// Uses user defined map to create the CutsMap
	CutsMap(Map<? extends Integer, ? extends List<DetectorCut>> UserMap){
		Map.clear();
		Map.putAll(UserMap);
	}	

	/**
	 * This function print the Map of PID and DetectorCut loaded in the CutsMap object
	 */
	public void show() {
		for (Integer name: Map.keySet()) { System.out.println(" PID: " + name);
		for(int ii=0;ii< Map.get(name).size(); ii++) {
			System.out.print(" | " + Map.get(name).get(ii).toString());
		}
		System.out.println(" ");
		}
	}
	

	/** Allow the user to substitute the pre-existing Map
	 * The user should provide a map in the form: Key Int, Value is a List<DetectorCut>
	 */
	public void	SetMap( Map<Integer, List<DetectorCut>> UserMap){
		Map.clear();
		Map.putAll(UserMap);
	}
	/**
	 * Add an entry in the Map 
	 * @param i the PID (Key of the map)
	 * @param CutsList a List<DetectorCut> representing the value in the map
	 */
	public CutsMap add(int i, List<DetectorCut> CutsList) {
		Map.put(i, CutsList); return this;
	}
	
	public static CutsMap getDefault() {
		CutsMap map = new CutsMap();
		return map;
	}
	//
	/**	
	 * 
	 * @return the Map<Integer, List<DetectorCut>> saved in the Objec CutsMap
	 */
	public Map GetMap() {
		return Map;
	}

	/**
	 *  Get the List<Detector> associated with a PID 
	 * @param pid is the key for the map (PID)
	 * @return the Value associated with a PID (null if not found)
	 */
	public List<DetectorCut> getPID(int pid) {
		for (Integer id: Map.keySet()) {
			if (id == pid) return Map.get(pid);
		}
		return new ArrayList<DetectorCut>(); //if you dont find that pid 
	}
	
	//Aggiugnere i printout delle sampling fractions qui. 
	public void Print(String time2, String workdirout) {
		// For all the Particles 
		for(Integer id: Map.keySet()) {
			// For all Detector I have cut on
			for( int i =0 ; i< Map.get(id).size(); i++) {
			// Get the Histograms 
				String Particella =id.toString() ;
				//this is a good one to study
			List<H2F> Histos_Pre=Map.get(id).get(i).Histograms_Pre();
			List<H2F> Histos_Aft=Map.get(id).get(i).Histograms_Aft();
			String s1=workdirout;
			new File(s1+"/Plots_FiducialCuts").mkdir();
			TDirectory directoryMain = new TDirectory();	
			String PathFolder = s1+"/Plots_FiducialCuts/"+time2;
			File dir = new File(PathFolder);dir.mkdir();
			new File (PathFolder).mkdir();
			System.out.println(" (FIDUCIAL CUTS) - > I am Plotting Multi Bins Histograms");
			directoryMain.mkdir("/main/");directoryMain.cd("/main/");
			EmbeddedCanvas FiducialCut_Pre = new EmbeddedCanvas();
			EmbeddedCanvas FiducialCut_Aft = new EmbeddedCanvas();
			FiducialCut_Pre.setSize(1600,1000);FiducialCut_Aft.setSize(1600,1000);
			FiducialCut_Pre.divide(3,3);FiducialCut_Aft.divide(3,3);
			FiducialCut_Pre.setAxisTitleSize(22); FiducialCut_Aft.setAxisTitleSize(22); 
			FiducialCut_Pre.setAxisFontSize(22);FiducialCut_Aft.setAxisFontSize(22);
			FiducialCut_Pre.setTitleSize(22);FiducialCut_Aft.setTitleSize(22);		
			for(int k=0; k<Histos_Pre.size() ; k++) {
			FiducialCut_Pre.cd(k);FiducialCut_Pre.draw(Histos_Pre.get(k));
			directoryMain.addDataSet(Histos_Pre.get(k));
			}
			 for(int k=0; k<Histos_Aft.size() ; k++) {
				FiducialCut_Aft.cd(k);FiducialCut_Aft.draw(Histos_Aft.get(k));
				directoryMain.addDataSet(Histos_Aft.get(k));
				}
			directoryMain.writeFile(PathFolder+"/"+Particella+"_"+i+"_FiducialCuts_Pre.hipo");
			String strg0 = String.format("%s/"+id+"_"+i+"_FiducialCut_Pre.png",PathFolder);
			FiducialCut_Pre.save(strg0);
		  System.out.println(" The particle ID is " + id +" Partiella is " + Particella);
			directoryMain.writeFile(PathFolder+"/"+Particella+"_"+i+"_FiducialCuts_Aft.hipo");
			String strg1 = String.format("%s/"+id+"_"+i+"FiducialCut_Aft.png",PathFolder);
			
			FiducialCut_Aft.save(strg1);
			}
		}	
	}
}
