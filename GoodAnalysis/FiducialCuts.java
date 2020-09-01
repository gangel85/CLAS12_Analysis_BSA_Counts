package GoodAnalysis;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FiducialCuts {
	CutsMap FC_Map; //This is the map used for the cuts 
	private Map<Integer, List<DetectorCut>> Map = new HashMap<Integer, List<DetectorCut>>();
	
	@SuppressWarnings("unchecked")
	public FiducialCuts(){
		this.FC_Map= new CutsMap();
		System.out.println("Map size " + FC_Map.GetMap().size());	
		Map.clear();
		System.out.println(" MAP IS " + FC_Map.GetMap().toString());
		Map.putAll(FC_Map.GetMap());
		}
	
	@SuppressWarnings("unchecked")
	FiducialCuts(CutsMap NewMap)
	{
		this.FC_Map = NewMap;
		Map.clear();
		Map.putAll(FC_Map.GetMap());
	}
	
	/** Return the Map used for the cuts
	 * @author gangel
	 */
  public CutsMap GetMap()
   {
	return this.FC_Map;
   }
   
  /**
   * This function check if a ParticleREC passes all the cuts loaded into the map.
   * @param testParticle is the ParticleREC we want to apply the cut
   * @return the boolean regarding the Cut (true = passed)
   */
	public boolean Status(ParticleREC testParticle) {
		int indicetrue =0;
		for ( Integer pid : Map.keySet()) {
			if(pid == testParticle.pid()) {
				this.FC_Map.getPID(pid);
				for(int i = 0 ; i< this.FC_Map.getPID(pid).size() ; i ++) {
					boolean resultFC=FC_Map.getPID(pid).get(i).Status(testParticle);
		
					if (resultFC==true) indicetrue++;
			    }
				if (indicetrue==this.FC_Map.getPID(pid).size() && this.FC_Map.getPID(pid).size()>0) {
					return true;
				}
				return false;
			}
								}
		return false;
	}

	/**
	 * This function save the histograms associated wqith fiducial cuts
	 * @param time2 is the time stump that will be used to save the folder
	 * @param workdirout is the ouptud dir path
	 */
public void Print(String time2, String workdirout) {
	this.FC_Map.Print(time2,workdirout);
	// TODO Auto-generated method stub
	
}
}