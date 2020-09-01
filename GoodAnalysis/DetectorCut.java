/**
 * 
 */
package GoodAnalysis;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jlab.groot.data.H2F;
import org.jlab.jnp.physics.Vector3;

/**
 * @author gangelini
 * This interface define a generic Detector Cut
 */
public interface DetectorCut {
 boolean Status(ParticleREC particle);
 int getSector(Vector3 hit);
void setName(String nameDetector);
List<H2F> Histograms_Aft();
String toString();
//	Map<Integer, List<String>> CutsMap = new HashMap<Integer, List<String>>();
	//List<String> CutsList = new ArrayList<String>();
List<H2F> Histograms_Pre();
 
}
