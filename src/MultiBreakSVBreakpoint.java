

import java.util.Arrays;


public class MultiBreakSVBreakpoint {
	public int id;
	public double[] supportsshort,supportslong,supportsthinned,supportslast,supportsfirst;
	public String cID;
	public boolean isCluster;
	
	public MultiBreakSVBreakpoint(int i,String c,boolean isc,int maxSupport) {
		id = i;
		supportsshort = new double[maxSupport+1];
		supportslong = new double[maxSupport+1];
		supportsthinned = new double[maxSupport+1];
		supportslast = new double[maxSupport+1];		
		supportsfirst = new double[maxSupport+1];
		Arrays.fill(supportsshort,0);
		Arrays.fill(supportslong,0);
		Arrays.fill(supportsthinned,0);
		Arrays.fill(supportslast,0);
		Arrays.fill(supportsfirst,0);
		cID = c;
		isCluster = isc;
	}
	
}
