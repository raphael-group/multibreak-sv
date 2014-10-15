

import java.util.HashMap;

public class Experiment {
	
	public double pseq;
	public double lambdaD;
	public String label;
	public double[] logprobcov;
	public HashMap<String,Double> deltaprobs;
	
	public Experiment(double p, double l,String lab) {
		pseq = p;
		lambdaD = l;
		deltaprobs = new HashMap<String,Double>();
		label = lab;
	}
	
}
