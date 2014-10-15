import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import java.util.Scanner;

/**
 * @author aritz
 */
public class MultiBreakSV {

	public String clusterfile,assignmentfile,experimentfile,
		matchfile,outputname;
	public static double perr;
	
	public HashMap<String,MultiBreakSVDiagram> diagrams; // <cID,MultiBreakSVDiagram>
	public static HashMap<String,String> esp2fragmentID; // <esp,fragmentID>
	public static HashMap<String,Integer> deltas; // <cID,delta_k>
	public HashMap<Double,HashMap<String,Double>> deltaprobs; // <lambdaD,<cID,Poiss(delta_k;lambda_d)>>
	public HashMap<String,String> esp2cid; // <ESP,cID>
	public HashMap<String,Integer> counts; // <summarystring,count>
	public HashMap<String,Double> logprobs; // <summarystring,logprob>
	public ArrayList<String> fragmentIDs;
	
	public static HashMap<String,MultiBreakSVBreakpoint> breakpoints; // <cID,breakpoint>
	public static HashMap<String,Fragment> fragments; // <fragmentID,Fragment> 
	public static HashMap<String,Experiment> experiments; // <experimentlabel,Experiment>
	
	public ArrayList<String> bpsWithNonSingletonSupport; 
	public BigInteger spacesize = new BigInteger("1");
	
	// commonly used variables
	public Fragment fragment;
	public FragmentAlignment fa;

	int[] moveCounts;
	int[] proposalCounts;
	HashMap<FragmentAlignment,Integer> alignmentHits;

	public MoveType moveType;
	// original: 123456
	public static Random rand = new Random(123456);
	public double mu = 124;
	public double sigma = 1;
	public double explambda = 10;
	public int burnin;
	public int numiters;
	public int skip;

	// variables to speed up computation.
	public int maxSupport;
	public static double[] prior;
	public static boolean SAVE_SPACE = false;
	public static boolean CLUSTERS_FIRST = false;
	
	// variables for enumerating vs. sampling
	public static boolean ENUMERATE = false;
	public static int MAX_NUM_ALIGN = -1;
	public int BASE = 10;

	// hard-coded variables
	public final boolean silent = true;
	public int intermediateJump = Integer.MAX_VALUE;

	public Iterator<String> string_iter;
	public Iterator<Integer> int_iter;

	public double Alogq,Anewlogq;
	public int prevassign,thisassign;
	public String naivefragment;

	/**
	 * Call with --help for usage.
	 * @param args
	 * @throws IOException
	 * @throws ClassNotFoundException
	 */
	public static void main(String[] args) throws IOException, ClassNotFoundException, Exception {
		long starttime = System.nanoTime();
		MultiBreakSV c2g = new MultiBreakSV(args);

		// parse assignmentfile, experimentfile, and clustersfile
		c2g.getFragmentAlignments(false);

		// construct cluster diagram
		c2g.constructClusterDiagrams(false);

		// precompute priors
		c2g.precomputePriors();

		// fill breakpoint array with nonsingleton clusters.
		c2g.fillNonSingletonArray();
		
		// Finally, run method. 
		if (ENUMERATE) { // enumerate all states.
			c2g.enumerate();
			long endtime = System.nanoTime();
			System.out.println("ENUMERATE TIME "+(endtime-starttime)/1000000000.0);
		} else { // run MCMC method.

		MultiBreakSVAssignmentMatrix A;
		if (SAVE_SPACE) 
			A = c2g.runMCMC(null);
		else 
			A = c2g.runMCMC(c2g.outputname+".mcmc.txt");	
		
		// write vertexIDs
		c2g.writeFinalFiles();

		// write state counts
		if (!SAVE_SPACE)
			c2g.writeStateCounts(A.keys);
		
		long endtime = System.nanoTime();
		System.out.println("SAMPLING TIME "+(endtime-starttime)/1000000000.0);
		}	
	}

	public MultiBreakSV() {}

	public MultiBreakSV(String[] args) {

		// print arguments and determine if help is specified.
		boolean printHelpAndExit = false;
		System.out.print("\nSystem Call:\njava mcmc.MultiBreakSV ");
		for(int i=0;i<args.length;i++) { 
			System.out.print(args[i] + " ");
			if (args[i].equals("--help") || args[i].equals("-h") || args[i].equals("help")) 
				printHelpAndExit = true;
		}
		System.out.println();
		
		// Print help and exit if --help is specified or if no arguments are specified.
		if (args.length == 0 || printHelpAndExit) {
			printUsage();
			System.out.println("Printing help and exiting.");
			System.exit(0);
		}
		
		// Try parsing arguments
		try {
			parseOptions(args);
		} catch (Exception e) { 
			printUsage();
			System.out.println("Error was");
			e.printStackTrace();
			System.exit(0);
		}
		// Get burnin counts
		burnin = numiters/10;
		skip = numiters/1000;

		// if numiters is small, burnin and skip should != 0.
		if(skip == 0)
			skip = 1;
		if(burnin == 0) 
			burnin = 10;

		// Initialize counters, HashMaps, and regular expressions.
		breakpoints = new HashMap<String,MultiBreakSVBreakpoint>();
		fragments = new HashMap<String,Fragment>();
		experiments = new HashMap<String,Experiment>();
		alignmentHits = new HashMap<FragmentAlignment,Integer>();
		esp2fragmentID = new HashMap<String,String>();
			
		moveCounts = new int[4];
		proposalCounts = new int[4];		
		Arrays.fill(moveCounts,0);
		Arrays.fill(proposalCounts,0);
		
		counts = new HashMap<String,Integer>();
		logprobs = new HashMap<String,Double>();
	}

	public void parseOptions(String[] args) throws IllegalArgumentException {

		clusterfile = null;
		assignmentfile = null;
		experimentfile = null;
		numiters = -1;
		perr = -1;
		outputname = "out";
		matchfile = null; 
		
		try {
		for (int i=0;i<args.length;i++) {
			if (args[i].equals("--clusterfile")) { clusterfile = args[i+1]; i++;} 
			else if (args[i].equals("--assignmentfile")) { assignmentfile = args[i+1]; i++;}
			else if (args[i].equals("--experimentfile")) { experimentfile = args[i+1]; i++;}
			else if (args[i].equals("--matchfile")) { matchfile = args[i+1]; i++;}
			else if (args[i].equals("--prefix")) { outputname = args[i+1]; i++;}
			else if (args[i].equals("--numiterations")) { numiters = Integer.parseInt(args[i+1]); i++;}
			else if (args[i].equals("--perr")) { perr = Double.parseDouble(args[i+1]); i++;}
			else if (args[i].equals("--savespace")) { SAVE_SPACE = true;}
			else if (args[i].equals("--enumerate")) { ENUMERATE = true; }
			else if (args[i].equals("--readclustersfirst")) { CLUSTERS_FIRST = true; }
			else 
				throw new IllegalArgumentException(args[i]+" is not a valid option.");
		}
		} catch (NumberFormatException e) {
			throw new IllegalArgumentException("--numiterations must be a positive integer, --perr must be a float between 0 and 1.");
		}
		
		if (clusterfile == null)
			throw new IllegalArgumentException("Cluster file is required.");
		if (assignmentfile == null)
			throw new IllegalArgumentException("Assignment file is required.");
		if (experimentfile == null)
			throw new IllegalArgumentException("Experiment file is required.");
		if (numiters == -1 && !ENUMERATE)
			throw new IllegalArgumentException("Number of iterations is required to sample.");
		if (perr == -1)
			throw new IllegalArgumentException("Mapping error probability (perr) is required.");
		if (numiters < 0 && !ENUMERATE) 
			throw new IllegalArgumentException("--numiterations must be a positive integer.");
		if (perr < 0 || perr > 1) 
			throw new IllegalArgumentException("--perr must be a value between 0 and 1.");
		
		System.out.println("\nclusterfile = " + clusterfile);
		System.out.println("assignmentfile = " + assignmentfile);
		System.out.println("experimentfile = " + experimentfile);
		if (matchfile != null)
			System.out.println("matchfile = " + matchfile);
		System.out.println("prefix = " + outputname);
		System.out.println("numiterations = " + numiters);
		System.out.println("perr = " + perr);
		System.out.println("save space? " + SAVE_SPACE);
		System.out.println("enumerate? " + ENUMERATE);
		System.out.println("clustersfirst? " + CLUSTERS_FIRST+"\n");
		
	}
	
	public void printUsage() {
		System.out.println("\n-------------------------------------------------------------------------------");
		System.out.println("MultiBreakSV.java: runs the MultiBreakSV method.");
		System.out.println("\nUsage: java -jar bin/MultiBreakSV.jar [ARGUMENTS] [OPTIONAL-ARGUMENTS]\n");
		System.out.println("  REQUIRED ARGUMENTS:\n" +
				"  --clusterfile STR\tGASV clusters file (run with maximal clusters)\n\n"+
				"  --assignmentfile STR\tFile of fragment assignments (with header).\n"+
				"\t\t\tColumn-delimited file where each line denotes a\n"+
				"\t\t\tsingle fragment alignment. Columns are\n"+
				"\t\t\t  1.FragmentID: Unique ID for each sequenced fragment\n"+
				"\t\t\t  2.DiscordantPairs: discordant pairs that indicate\n"+
				"\t\t\t    this alignment (these are clustered by GASV)\n"+
				"\t\t\t  3.NumErrors: number of errors in the alignment\n"+
				"\t\t\t  4.BasesInAlignment: the number of bases in the\n"+
				"\t\t\t    alignment (TargetEnd-TargetStart+1)\n"+
				"\t\t\t  5.ExperimentLabel: name of experiment/platform/run.\n"+
				"\t\t\tNote that the NumErrors and BasesInAlignment account\n"+
				"\t\t\tfor any concordantly-aligned pairs in the mapping\n"+
				"\t\t\t(for fragments in particular). BasesInAlignment is the\n"+
                "\t\t\tsum of the target lengths of the subalignments.\n\n"+
				"  --experimentfile STR\tFile of experiments/platforms/runs (with header).\n"+
				"\t\t\tColumn-delimited file where each line denotes a\n"+
				"\t\t\tsingle experiment.  Columns are\n"+
				"\t\t\t  1.ExperimentLabel: name of experiment/platform/run.\n"+
				"\t\t\t  2.Pseq: sequence error probability (e.g. 0.15 for\n"+
				"\t\t\t    PacBio, 0.01 for Illumina)\n"+
				"\t\t\t  3.LambdaD: expected number of fragments to support\n"+
				"\t\t\t    each SV (e.g. 3 for low-coverage PacBio)\n"+
				"\t\t\tThere are only multiple lines in this file if a hybrid\n"+
				"\t\t\texperiment is run.\n\n"+
				"  --numiterations INT\tNumber of iterations\n\n"+
				"  --perr FLOAT\t\tProbability that a fragment is unmapped\n"+
				"\t\t\t(e.g. 0.001)\n\n"+
				"  OPTIONAL ARGUMENTS:\n"+
				"  --prefix STR\t\tprefix to append to output files. Default is 'out'.\n\n"+
				"  --matchfile STR\tFile to initialize MCMC method to. The matchfile\n"+
				"\t\t\tis a single column of all ESP alignments to be set.\n"+
				"\t\t\tDefault is none.\n\n"+
				"  --savespace\t\tDoes not write mcmc.txt and .samplingprobs files;\n"+
				"\t\t\tthese tend to be large. Default is false.\n"+
				"  --enumerate\t\tEnumerates all possible combinations and explicitly\n" +
				"\t\t\tcomputes each state.  Recommended only for problems with < 100\n" +
				"\t\t\tstates.  Default is false.\n"+
				"  --readclustersfirst\t\tReads clusters file first rather than\n"+
				"\t\t\tassignmentsfile. Useful if jobs are split across clusters.\n"+
				"\t\t\tDefault is false.");
		System.out.println("-------------------------------------------------------------------------------\n");
	}
	
	public void readExperimentFile() throws FileNotFoundException {
		/** Read Experiment File **/
		Scanner scan = new Scanner(new File(experimentfile));
		scan.nextLine(); // skip header
		String explabel;
		Experiment exp;
		while(scan.hasNext()) {
			explabel = scan.next();
			exp = new Experiment(scan.nextDouble(),scan.nextDouble(),explabel);
			experiments.put(explabel, exp);
		}
		scan.close();
	}
	
	public HashSet<String> readAssignmentFile(HashSet<String> allesps) throws Exception {
		String explabel;
		/** Read Assignment File **/		
		Scanner scan = new Scanner(new File(assignmentfile));
		scan.nextLine(); // skip header
		String fragmentID;
		int numerrors;
		int length;
		String[] esps;
		while(scan.hasNext()) {		
			// process next line
			try {
				fragmentID = scan.next();
				esps = scan.next().split(",");
				numerrors = scan.nextInt();
				length = scan.nextInt();
				explabel = scan.next();
			} catch (Exception e) {
				System.out.println("Error parsing line. Check assignmentfile format.");
				throw e;
			}

			if (!experiments.containsKey(explabel)) {
				System.out.println("ERROR: Experiment Label " + explabel + " not in Experiments file.");
				System.exit(0);
			}
			
			if (CLUSTERS_FIRST) {
				// only add fragment if it's part of some cluster.
				boolean espInClust = false;
				for(int i=0;i<esps.length;i++) { 
					if (allesps.contains(esps[i])) {
						espInClust = true;
						break;
					}
				}
				
				// skip this fragment ID if it's not in a cluster.
				if (!espInClust) 
					continue;
			} 

			// If this is the first alignment for this fragment, create a new Fragment object.
			if (!fragments.containsKey(fragmentID))
				fragments.put(fragmentID,new Fragment(fragmentID));

			// Set fragment alignments.
			fa = new FragmentAlignment(experiments.get(explabel),esps,numerrors,length);
			fragments.get(fragmentID).addFragmentAlignment(fa);

			// populate allesps
			for(int i=0;i<esps.length;i++) {
				allesps.add(esps[i]); // will already be in allesps if CLUSTER_FIRST=true; doesn't hurt to re-add.
				esp2fragmentID.put(esps[i],fragmentID);
			}
		}
		scan.close();

		string_iter = fragments.keySet().iterator();
		while(string_iter.hasNext()) {
			fragmentID = string_iter.next();
			
			// initialize alignmentHits hash map.
			for(int i=0;i<fragments.get(fragmentID).fragmentalignments.size();i++) 
				alignmentHits.put(fragments.get(fragmentID).fragmentalignments.get(i),0);		
		}		
		return allesps;
	}
	
	public HashSet<String> readClustersFile(HashSet<String> allesps) throws FileNotFoundException {
		/** Read Clusters File **/
		/** Set All Breakpoints and MaxSupport Here **/
		/** Set clustermembership here too **/
		String[] esplist;
		String[] row;
		String cID;
		maxSupport = 0;
		HashMap<String,Boolean> clust = new HashMap<String,Boolean>(); // <cID,isCluster?>
		esp2cid = new HashMap<String,String>();
		Scanner scan = new Scanner(new File(clusterfile));
		boolean emitwarning = false;
		while(scan.hasNext()) {
			row = scan.nextLine().split("\t");

			if(row[0].equals("cID")) // skip header
				continue;

			cID = row[0];
			
			// skip maximal cluster
			if (cID.contains("."))
				continue;

			if(maxSupport < Integer.parseInt(row[1]))
				maxSupport = Integer.parseInt(row[1]);

			// Check to make sure there is at least one ESP in the ESPlist from above.
			esplist = row[4].split(",");
			boolean someESP = false;
			for(int i=0;i<esplist.length;i++) {
				esplist[i] = esplist[i].trim();
				
				// if there is an ESP in the clusterfile that is NOT in the ESPfiles, note this.
				if(!CLUSTERS_FIRST && !allesps.contains(esplist[i])) { 
					emitwarning = true;
					continue;
				}

				someESP = true; // at least one ESP is in the ESPlist.	
				// add esp to clustermembership.
				esp2cid.put(esplist[i],cID);	
				if (CLUSTERS_FIRST) // populate allesps
					allesps.add(esplist[i]);
			}

			if(!CLUSTERS_FIRST && !someESP) // Ignore if this cluster doesn't contain any ESPs in ESPlist.
				continue;
				
			// Add cluster to clust HashMap.  False = Conn Comp; True = Cluster
			if(row[2].equals("-1")) // conn. comp.
				clust.put(cID,false);
			else 
				clust.put(cID,true);
		}
		scan.close();

		if(emitwarning)
			System.out.println("WARNING! Some fragments in clusters are NOT in ESP file");
			

		/** Add Breakpoints for each Cluster **/
		string_iter = clust.keySet().iterator();
		int num = 1;
		while(string_iter.hasNext()) {
			cID = string_iter.next();
			breakpoints.put(cID,new MultiBreakSVBreakpoint(num,cID,clust.get(cID),maxSupport));
			num++;
		}		

		return allesps;
	}

	public void intersectClustersAndFragments() {
		System.out.println("Making sure that clusters and fragments agree with each other...");
		System.out.println("  Starting: " + breakpoints.size() + " breakpoints (clusters) and " + fragments.size() + " fragments.");
		/** Remove Fragments if they don't appear in any cluster **/
		string_iter = fragments.keySet().iterator();
		String fragmentname;
		HashSet<String> toRemove = new HashSet<String>();
		while(string_iter.hasNext()) {
			fragmentname = string_iter.next();
			ArrayList<FragmentAlignment> cur = fragments.get(fragmentname).fragmentalignments;
			ArrayList<FragmentAlignment> tmp = new ArrayList<FragmentAlignment>();
			for (int i=0;i<cur.size();i++) {
				int tmpnum = 0;
				for (int j=0;j<cur.get(i).alignments.size();j++) {						
					if (esp2cid.containsKey(cur.get(i).alignments.get(j))) { 
						System.out.println(cur.get(i).alignments.get(j) + " FOUND");
						tmpnum++;
					} else { 
						System.out.println(cur.get(i).alignments.get(j) + " MISSING");
					}
		
				}
				//System.out.println("  Assignment #"+i+": tmpnum = " + tmpnum + " and size = " + cur.get(i).alignments.size());
				if (tmpnum == cur.get(i).alignments.size())
					tmp.add(cur.get(i));
			}
			if(tmp.size()==0)
				toRemove.add(fragmentname);
			else {
				fragments.get(fragmentname).fragmentalignments = tmp;			
			}
		}

		if (toRemove.size() > 0) {
			System.out.println("  Removing " + toRemove.size() + " of " + fragments.size() +" fragments that are NOT in clusters file.");
			string_iter = toRemove.iterator();
			while(string_iter.hasNext()) {
				String fragment = string_iter.next();
				//System.out.println("  removing " + fragment);
				fragments.remove(fragment);
			}
		}
		 
		
		// Get Fragment IDs (arraylist), since we have trimmed fragments.
		getFragmentIDs();
				
		/** Remove clusters that don't contain any fragments **/
		string_iter = breakpoints.keySet().iterator();
		String bpname;
		String cid;
		toRemove = new HashSet<String>();
		while(string_iter.hasNext()) {
			bpname = string_iter.next();
			cid = breakpoints.get(bpname).cID;
			//System.out.println("checking " + cid);
			boolean keep = false;
			for(int i=0;i<fragmentIDs.size();i++) {
				ArrayList<String> tmp = fragments.get(fragmentIDs.get(i)).getAllAlignments();
				for(int j=0;j<tmp.size();j++) {
					if (esp2cid.get(tmp.get(j)).equals(cid))
						keep = true;
				}
			}
			if (!keep) {
				toRemove.add(bpname);
			}
		}
		if (toRemove.size() > 0) {
			System.out.println("Removing " + toRemove.size() + " breakpoints that do NOT contain a fragment in the assignmentfile.");
			string_iter = toRemove.iterator();
			while(string_iter.hasNext())
				breakpoints.remove(string_iter.next());
		}

	}
	public void getFragmentAlignments(boolean verbose) throws FileNotFoundException,Exception {
		/** Read Experiment File **/
		readExperimentFile();
		
		if (CLUSTERS_FIRST) {
			/** Read Clusters File **/
			HashSet<String> allesps = readClustersFile(new HashSet<String>());
			
			/** Read Assignment File **/	
			readAssignmentFile(allesps);	
		
		} else {
		
			/** Read Assignment File **/	
			HashSet<String> allesps = readAssignmentFile(new HashSet<String>());
				
			/** Read Clusters File **/
			readClustersFile(allesps);
		
		}
		
		// trim fragments and breakpoints data structures to 
		// only contain relevant fragments and breakpoints.
		intersectClustersAndFragments();
		
		System.out.println("\n==================\nFile Statistics:");
		System.out.println("  " + breakpoints.size() + " clusters (# of breakpoints).");
		string_iter = breakpoints.keySet().iterator();
		System.out.print("  ");
		while(string_iter.hasNext()) 
			System.out.print(" "+string_iter.next());
		System.out.println();
		
		// count ambiguous and space size.
		int ambig = 0;
		string_iter = fragments.keySet().iterator();
		String fragmentname;
		while(string_iter.hasNext()) {
			fragmentname = string_iter.next();
			// number of possible mappings is # of alignments plus one (no alignment).
			spacesize = spacesize.multiply(new BigInteger(""+(fragments.get(fragmentname).fragmentalignments.size()+1)));
			if(fragments.get(fragmentname).fragmentalignments.size() > 1) 
			ambig++;
		}
		System.out.println("  "+fragments.size() + " fragments added (" + ambig + " fragments have ambiguous alignments)");
		System.out.println("  The solution space has size " + spacesize);
		System.out.println("==================");
		//if(spacesize.equals(new BigInteger("1"))) {
		//	System.out.println("There is only 1 possible state. Exiting.");
		//	System.exit(0);
		//}
	}

	public void precomputePriors() {
		System.out.println("\nPrecomputing Priors");
			
		// compute all possible log probabilities of coverage.
		string_iter = experiments.keySet().iterator();
		String explabel;
		while(string_iter.hasNext()) {
			explabel = string_iter.next();
			experiments.get(explabel).logprobcov = new double[maxSupport+1];
			for(int i=1;i<maxSupport+1;i++)
				experiments.get(explabel).logprobcov[i] = logPoisson(i,experiments.get(explabel).lambdaD);	
		}

		// pre-compute prior.
		// note that this isn't quite the prior - this is a prior on the number of clust/conn comps.
		prior = new double[(int)Math.max(breakpoints.size()*2,100)];
		for(int i=0;i<prior.length;i++) {
			// gaussian prior with mean mu and standard deviation sigma
			//prior[i] = -Math.log(sigma*Math.sqrt(2*Math.PI))-Math.pow(i-mu,2)/(2*sigma*sigma);

			// exponential with lambda = log ( lambda * e ^ (-lambda * x) ) 
			//prior[i] = -explambda*i + Math.log(explambda);

			// uniform prior
			prior[i] = 0;
		}

		// Last but not least, compute the probability of delta_k ESPs supporting each variant.
		System.out.println("computing probability of delta_k supports...");
		deltas = new HashMap<String,Integer>();
		string_iter = diagrams.keySet().iterator();
		String cID;
		HashSet<String> seen = new HashSet<String>();
		String fragmentID;
		Iterator<String> leaf_iter;
		Iterator<String> exp_iter;
		while(string_iter.hasNext()) {
			cID = string_iter.next();
			seen.clear();
			leaf_iter = diagrams.get(cID).leaves.keySet().iterator();
			while(leaf_iter.hasNext()) {
				fragmentID = esp2fragmentID.get(diagrams.get(cID).leaves.get(leaf_iter.next()));
				seen.add(fragmentID);
			}
			deltas.put(cID,seen.size());
			exp_iter = experiments.keySet().iterator();
			while(exp_iter.hasNext()) {
				explabel = exp_iter.next();
				experiments.get(explabel).deltaprobs.put(cID,poisson(seen.size(),experiments.get(explabel).lambdaD));
			}
			//System.out.println(" cID has " + seen.size() + " fragments representing " + diagrams.get(cID).leaves.size() + " esps.");
		}
	}

	public void constructClusterDiagrams(boolean verbose) throws IOException {
		System.out.println("\nConstructing Cluster Diagrams from cluster file...");

		diagrams = new HashMap<String,MultiBreakSVDiagram>();

		Scanner scan = new Scanner(new File(clusterfile));
		String cID,local,names;
		Iterator<String> iter;

		int numclust = 0;
		int numconncomps = 0;
		String[] row = advance(scan);
		while(row != null) {
			
			// skip header
			if(row[0].equals("cID")) {
				row = advance(scan);
				continue;
			}

			cID = row[0];
			local = row[2];
			names = row[4]; 

			// if it's not in breakpoints, skip.
			if (!breakpoints.containsKey(cID)) {
				row = advance(scan);
				continue;
			}

			// create diagram
			MultiBreakSVDiagram diagram = new MultiBreakSVDiagram();
			if (!local.equals("-1")) { // cluster
				if(verbose)
					System.out.println("constructing "+cID+" from cluster");
				diagram.constructFromCluster(cID,names.split(", "));
				numclust++;
				
				row = advance(scan);
				
			} else { // connected component
				if(verbose)
					System.out.println("constructing "+cID+" from connected component");
				ArrayList<String[]> conncomps = new ArrayList<String[]>();
				ArrayList<String> cIDs = new ArrayList<String>();
				
				row = advance(scan);				
				while(row != null && row[0].startsWith(cID)) {
					conncomps.add(row[4].split(", "));
					cIDs.add(row[0]);
					row = advance(scan);
				}
		
				diagram.constructFromConnComp(cID,conncomps,cIDs,true,verbose);	
				numconncomps++;
			}

			// add diagram to candidates.

			if(verbose)
				System.out.println("adding diagram to candidates");

			// update hash tables
			iter = diagram.leaves.keySet().iterator();
			String esp;
			while(iter.hasNext()) {
				esp = iter.next();
				if(esp2cid.containsKey(esp))
					esp2cid.put(esp,cID);
			}

			// add diagram file name to cluster array
			if(verbose)
				System.out.println("adding diagram to cluster array");
			diagrams.put(cID,diagram);
			
		}
		scan.close();

		System.out.println("Read in " + numclust + " clusters and " + numconncomps + " connected components from the GASV file.");
		System.out.println("This leaves "+diagrams.size()+" clusters/conncomps ");
	}
	
	public String[] advance(Scanner scan) {
		if(scan.hasNext())
			return scan.nextLine().split("\t");
		else 
			return null;
	}

	public void fillNonSingletonArray() {
		bpsWithNonSingletonSupport = new ArrayList<String>();
		string_iter = breakpoints.keySet().iterator();
		Iterator<String> leaf_iter;
		String cID;
		String fragmentID;
		int uniqueCounts;
		while(string_iter.hasNext()) {
			uniqueCounts = 0;
			cID = string_iter.next();
			leaf_iter = diagrams.get(cID).leaves.keySet().iterator();
			while(leaf_iter.hasNext()) {
				fragmentID = esp2fragmentID.get(diagrams.get(cID).leaves.get(leaf_iter.next()));
				if(fragments.containsKey(fragmentID) && fragments.get(fragmentID).fragmentalignments.size() == 1)
					uniqueCounts++;
			}
			
			// if there are 0 unique counts, then this move won't change anything.
			// if there is only 1 unique count, it can be reached via a naive move.
			// Ignore if this is the case.
			if(uniqueCounts > 1) 
				bpsWithNonSingletonSupport.add(cID);
		}
	}
	
	public MultiBreakSVAssignmentMatrix runMCMC(String file) throws IOException, ClassNotFoundException {
		System.out.println("\nRunning MCMC for " + numiters + " iterations");
		System.out.println("  Short Burnin = " + burnin);
		System.out.println("  Long Burnin = " + (numiters-burnin));
		System.out.println("  Thinning every " + skip + " iterations");	

		BufferedWriter writer = null;
		if (file != null) { 
			writer = new BufferedWriter(new FileWriter(file));
			writer.write("Iteration\tLogProb\tNumBreakpoints\tMadeMove?\n");
		}
		
		// initialize
		System.out.println("=================\nInitializing Matrix...");
		MultiBreakSVAssignmentMatrix A = initializeMatrix();

		System.out.println("=================\nInitialized Log Probability is "+A.logprob+"\n");
		A.printSupport();

		// ITERATIONS
		int divide = numiters/25;
		MultiBreakSVAssignmentMatrix prevA;
		String fragmentID;
		for(int iter=0;iter<numiters;iter++) {
			if(iter%divide==0)
				System.out.println("Iteration " + iter + " of " + numiters+": prob is " + A.logprob + ", "+counts.size() + " states have been sampled, "+logprobs.size()+" state log probs have been computed");

			// iterate
			prevA = A;
			A = iterate(A,false,iter);
			if(writer != null) 
				writer.write(iter+"\t"+A.logprob);

			// Increment Counters 
			incrementAllCounters(A,iter);
			
			// write number of breakpoints
			if(iter % skip == 0) {// only write every 'skip' times. 
				string_iter = A.supports.keySet().iterator();
				int numbps = 0;
				while(string_iter.hasNext()) {
					numbps+=A.supports.get(string_iter.next()).size();
				}
				if(writer != null)
					writer.write("\t"+numbps);
			}
			else 
				if(writer != null)
					writer.write("\t0");
			
			if(writer != null) {
			if(A.sameAs(prevA))
				writer.write("\t0\n");
			else
				writer.write("\t1\n");
			}

			// if iter is divisible by 1 million AND iter is bigger than burnin, write intermediate file.
			if(iter%intermediateJump == 0 && iter > 0 && iter > burnin+100000) {
				writeIntermediateFiles(iter);
			}
			//if(iter == 0) {
			//	writeFirstFiles();
			//}
		}

		System.out.println("  After "+numiters+" iterations, probability is "+A.logprob+"\n");
		A.printSupport();

		System.out.println(moveCounts[0]+" of " + proposalCounts[0] + "("+(moveCounts[0]/(double)proposalCounts[0])+") Naive Moves");
		System.out.println(moveCounts[1]+" of " + proposalCounts[1] + "("+(moveCounts[1]/(double)proposalCounts[1])+") AlterSV Moves");
		
		// get average assignment counts
		string_iter = fragments.keySet().iterator();
		while(string_iter.hasNext()) {
			fragmentID = string_iter.next(); 

			fragments.get(fragmentID).avgNoassignmentshort/=(numiters-burnin);
			for(int j=0;j<fragments.get(fragmentID).avgAssignmentsshort.size();j++) {
				double count = fragments.get(fragmentID).avgAssignmentsshort.get(j);
				fragments.get(fragmentID).avgAssignmentsshort.set(j,count/(numiters-burnin));
			}

			fragments.get(fragmentID).avgNoassignmentlong/=burnin;
			for(int j=0;j<fragments.get(fragmentID).avgAssignmentslong.size();j++) {
				double count = fragments.get(fragmentID).avgAssignmentslong.get(j);
				fragments.get(fragmentID).avgAssignmentslong.set(j,count/burnin);
			}

			fragments.get(fragmentID).avgNoassignmentthinned/=(numiters/skip);
			for(int j=0;j<fragments.get(fragmentID).avgAssignmentsthinned.size();j++) {
				double count = fragments.get(fragmentID).avgAssignmentsthinned.get(j);
				fragments.get(fragmentID).avgAssignmentsthinned.set(j,count/(numiters/skip));
			}

			// don't need to average first A.
			// don't need to average last A. 
		}

		if(writer != null)
			writer.close();
		
		return A;

	}
	
	public void incrementAllCounters(MultiBreakSVAssignmentMatrix A,int iter) {
		// increment counts
		String str = A.summaryString();
		if(!counts.containsKey(str)) {
			counts.put(str,1);
			logprobs.put(str,A.logprob);
		} else {
			if(logprobs.get(str) != A.logprob)
				throw new Error("Error! Stored prob is different from this prob.");
			counts.put(str,counts.get(str)+1);
		}
		
		// increment breakpoint supports
		string_iter = breakpoints.keySet().iterator();
		String cID;
		while(string_iter.hasNext()) {
			cID = string_iter.next();
			if(A.supports.get(cID).size() == 0) { // no support. increment 0.
				if(iter == 0) 
					breakpoints.get(cID).supportsfirst[0]++;

				if(iter > burnin)
					breakpoints.get(cID).supportsshort[0]++;

				if(iter > numiters-burnin)
					breakpoints.get(cID).supportslong[0]++;

				if(iter % skip == 0) 
					breakpoints.get(cID).supportsthinned[0]++;

				if(iter == numiters-1) 
					breakpoints.get(cID).supportslast[0]++;
				continue;
			}

			// otherwise, increment appropriate supports.
			Iterator<String> exp_iter = experiments.keySet().iterator();
			Experiment exp;
			while(exp_iter.hasNext()) {
				exp = experiments.get(exp_iter.next());
				for(int j=0;j<A.supports.get(cID).get(exp).size();j++) {
					if(iter == 0) 
						breakpoints.get(cID).supportsfirst[A.supports.get(cID).get(exp).get(j)]++;

					if(iter > burnin)
						breakpoints.get(cID).supportsshort[A.supports.get(cID).get(exp).get(j)]++;

					if(iter > numiters-burnin)
						breakpoints.get(cID).supportslong[A.supports.get(cID).get(exp).get(j)]++;

					if(iter % skip == 0) 
						breakpoints.get(cID).supportsthinned[A.supports.get(cID).get(exp).get(j)]++;

					if(iter == numiters-1) 
						breakpoints.get(cID).supportslast[A.supports.get(cID).get(exp).get(j)]++;
				}
			}
		}

		// increment fragment assignment counts
		string_iter = fragments.keySet().iterator();
		String fragmentID;
		while(string_iter.hasNext()) {
			fragmentID = string_iter.next(); 
			int ind = A.assignments.get(fragmentID);

			if(iter == 0) {
				if(ind == -1) 
					fragments.get(fragmentID).avgNoassignmentfirst++;
				else {
					fragments.get(fragmentID).avgAssignmentfirst.set(ind,1.0);
				}
			}

			if(iter > burnin) {
				if(ind == -1) 
					fragments.get(fragmentID).avgNoassignmentshort++;
				else {
					double newcount = fragments.get(fragmentID).avgAssignmentsshort.get(ind)+1;
					fragments.get(fragmentID).avgAssignmentsshort.set(ind,newcount);
				}
			}

			if(iter > numiters-burnin) {
				if(ind == -1) 
					fragments.get(fragmentID).avgNoassignmentlong++;
				else {
					double newcount = fragments.get(fragmentID).avgAssignmentslong.get(ind)+1;
					fragments.get(fragmentID).avgAssignmentslong.set(ind,newcount);
				}
			}

			if(iter % skip == 0) {
				if(ind == -1) 
					fragments.get(fragmentID).avgNoassignmentthinned++;
				else {
					double newcount = fragments.get(fragmentID).avgAssignmentsthinned.get(ind)+1;
					fragments.get(fragmentID).avgAssignmentsthinned.set(ind,newcount);
				}
			}

			if(iter == numiters-1) {
				if(ind == -1) 
					fragments.get(fragmentID).avgNoassignmentlast++;
				else {
					fragments.get(fragmentID).avgAssignmentlast.set(ind,1.0);
				}
			}
		}
	}

	public MultiBreakSVAssignmentMatrix initializeMatrix() throws IOException, ClassNotFoundException {
		Fragment f;
		FragmentAlignment fa;
		String fragmentID;
		MultiBreakSVAssignmentMatrix A = new MultiBreakSVAssignmentMatrix();

		if (matchfile == null) { // no matchfile = initialize uniq alignments only.

			System.out.println("  Initialize A by using random alignments.");
			string_iter = fragments.keySet().iterator();
			
			int num = 0;
			int skipped = 0;
			while(string_iter.hasNext()) {
				fragmentID = string_iter.next(); 
				f = fragments.get(fragmentID);

				// SELECT EACH ONE UNIFORMLY
				if(rand.nextDouble() < perr) {
					skipped++;
					A.assignments.put(fragmentID,-1);
				} else {
					num++;
					A.assignments.put(fragmentID,rand.nextInt(f.fragmentalignments.size()));
				}
			}
			/*
			// Code to initialize to All errors
			System.out.println("initialize using ALL errors");
			int_iter = fragments.keySet().iterator();
			String fragmentID;
			while(int_iter.hasNext()) {
				fragmentID = int_iter.next(); 
				s = fragments.get(fragmentID);

				A.assignments.put(fragmentID,-1);
			}
			 */
			System.out.println("  "+num+" esps recorded, " + skipped + " esps skipped.");

		} else { // matchfile: set as initial assignment
			// matchfile is a single column of all ESP alignments to be set.
			System.out.println("  Set file " + matchfile + " as initial A.");
			Scanner scan = new Scanner(new File(matchfile));
			ArrayList<String> espstoset = new ArrayList<String>();
			while(scan.hasNext()) {
				String l = scan.next();
				espstoset.add(l);
			}
			scan.close();

			string_iter = fragments.keySet().iterator();
			int numset = 0;
			int numbad = 0;
			int numerr = 0;
			boolean foundall;
			boolean foundany;
			int tot = 0;
			while(string_iter.hasNext()) {
				tot++;
				fragmentID = string_iter.next(); 
				f = fragments.get(fragmentID);

				foundall = false;
				foundany = false;
				for(int j=0;j<f.fragmentalignments.size();j++) {
					fa = f.fragmentalignments.get(j);	
					int thisval = 0;
					for(int l=0;l<espstoset.size();l++) {
						if (fa.contains(espstoset.get(l))) {
							thisval++;
							foundany = true;
							break;
						}
					}
					if(thisval == fa.alignments.size()) {
						A.assignments.put(fragmentID,j);
						numset++;
						foundall = true;
						break;
					} 
				}
				if (!foundany && !foundall) 
					numerr++;
				
				if(foundany && !foundall)  {
					System.out.println("WARNING: fragment " + fragmentID + " was not set because there wasn't an alignment were all ESPs were in matchfile.");
					numbad++;
				}
		}

			System.out.println(espstoset.size()+" esps recorded.");
			System.out.println(tot+ " total fragments");
			System.out.println(numset +" fragments set in A, and " + numerr + " fragments set to eror.");
			System.out.println(numbad + " fragments had an ESP in the matchfile but no assignment was set.");
		}
		System.out.println("  Setting Breakpoint According to Diagrams...");
		A.setBreakpoints(diagrams);
		
		A.printSupport();
		A.printAssignments();
		
		System.out.println("  Computing Log Prob...");
		A.computeLogProb(true);
		String str = A.summaryString();
		logprobs.put(str,A.logprob);
		System.out.println("  Done initializing A");
		
		return A;
	}

	public MultiBreakSVAssignmentMatrix iterate(MultiBreakSVAssignmentMatrix A,boolean verbose,int iter) throws IOException, ClassNotFoundException {		
		Alogq = Double.NaN;
		Anewlogq = Double.NaN;
		MultiBreakSVAssignmentMatrix Anew;
		double ratio;
		double r = rand.nextDouble();

		// IMPLEMENT LAZY CHAIN: W/ prob 1/2, stay.
		if(r < 0.5) {
			if(!silent)
				System.out.println("#"+iter+": LAZY CHAIN - returning A.");
			return A;
		} 
		// draw new random variable
		r = rand.nextDouble();

		// Make Naive Move 90% of the time (100% of the time if there are no SVs to pick from )
		if (r < 0.9 || bpsWithNonSingletonSupport.size() == 0) {
		moveType = MoveType.Naive;
		Anew = naiveMove(A,verbose);
		proposalCounts[0]++;
		if (A.sameAs(Anew)) 
			throw new Error("Error! Now A != Anew for move Naive");
		
		} else { // Otherwise, make SVmove 10% of the time
			moveType = MoveType.AlterSV;
			Anew = svMove(A,verbose);
			proposalCounts[1]++;
			if (A.sameAs(Anew)) 
				throw new Error("Error! Now A != Anew for move AlterSV");
		}
		
		String str = Anew.summaryString();
		if(logprobs.containsKey(str))
			Anew.logprob = logprobs.get(str);
		else
			Anew.computeLogProb(false);	
		
	
		

		if(Double.isNaN(Alogq) || Double.isNaN(Anewlogq) ||
				Double.isInfinite(Alogq) || Double.isInfinite(Anewlogq))
			throw new Error("Alogq = " + Alogq + " and Anewlogq = " + Anewlogq + " after move type "+ moveType);

		ratio = Math.min(1,Math.exp((Anew.logprob+Alogq)-(A.logprob+Anewlogq)));

		r = rand.nextDouble();
		if(ratio > r) { // move
			if(!silent) {
				System.out.println("#"+iter+": MAKING MOVE OF TYPE " + moveType);
				System.out.println("\tFragment "+ naivefragment+": Moving from " + prevassign + " to " + thisassign);
				System.out.println("\told log prob = " + A.logprob);
				System.out.println("\tnew log prob = " + Anew.logprob);
				System.out.println("\tAlogq = " + Alogq);
				System.out.println("\tAnewlogq = " + Anewlogq);
				System.out.println("\tratio = " + ratio);
				System.out.println("\trandom() = " + r);
				System.out.println();
			}
			if(moveType == MoveType.Naive)
			moveCounts[0]++;
			else
				moveCounts[1]++;
			
			return Anew;

		} else {
			if(!silent) {
				System.out.println("#"+iter+": STAYING PUT FOR TYPE " + moveType);
				System.out.println("\told log prob = " + A.logprob);
				System.out.println("\tnew log prob = " + Anew.logprob);
				System.out.println("\tAlogq = " + Alogq);
				System.out.println("\tAnewlogq = " + Anewlogq);
				System.out.println("\tratio = " + ratio);
				System.out.println("\trandom() = " + r);
				System.out.println();
			}

			return A;
		}
	}
	
	public MultiBreakSVAssignmentMatrix naiveMove(MultiBreakSVAssignmentMatrix A, boolean verbose) throws IOException, ClassNotFoundException {
		MultiBreakSVAssignmentMatrix Anew = new MultiBreakSVAssignmentMatrix(A);

		// STEP 1: SAMPLE STROBE
		String fragmentID = sampleFragmentUniformly();
		if (verbose)
			System.out.println("CHOOSING FRAGMENT "+fragmentID + " : " + fragments.get(fragmentID).fragmentalignments.size() + " fragment alignments.");
		naivefragment = fragmentID;
		ArrayList<FragmentAlignment> sa = fragments.get(fragmentID).fragmentalignments;
		prevassign = A.assignments.get(fragmentID);	
				
		if(verbose)
			System.out.println("--> previous assignment is " + prevassign);
		
		double r = rand.nextDouble();		
		
		if(prevassign == -1) { // if previous assignment was an error
			
			thisassign = 0;
				
			// get index w/ cumulative probability > r.
			double sum = (double)1/sa.size();
			while(thisassign<sa.size()-1 && r > sum) {
				thisassign++;
				sum+=(double)1/sa.size();
			}
			
			// set assignment
			Anew.assignments.put(fragmentID,thisassign);
			
			if(verbose)
				System.out.println("--> new assignment is " + thisassign);
			
			// set Alogq
			if(sa.size() == 1) { // only one alignment - prob is 1
				Alogq = -Math.log(fragments.size()) + Math.log(1);
				Anewlogq = -Math.log(fragments.size()) + Math.log(1);
			} else {
				Alogq = - Math.log(fragments.size()) + Math.log(perr);
				Anewlogq = - Math.log(fragments.size()) - Math.log(sa.size());
			}

		} else { // if previous assignment was not an error

			if(r < perr || sa.size() == 1) { 
				
				// set assignment to an error
				thisassign = -1;
				Anew.assignments.put(fragmentID,thisassign);
			
				if(verbose)
					System.out.println("--> new assignment is " + thisassign);
				
				// set Anewlogq
				if(sa.size() == 1) { // if there is only one alignment, then this prob is 1
					Alogq = -Math.log(fragments.size()) + Math.log(1);
					Anewlogq = -Math.log(fragments.size()) + Math.log(1);
				} else {
					Alogq = -Math.log(fragments.size()) - Math.log(sa.size());
					Anewlogq = -Math.log(fragments.size()) + Math.log(perr);
				}
			} else { // Otherwise, set to an alignment. Normalize by (1-prevassignment prob).
				
				thisassign = 0;
				r = rand.nextDouble();

				// get index (that is not the previous index) w/ cumulative probability > r.
				double sum = (double)1/(sa.size()-1);
				while(thisassign<sa.size()-1 && (thisassign == prevassign || r > sum)){
					thisassign++;
					if(thisassign != prevassign) // skip prev assignment.
						sum+=(double)1/(sa.size()-1);
				}
				
				if(thisassign >= sa.size()) {
					System.out.println("ERROR! alignment not selected for fragment " + fragmentID);
					throw new Error();
				}
				
				if(verbose)
					System.out.println("--> new assignment is " + thisassign);
				
				// set assignment 
				Anew.assignments.put(fragmentID,thisassign);	
				
				// set Anewlogq & Alog
				Anewlogq = -Math.log(fragments.size()) + Math.log(1-perr) - Math.log(sa.size()-1);
				Alogq = -Math.log(fragments.size()) + Math.log(1-perr) - Math.log(sa.size()-1);
			}
		}

		if(!silent)
			System.out.println(fragmentID+" alignment moved from " + A.assignments.get(fragmentID) + " to " + Anew.assignments.get(fragmentID));

		// only need to recompute breakpoints for connected components that have changed.
		ArrayList<String> changedFragmentIDs = new ArrayList<String>();
		changedFragmentIDs.add(fragmentID);
		HashMap<MultiBreakSVDiagram,Boolean> changedClusters = getChangedClusters(A,Anew,changedFragmentIDs);

		// reset breakpoints for Anew.
		Iterator<MultiBreakSVDiagram> iter2 = changedClusters.keySet().iterator();
		while(iter2.hasNext()) 
			Anew.setBreakpoints(iter2.next());
		
		return Anew;
	}

	public MultiBreakSVAssignmentMatrix svMove(MultiBreakSVAssignmentMatrix A, boolean verbose) throws IOException, ClassNotFoundException {
		MultiBreakSVAssignmentMatrix Anew = new MultiBreakSVAssignmentMatrix(A);

		// STEP 1: SAMPLE SV with TOTAL SUPPORT > 1
		String cID = sampleBreakpointUniformly();
		if (verbose)
			System.out.println("CHOOSING CLUSTER "+cID + " : " + deltas.get(cID) + " support (unique + ambig alignments).");
		
		// only need to recompute breakpoints for connected components that have changed.
		ArrayList<String> inds = new ArrayList<String>();
		
		
		Iterator<String> leaf_iter = diagrams.get(cID).leaves.keySet().iterator();
		String fragmentID;
		int numchanged = 0;
		while(leaf_iter.hasNext()) {
			
			fragmentID = esp2fragmentID.get(diagrams.get(cID).leaves.get(leaf_iter.next()));
			if(fragments.containsKey(fragmentID) && 
					fragments.get(fragmentID).fragmentalignments.size() == 1) {
				numchanged++;
				prevassign = A.assignments.get(fragmentID);	
			
				
				if(prevassign == -1) { // if previous assignment was an error
					thisassign = 0;
				} else { // if previous assignment was not an error
					thisassign = -1;
				}

				Anew.assignments.put(fragmentID,thisassign);
				
				if(verbose) {
					System.out.println("Moving fragment "+fragmentID+": original assignment is " + prevassign + " and new assignment is " + thisassign);
				}
				
				// add fragment to set of changed indices.
				inds.add(fragmentID);
			}
		}
		
		if(!silent) { 
			System.out.println("changing " + numchanged + " unique alignments for " + cID);
			string_iter = A.assignments.keySet().iterator();
			while(string_iter.hasNext()) {
				fragmentID = string_iter.next();
				System.out.println("old="+A.assignments.get(fragmentID) + " and new=" + Anew.assignments.get(fragmentID));
			}
		}

		// set Alogq
		Alogq = -Math.log(bpsWithNonSingletonSupport.size());
		Anewlogq = -Math.log(bpsWithNonSingletonSupport.size());
		
		HashMap<MultiBreakSVDiagram,Boolean> changedClusters = getChangedClusters(A,Anew,inds);

		// reset breakpoints for Anew.
		Iterator<MultiBreakSVDiagram> iter2 = changedClusters.keySet().iterator();
		while(iter2.hasNext()) 
			Anew.setBreakpoints(iter2.next());
		
		return Anew;
	}

	public ArrayList<Double> normalize(ArrayList<Double> list) {

		double sum = 0;
		for(int i=0;i<list.size();i++) 
			sum+=list.get(i);

		for(int i=0;i<list.size();i++)
			list.set(i,list.get(i)/sum);

		return list;
	}


	public ArrayList<String> getUnion(String cID1,String cID2) {
		MultiBreakSVDiagram diagram1 = diagrams.get(cID1);
		MultiBreakSVDiagram diagram2 = diagrams.get(cID2);

		Iterator<String> iter1 = diagram1.leaves.keySet().iterator();
		ArrayList<String> list1 = new ArrayList<String>();
		ArrayList<String> strlist1 = new ArrayList<String>();
		while (iter1.hasNext())  {
			strlist1.add(iter1.next());
			list1.add(getFragmentID(strlist1.get(strlist1.size()-1)));
		}

		Iterator<String> iter2 = diagram2.leaves.keySet().iterator();
		ArrayList<String> list2 = new ArrayList<String>();
		ArrayList<String> strlist2 = new ArrayList<String>();
		while (iter2.hasNext())  {
			strlist2.add(iter2.next());
			list2.add(getFragmentID(strlist2.get(strlist2.size()-1)));
		}

		ArrayList<String> union = new ArrayList<String>();
		for(int i=0;i<list1.size();i++) {
			if(list2.contains(list1.get(i))) {
				int j = list2.indexOf(list1.get(i));
				union.add(strlist1.get(i)+" "+strlist2.get(j));
			}

		}
		return union;
	}

	public String sampleFragmentUniformly() {
		int j= (int)(rand.nextDouble()*fragments.size());
		return (String) fragments.keySet().toArray()[j];
	}

	// NOTE: only samples bps with >1 unique alignments.
	public String sampleBreakpointUniformly() {
		int j = (int)(rand.nextDouble()*bpsWithNonSingletonSupport.size());
		return bpsWithNonSingletonSupport.get(j);
	}
	
	/** Enumerates all states rather than sampling. **/
	public void enumerate() throws IOException, ClassNotFoundException {
		// get max num align
		Iterator<Fragment> iter = fragments.values().iterator();
		Fragment s;
		MAX_NUM_ALIGN = 0;
		while(iter.hasNext()) {
			s = iter.next();
			if(s.fragmentalignments.size() > MAX_NUM_ALIGN)
				MAX_NUM_ALIGN = s.fragmentalignments.size();
		}
		System.out.println("Maximum number of alignments are " + MAX_NUM_ALIGN);
		
		// get spacesize.
		spacesize = new BigInteger(""+(MAX_NUM_ALIGN+1)).pow(fragments.size());
		System.out.println("State space size is " + spacesize);
		
		// now, enumerate.
		MultiBreakSVAssignmentMatrix A = new MultiBreakSVAssignmentMatrix();

		string_iter = fragments.keySet().iterator();
		String fragmentID;
		while(string_iter.hasNext()) {
			fragmentID = string_iter.next(); 
			A.assignments.put(fragmentID,-1);
		}
		
		char[] assignment = new char[fragments.size()];
		boolean ok = true;
		int n;
		long i=0;
		ArrayList<Double> logposterior = new ArrayList<Double>();
		ArrayList<String> assignmentstrings = new ArrayList<String>();
		ArrayList<int[]> assignments = new ArrayList<int[]>();
		while(i < spacesize.longValue()) {
			if(i % 100 == 0) { 
				System.out.println(i + " of " + spacesize.floatValue());
			}
			
			// get quaternary representation
			assignment = quaternaryToAssignment(new BigInteger(""+i));
			//System.out.println("i="+i + ": " + new String(assignment));
			
			// set alignment
			ok = true;
			for(int j=0;j<assignment.length;j++) {
				n = Character.getNumericValue(assignment[j]);
				//if(n == fragments.get(fragmentIDs.get(j)).fragmentalignments.size())
				if(n == MAX_NUM_ALIGN)
					A.assignments.put(fragmentIDs.get(j),-1);
				else {
					if(n >= fragments.get(fragmentIDs.get(j)).fragmentalignments.size()) {
						ok = false;
					}
					A.assignments.put(fragmentIDs.get(j),n);
				}
			}
			
			if(!ok) {
				// advance i 
				long oldi = i;
				int j= 0;
				n = Character.getNumericValue(assignment[j]);
				while(n == MAX_NUM_ALIGN || n < fragments.get(fragmentIDs.get(j)).fragmentalignments.size()) {
					j++;
					if(j == fragmentIDs.size()-1)
						break;
					n = Character.getNumericValue(assignment[j]);
				}
				
				assignment[j] = Character.forDigit(MAX_NUM_ALIGN,BASE);
				j++;
				while(j < fragmentIDs.size()) {
					assignment[j] = Character.forDigit(0,BASE);
					j++;
				}
				i = assignmentToQuaternary(assignment).longValue();
				if (i == oldi)
					i++;
				
				if(i < 0) {
					System.out.println("NEG: " + oldi + " -> " + i + ": " + new String(assignment));
					System.exit(0);
				}
				continue;
			}
			
			//System.out.println(new String(assignment) + " --> " + concatAlignments(A));
			A.setBreakpoints(diagrams);
			A.computeLogProb(false);
			logposterior.add(A.logprob);
			assignmentstrings.add(concatAlignments(A));
			int[] tmp = new int[fragmentIDs.size()];
			for (int j=0;j<fragmentIDs.size();j++) {
				tmp[j] = A.assignments.get(fragmentIDs.get(j));
			}
			assignments.add(tmp);				
			i++;
		}
		
		double normalization = 0;
		for (int j=0;j<logposterior.size();j++) {
			normalization += Math.exp(logposterior.get(j));
		}
		System.out.println("Normalization is "+ normalization);
		ArrayList<Double> normposterior = new ArrayList<Double>();
		
		String statefile = outputname+".allstates";
		BufferedWriter writer = new BufferedWriter(new FileWriter(statefile));
		writer.write("#\tLogPosterior\tPosterior\tNormalizedPosterior");
		for(int j=0;j<fragmentIDs.size();j++) {
			writer.write("\t"+fragmentIDs.get(j));
		}
		writer.write("\n");
		for (int j=0;j<logposterior.size();j++) {
			normposterior.add(Math.exp(logposterior.get(j))/normalization);
			writer.write(j+"\t"+logposterior.get(j)+"\t"+Math.exp(logposterior.get(j))+"\t"+Math.exp(logposterior.get(j))/normalization+"\t"+assignmentstrings.get(j)+"\n");
		}
		writer.close();
		
		System.out.println("Results written to " + statefile);
		
		String assignfile = outputname+".assignment.posteriors";
		writer = new BufferedWriter(new FileWriter(assignfile));
		writer.write("FragmentID\tAssignmentIndex\tNormPosterior\tDiscordantPairs\n");
		for(int k=0;k<fragmentIDs.size();k++) {
			fragmentID = fragmentIDs.get(k);
			writer.write(fragments.get(fragmentID).fragmentID+"\t-1\t"+
					sumStates(assignments,normposterior,k,-1)+"\terror\n");
			for(int j=0;j<fragments.get(fragmentID).fragmentalignments.size();j++) {
				writer.write(fragments.get(fragmentID).fragmentID+"\t"+
						j+"\t"+
						sumStates(assignments,normposterior,k,j)+"\t"+
						fragments.get(fragmentID).concatenateLabels(j)+"\n");
			}
		}
		writer.close();
		System.out.println("File for posterior assignents is " + assignfile);
	}
	
	public double sumStates(ArrayList<int[]> assignments, 
			ArrayList<Double> normposterior, 
			int fragmentIndex, int assignmentIndex) {
		double sum = 0;

		for (int i=0;i<assignments.size();i++) {
			if (assignments.get(i)[fragmentIndex] == assignmentIndex)
				sum+= normposterior.get(i);
		}
		return sum;		
	}
	
	public char[] quaternaryToAssignment(BigInteger b) {
		char[] assignment = new char[fragments.size()];
		char[] tmp;
		int shift;

		Arrays.fill(assignment,'0');
		tmp = b.toString(MAX_NUM_ALIGN+1).toCharArray(); // convert to base MAX_NUM_ALIGN+1 (for error)s
		shift = assignment.length-1;
		for(int j=tmp.length-1;j>=0;j--) {
			assignment[shift] = tmp[j];
			shift--;
		}
		return assignment;
	}
	
	public BigInteger assignmentToQuaternary(MultiBreakSVAssignmentMatrix A) {
		char[] assignment = new char[fragmentIDs.size()];
		for(int i=0;i<fragmentIDs.size();i++) {
			if(A.assignments.get(fragmentIDs.get(i)) == -1)
				assignment[i] = Character.forDigit(MAX_NUM_ALIGN,BASE);
			else
				assignment[i] = Character.forDigit(A.assignments.get(fragmentIDs.get(i)),BASE);	
		}
		
		return new BigInteger(new String(assignment),BASE);
	}
	
	public BigInteger assignmentToQuaternary(char[] assignment) {
		return new BigInteger(new String(assignment),MAX_NUM_ALIGN+1);
	}

	public String concatAlignments(MultiBreakSVAssignmentMatrix A) {
		String str = ""+A.assignments.get(fragmentIDs.get(0));
		for(int i=1;i<fragmentIDs.size();i++) {
			str+="\t"+A.assignments.get(fragmentIDs.get(i));
		}
		return str;
	}

	public HashMap<MultiBreakSVDiagram,Boolean> getChangedClusters(MultiBreakSVAssignmentMatrix A, 
			MultiBreakSVAssignmentMatrix Anew, ArrayList<String> inds) {
		HashMap<MultiBreakSVDiagram,Boolean> changedClusters = new HashMap<MultiBreakSVDiagram,Boolean>();
		String cID;
		String fragmentID;
		FragmentAlignment fa;
		for(int k=0;k<inds.size();k++) {
			fragmentID = inds.get(k);

			if (A.assignments.get(fragmentID) == Anew.assignments.get(fragmentID))
				continue;

			if(A.assignments.get(fragmentID) != -1) {
				fa = fragments.get(fragmentID).fragmentalignments.get(A.assignments.get(fragmentID));
				for(int i=0;i<fa.alignments.size();i++) {
					if(!esp2cid.containsKey(fa.alignments.get(i))) {
						System.out.println("KEY: "+fa.alignments.get(i)+" NOT FOUND for fragment ID " + fragmentID + " and index = " + i);
						throw new Error();
					}
					cID = esp2cid.get(fa.alignments.get(i));
					if(!cID.equals(""))
						changedClusters.put(diagrams.get(cID),true);
				}
			}

			if(Anew.assignments.get(fragmentID) != -1) {
				fa = fragments.get(fragmentID).fragmentalignments.get(Anew.assignments.get(fragmentID));
				for(int i=0;i<fa.alignments.size();i++) {
					if(!esp2cid.containsKey(fa.alignments.get(i))) {
						System.out.println("KEY: "+fa.alignments.get(i)+" NOT FOUND for fragment ID " + fragmentID + " and index = " + i);
						throw new Error();
					}
					cID = esp2cid.get(fa.alignments.get(i));
					if(!cID.equals(""))
						changedClusters.put(diagrams.get(cID),true);
				}
			}

		}
		return changedClusters;
	}


	public double logPoisson(int k,double lambdaD) {
		double val = -lambdaD + k*Math.log(lambdaD);
		for(int i=1;i<=k;i++)
			val-=Math.log(i);
		return val;
	}

	// (lambda^k) e^(-lambda) / k! 
	public double poisson(int k, double lambdaD) {
		double val =  Math.exp(-lambdaD)*Math.pow(lambdaD,k);
		for(int i=1;i<=k;i++) {
			val/=i;
		}
		return val;
	}
	
	// (n choose k) p^k (1-p)^{n-k}
	public static double logBinomial(int k, int n, double p) {
		double val;
		//System.out.println("k="+k+",n="+n+",p="+p);
		if((double)n*p <= 5 || (double)n*(1-p) <= 5) // compute exact binomial
			val = logbinomialcoefficient(n,k)+k*Math.log(p)+(n-k)*Math.log(1-p);
			//System.out.println("log binomial(k;n,p) = "+val);
		else // compute normal approximation: mu = np, sigmasq = np(1-p)
			val = logNormal(k,(double)n*p,(double)n*p*(1-p));
			//System.out.println("log normal(k;np,np(1-p)) = "+val);
			//System.out.println();
		return val;
	}

	// from http://www.lykkenborg.no/java/2006/03/binomial-coefficient.html
	// turned into logspace
	public static double logbinomialcoefficient(int n, int r) {
		double t = 0;

		int m = n - r; // r = Math.max(r, n - r);
		if (r < m) {
			r = m;
		}

		for (int i = n, j = 1; i > r; i--, j++) {
			t = t + Math.log(i) - Math.log(j);
		}

		return t;
	}

	// 1/(sigma * sqrt(2pi)) e^{-(x-mu)^2/(2sigmasq)}
	public static double logNormal(double x, double mu, double sigmasq) {
		double t = -Math.pow(x-mu,2)/(2*sigmasq) - (Math.log(2)+Math.log(Math.PI)+Math.log(sigmasq))/2;
		return t;
	}

	public static String getFragmentID(String esp) {
		if (esp2fragmentID.containsKey(esp))
			return esp2fragmentID.get(esp);
		return null;
	}
	
	public void getFragmentIDs() {
		fragmentIDs = new ArrayList<String>();
		string_iter = fragments.keySet().iterator();
		while(string_iter.hasNext())
			fragmentIDs.add(string_iter.next());
		Collections.sort(fragmentIDs);
	}


	public void writeFinalFiles() throws IOException {
		System.out.println("writing output files with prefix " + outputname);
		String bpfile;
		BufferedWriter writer;
		String assignfile;
		String cID;
		
		/*
		/// WRITE SHORT BURNIN
		bpfile = outputname+".finalbreakpoints.shortburnin";
		//System.out.println("writing breakpoints to " + bpfile);
		writer = new BufferedWriter(new FileWriter(bpfile));
		writer.write("ClusterID\tNumIters\tSupport0\tSupport1\tSupport2...\n");
		string_iter = breakpoints.keySet().iterator();
	
		while(string_iter.hasNext()) {
			cID = string_iter.next();
			writer.write(breakpoints.get(cID).cID+"\t"+(numiters-burnin));
			for(int j=0;j<breakpoints.get(cID).supportsshort.length;j++)
				writer.write("\t"+breakpoints.get(cID).supportsshort[j]);
			writer.write("\n");
		}
		writer.close();

		assignfile = outputname+".finalassignment.shortburnin";
		//System.out.println("writing fragment assignments to " + assignfile);
		writer = new BufferedWriter(new FileWriter(assignfile));
		writer.write("FragmentID\tAssignmentIndex\tAvgAssignment\tDiscordantPairs\n");
		int_iter = fragments.keySet().iterator();
		while(int_iter.hasNext()) {
			fragmentID = int_iter.next();
			writer.write(fragments.get(fragmentID).id+"\t-1\t"+fragments.get(fragmentID).avgNoassignmentshort+"\terror\n");
			for(int j=0;j<fragments.get(fragmentID).fragmentalignments.size();j++) {
				writer.write(fragments.get(fragmentID).id+"\t"+
						j+"\t"+
						fragments.get(fragmentID).avgAssignmentsshort.get(j)+"\t"+
						fragments.get(fragmentID).concatenateLabels(j)+"\n");
			}
		}
		writer.close();	
	*/
		
		/// WRITE LONG BURNIN
		bpfile = outputname+".finalbreakpoints.longburnin";
		//System.out.println("writing breakpoints to " + bpfile);
		writer = new BufferedWriter(new FileWriter(bpfile));
		writer.write("ClusterID\tNumIters\tSupport0\tSpport1\tSupport2...\n");
		string_iter = breakpoints.keySet().iterator();
		while(string_iter.hasNext()) {
			cID = string_iter.next();
			writer.write(breakpoints.get(cID).cID+"\t"+burnin);
			for(int j=0;j<breakpoints.get(cID).supportslong.length;j++)
				writer.write("\t"+breakpoints.get(cID).supportslong[j]);
			writer.write("\n");
		}
		writer.close();

		assignfile = outputname+".finalassignment.longburnin";
		//System.out.println("writing fragment assignments to " + assignfile);
		writer = new BufferedWriter(new FileWriter(assignfile));
		writer.write("FragmentID\tAssignmentIndex\tAvgAssignment\tDiscordantPairs\n");
		for(String fragmentID : fragmentIDs) {
			writer.write(fragments.get(fragmentID).fragmentID+"\t-1\t"+fragments.get(fragmentID).avgNoassignmentlong+"\terror\n");
			for(int j=0;j<fragments.get(fragmentID).fragmentalignments.size();j++) {
				writer.write(fragments.get(fragmentID).fragmentID+"\t"+
						j+"\t"+
						fragments.get(fragmentID).avgAssignmentslong.get(j)+"\t"+
						fragments.get(fragmentID).concatenateLabels(j)+"\n");
			}
		}
		writer.close();	

		/*
		/// WRITE THINNING
		bpfile = outputname+".finalbreakpoints.thinned";
		//System.out.println("writing breakpoints to " + bpfile);
		writer = new BufferedWriter(new FileWriter(bpfile));
		writer.write("ClusterID\tNumIters\tSupport0\tSpport1\tSupport2...\n");
		string_iter = breakpoints.keySet().iterator();
		while(string_iter.hasNext()) {
			cID = string_iter.next();
			writer.write(breakpoints.get(cID).cID+"\t"+(numiters / skip));
			for(int j=0;j<breakpoints.get(cID).supportsthinned.length;j++)
				writer.write("\t"+breakpoints.get(cID).supportsthinned[j]);
			writer.write("\n");
		}
		writer.close();

		assignfile = outputname+".finalassignment.thinned";
		//System.out.println("writing fragment assignments to " + assignfile);
		writer = new BufferedWriter(new FileWriter(assignfile));
		writer.write("FragmentID\tAssignmentIndex\tAvgAssignment\tDiscordantPairs\n");
		int_iter = fragments.keySet().iterator();
		while(int_iter.hasNext()) {
			fragmentID = int_iter.next();
			writer.write(fragments.get(fragmentID).id+"\t-1\t"+fragments.get(fragmentID).avgNoassignmentthinned+"\terror\n");
			for(int j=0;j<fragments.get(fragmentID).fragmentalignments.size();j++) {
				writer.write(fragments.get(fragmentID).id+"\t"+
						j+"\t"+
						fragments.get(fragmentID).avgAssignmentsthinned.get(j)+"\t"+
						fragments.get(fragmentID).concatenateLabels(j)+"\n");
			}
		}
		writer.close();	

		/// WRITE LAST
		bpfile = outputname+".finalbreakpoints.last";
		//System.out.println("writing breakpoints to " + bpfile);
		writer = new BufferedWriter(new FileWriter(bpfile));
		writer.write("ClusterID\tNumIters\tSupport0\tSpport1\tSupport2...\n");
		string_iter = breakpoints.keySet().iterator();
		while(string_iter.hasNext()) {
			cID = string_iter.next();
			writer.write(breakpoints.get(cID).cID+"\t1");
			for(int j=0;j<breakpoints.get(cID).supportslast.length;j++)
				writer.write("\t"+breakpoints.get(cID).supportslast[j]);
			writer.write("\n");
		}
		writer.close();

		assignfile = outputname+".finalassignment.last";
		//System.out.println("writing fragment assignments to " + assignfile);
		writer = new BufferedWriter(new FileWriter(assignfile));
		writer.write("FragmentID\tAssignmentIndex\tAvgAssignment\tDiscordantPairs\n");
		int_iter = fragments.keySet().iterator();
		while(int_iter.hasNext()) {
			fragmentID = int_iter.next();
			writer.write(fragments.get(fragmentID).id+"\t-1\t"+fragments.get(fragmentID).avgNoassignmentlast+"\terror\n");
			for(int j=0;j<fragments.get(fragmentID).fragmentalignments.size();j++) {
				writer.write(fragments.get(fragmentID).id+"\t"+
						j+"\t"+
						fragments.get(fragmentID).avgAssignmentlast.get(j)+"\t"+
						fragments.get(fragmentID).concatenateLabels(j)+"\n");
			}
		}
		writer.close();	
		*/
	}

	// just write short burnin.
	public void writeIntermediateFiles(int iter) throws IOException {
		System.out.println("writing intermediate files for iter " + iter);
		String bpfile = outputname+".iter"+iter+".finalbreakpoints.shortburnin";
		//System.out.println("writing breakpoints to " + bpfile);
		BufferedWriter writer = new BufferedWriter(new FileWriter(bpfile));
		writer.write("ClusterID\tNumIters\tSupport0\tSupport1\tSupport2...\n");
		string_iter = breakpoints.keySet().iterator();
		String cID;
		while(string_iter.hasNext()) {
			cID = string_iter.next();
			writer.write(breakpoints.get(cID).cID+"\t"+(iter-burnin));
			for(int j=0;j<breakpoints.get(cID).supportsshort.length;j++)
				writer.write("\t"+breakpoints.get(cID).supportsshort[j]);
			writer.write("\n");
		}
		writer.close();

		String assignfile = outputname+".iter"+iter+".finalassignment.shortburnin";
		System.out.println("writing fragment assignments to " + assignfile);
		writer = new BufferedWriter(new FileWriter(assignfile));
		//writer.write("FragmentID\tAssignmentIndex\tAvgAssignment\tDiscordantPairs\n");
		string_iter = fragments.keySet().iterator();
		String fragmentID;
		while(string_iter.hasNext()) {
			fragmentID = string_iter.next();
			writer.write(fragments.get(fragmentID).fragmentID+"\t-1\t"+(double)fragments.get(fragmentID).avgNoassignmentshort/(iter-burnin)+"\terror\n");
			for(int j=0;j<fragments.get(fragmentID).fragmentalignments.size();j++) {
				writer.write(fragments.get(fragmentID).fragmentID+"\t"+
						j+"\t"+
						(double)fragments.get(fragmentID).avgAssignmentsshort.get(j)/(iter-burnin)+"\t"+
						fragments.get(fragmentID).concatenateLabels(j)+"\n");
			}
		}
		writer.close();	
	}

	public void writeFirstFiles() throws IOException {
		String cID;
		String fragmentID;
		
		String bpfile = outputname+".finalbreakpoints.first";
		//System.out.println("writing breakpoints to " + bpfile);
		BufferedWriter writer = new BufferedWriter(new FileWriter(bpfile));
		writer.write("ClusterID\tNumIters\tSupport0\tSpport1\tSupport2...\n");
		string_iter = breakpoints.keySet().iterator();
		while(string_iter.hasNext()) {
			cID = string_iter.next();
			writer.write(breakpoints.get(cID).cID+"\t1");
			for(int j=0;j<breakpoints.get(cID).supportsfirst.length;j++)
				writer.write("\t"+breakpoints.get(cID).supportsfirst[j]);
			writer.write("\n");
		}
		writer.close();

		String assignfile = outputname+".finalassignment.first";
		//System.out.println("writing fragment assignments to " + assignfile);
		writer = new BufferedWriter(new FileWriter(assignfile));
		writer.write("FragmentID\tAssignmentIndex\tAvgAssignment\tDiscordantPairs\n");
		string_iter = fragments.keySet().iterator();
		while(string_iter.hasNext()) {
			fragmentID = string_iter.next();
			writer.write(fragments.get(fragmentID).fragmentID+"\t-1\t"+fragments.get(fragmentID).avgNoassignmentfirst+"\terror\n");
			for(int j=0;j<fragments.get(fragmentID).fragmentalignments.size();j++) {
				writer.write(fragments.get(fragmentID).fragmentID+"\t"+
						j+"\t"+
						fragments.get(fragmentID).avgAssignmentfirst.get(j)+"\t"+
						fragments.get(fragmentID).concatenateLabels(j)+"\n");
			}
		}
		writer.close();	
	}

	public void writeStateCounts(ArrayList<String> keys) throws IOException {
		String countfile = outputname+".samplingprobs";
		BufferedWriter writer = new BufferedWriter(new FileWriter(countfile));
		writer.write("SampledProb\tLogProb");
		for(int i=0;i<keys.size();i++) {
			writer.write("\t"+keys.get(i));
		}
		writer.write("\n");
		string_iter = counts.keySet().iterator();
		String key;
		while(string_iter.hasNext()) {
			key = string_iter.next();
			writer.write(counts.get(key)/(double)numiters+"\t"+logprobs.get(key)+"\t"+key+"\n");
		}
		writer.close();
	}
	
	
}