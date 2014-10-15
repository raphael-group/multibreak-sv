/*
Copyright 2014 Brown University, Providence, RI.

                         All Rights Reserved

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose other than its incorporation into a
commercial product is hereby granted without fee, provided that the
above copyright notice appear in all copies and that both that
copyright notice and this permission notice appear in supporting
documentation, and that the name of Brown University not be used in
advertising or publicity pertaining to distribution of the software
without specific, written prior permission.

BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
http://cs.brown.edu/people/braphael/software.html

If you use MultiBreak-SV in your research, please cite:

A. Ritz, A. Bashir, S. Sindi, D. Hsu, I. Hajirasouliha, and B. J. Raphael. Characterization of Structural Variants with Single Molecule and Hybrid Sequencing Approaches, under review.
*/
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

public class MultiBreakSVAssignmentMatrix {
	public HashMap<String,HashMap<Experiment,ArrayList<Integer>>> supports; // <cid,<Experiment,supports>>
	public HashMap<String,Integer> assignments; // <fragment,assignment>
	public HashMap<String,Double> logprobcov; // <cID,logprobcov>
	public HashMap<Experiment,Integer> numerrors; // Experiment, numerors
	public HashMap<Experiment,Integer> len; // Experiment, len

	public ArrayList<String> keys; // sorted strobe indices

	public double prior;
	public double logprob;

	public Iterator<String> string_iter;
	public Iterator<Integer> int_iter;

	public MultiBreakSVAssignmentMatrix() {

		supports = new HashMap<String,HashMap<Experiment,ArrayList<Integer>>>();
		logprobcov = new HashMap<String,Double>();
		string_iter = MultiBreakSV.breakpoints.keySet().iterator();
		String cID;
		String explabel;
		Iterator<String> exp_iter;
		while(string_iter.hasNext()) {
			cID = string_iter.next();
			supports.put(cID,new HashMap<Experiment,ArrayList<Integer>>());
			exp_iter = MultiBreakSV.experiments.keySet().iterator();
			while(exp_iter.hasNext()) {
				explabel = exp_iter.next();
				supports.get(cID).put(MultiBreakSV.experiments.get(explabel),new ArrayList<Integer>());
			}
			logprobcov.put(cID,0.0);
		}

		assignments = new HashMap<String,Integer>();
		string_iter = MultiBreakSV.fragments.keySet().iterator();
		String fragmentID;
		while(string_iter.hasNext()) {
			fragmentID = string_iter.next();
			assignments.put(fragmentID,-1);
		}

		string_iter = assignments.keySet().iterator();
		keys = new ArrayList<String>();
		while(string_iter.hasNext()) {
			keys.add(string_iter.next());
		}
		Collections.sort(keys);

		numerrors  = new HashMap<Experiment,Integer>();
		len = new HashMap<Experiment,Integer>();
	}

	public MultiBreakSVAssignmentMatrix(MultiBreakSVAssignmentMatrix A) {
		prior = A.prior;
		logprob = A.logprob;

		supports = new HashMap<String,HashMap<Experiment,ArrayList<Integer>>>();
		logprobcov = new HashMap<String,Double>();
		string_iter = MultiBreakSV.breakpoints.keySet().iterator();
		String cID;
		Experiment exp;
		Iterator<String> exp_iter;
		while(string_iter.hasNext()) {
			cID = string_iter.next();
			supports.put(cID,new HashMap<Experiment,ArrayList<Integer>>());
			exp_iter = MultiBreakSV.experiments.keySet().iterator();
			while(exp_iter.hasNext()) {
				exp = MultiBreakSV.experiments.get(exp_iter.next());
				supports.get(cID).put(exp,new ArrayList<Integer>());
				for(int j=0;j<A.supports.get(cID).get(exp).size();j++)
					supports.get(cID).get(exp).add(A.supports.get(cID).get(exp).get(j));
			}
			logprobcov.put(cID,A.logprobcov.get(cID));
		}


		assignments = new HashMap<String,Integer>();
		string_iter = MultiBreakSV.fragments.keySet().iterator();
		String fragmentID;
		while(string_iter.hasNext()) {
			fragmentID = string_iter.next();
			assignments.put(fragmentID,A.assignments.get(fragmentID));
		}

		string_iter = assignments.keySet().iterator();
		keys = new ArrayList<String>();
		while(string_iter.hasNext()) {
			keys.add(string_iter.next());
		}
		Collections.sort(keys);

		numerrors = A.numerrors;
		len = A.len;
	}

	public void computeLogProb(boolean verbose) {
		logprob = 0;

		// compute log prob of SV coverage
		if(verbose)
			System.out.println("Log Prob of SV Coverage:");

		string_iter = supports.keySet().iterator();
		String cID;
		int numbreakpoints = 0;
		int thissupport;
		Iterator<String> exp_iter;
		Experiment exp;
		while(string_iter.hasNext()) {

			cID = string_iter.next();

			logprobcov.put(cID,0.0);
			int maxbreakpoints = 0;
			exp_iter = MultiBreakSV.experiments.keySet().iterator();
			while(exp_iter.hasNext()) {

				exp = MultiBreakSV.experiments.get(exp_iter.next());

				if(supports.get(cID).get(exp).size() > maxbreakpoints)
					maxbreakpoints = supports.get(cID).get(exp).size();

				for(int j=0;j<supports.get(cID).get(exp).size();j++) {
					thissupport = supports.get(cID).get(exp).get(j);
					if(verbose) 
						System.out.println("ExperimentLabel: " + exp.label + " with LambdaD = " + 
								exp.lambdaD+ " :" + cID + " #" + j + "\t: support "+
								thissupport + "\t probability " + exp.logprobcov[thissupport]);


					// Compute LogProb Coverage with no combinatorial coefficient
					logprobcov.put(cID,logprobcov.get(cID) + 
							exp.logprobcov[thissupport]);
				}
			}
			numbreakpoints+=maxbreakpoints;
			logprob+=logprobcov.get(cID);
			if(logprobcov.get(cID) > 1) {
				System.out.println("ERROR! logprobcov.get("+cID+") = "+logprobcov.get(cID));
				throw new Error();
			}
		}

		// compute prior
		if(verbose)
			System.out.println("\nLogPrior\t" + MultiBreakSV.prior[numbreakpoints]+"\n");

		prior = MultiBreakSV.prior[numbreakpoints];
		logprob+=prior;
		if(prior > 1) {
			System.out.println("ERROR! log prior = "+prior);
			throw new Error();
		}

		// compute log prob of alignments
		if(verbose)
			System.out.println("Log Prob of Alignments");

		numerrors  = new HashMap<Experiment,Integer>();
		len = new HashMap<Experiment,Integer>();
		exp_iter = MultiBreakSV.experiments.keySet().iterator();
		String explabel;
		while(exp_iter.hasNext()) {
			explabel = exp_iter.next();
			numerrors.put(MultiBreakSV.experiments.get(explabel),0);
			len.put(MultiBreakSV.experiments.get(explabel),0);
		}

		string_iter = assignments.keySet().iterator();
		String fragmentID;
		Fragment f;
		FragmentAlignment fa;
		int numAssignedToError = 0;	
		while(string_iter.hasNext()) {
			fragmentID = string_iter.next();
			if(assignments.get(fragmentID) == -1) {
				numAssignedToError++;
				f = MultiBreakSV.fragments.get(fragmentID);
				fa = f.fragmentalignments.get(0); // all strobes have at least one alignment.

			} else {
				f = MultiBreakSV.fragments.get(fragmentID);
				fa = f.fragmentalignments.get(assignments.get(fragmentID));
				numerrors.put(fa.exp,numerrors.get(fa.exp)+fa.numerrors);
				len.put(fa.exp,len.get(fa.exp)+fa.len);		
			}
		}

		logprob+= numAssignedToError*Math.log(MultiBreakSV.perr);
		
		exp_iter = MultiBreakSV.experiments.keySet().iterator();
		while(exp_iter.hasNext()) {
			explabel = exp_iter.next();
			logprob+=MultiBreakSV.logBinomial(numerrors.get(MultiBreakSV.experiments.get(explabel)),
					len.get(MultiBreakSV.experiments.get(explabel)),MultiBreakSV.experiments.get(explabel).pseq);
		}
		
		if(verbose) {
			System.out.println(numAssignedToError+" assigned to errors\n");
			exp_iter = MultiBreakSV.experiments.keySet().iterator();
			while(exp_iter.hasNext()) {
				explabel = exp_iter.next();
				System.out.println(explabel+": lambdaD=" + MultiBreakSV.experiments.get(explabel).lambdaD + 
						", pseq=" + MultiBreakSV.experiments.get(explabel).pseq +  
						", numerrors="+numerrors.get(MultiBreakSV.experiments.get(explabel))+
						", length="+len.get(MultiBreakSV.experiments.get(explabel))+
						"\n  log B(numerrors,length,pseq) = " + MultiBreakSV.logBinomial(numerrors.get(MultiBreakSV.experiments.get(explabel)),
								len.get(MultiBreakSV.experiments.get(explabel)),MultiBreakSV.experiments.get(explabel).pseq));
			}
			System.out.println();
		}
	}

	public void setBreakpoints(HashMap<String,MultiBreakSVDiagram> diagrams) {
		MultiBreakSVDiagram d;
		Iterator<String> iter = diagrams.keySet().iterator();
		while(iter.hasNext()) {
			d = diagrams.get(iter.next());
			setBreakpoints(d,false);
		}
		d = null;
	}

	public void setBreakpoints(MultiBreakSVDiagram diagram) {
		setBreakpoints(diagram,false);
	}

	// sets the breakpoints for selected alignments in diagram.
	// uses approximation of set cover.
	public void setBreakpoints(MultiBreakSVDiagram diagram,boolean verbose){

		if(verbose && diagram.isCluster)
			System.out.println("Setting breakpoints for cluster " + diagram.cIDroot);
		if(verbose && !diagram.isCluster) {
			System.out.println("Setting breakpoints for conn comp "+ diagram.cIDroot);
			diagram.print();
		}

		////////////////////////////
		// Step 1: Get ESPs that are selected.
		////////////////////////////
		HashMap<Experiment,ArrayList<String>> esps = new HashMap<Experiment,ArrayList<String>>();
		String explabel;
		Iterator<String> exp_iter = MultiBreakSV.experiments.keySet().iterator();
		while(exp_iter.hasNext()) {
			explabel = exp_iter.next();
			//System.out.println("LAMBDAD: "+MultiBreakSV.lambdaDs.get(l));
			esps.put(MultiBreakSV.experiments.get(explabel),new ArrayList<String>());
		}

		String key;
		FragmentAlignment fa;
		String fragmentID;
		Iterator<String> iter = diagram.leaves.keySet().iterator();
		//diagram.print();
		int numselected = 0;
		while(iter.hasNext()) {
			key = iter.next();

			// first, get index in strobe array for this subread alignment.
			fragmentID = MultiBreakSV.getFragmentID(key);

			// if strobeID is not present, it must have been removed during graph simp.
			if(!MultiBreakSV.fragments.containsKey(fragmentID))
				continue;

			// if this strobe has an alignment, see if this subread alignment is in it.
			if (assignments.get(fragmentID) != -1) {
				//System.out.println(strobeID+":"+MultiBreakSV.strobes.get(strobeID).strobealignments.size() + " and " + assignments.get(strobeID));
				//System.out.flush();
				fa = MultiBreakSV.fragments.get(fragmentID).fragmentalignments.get(assignments.get(fragmentID));
				for(int j=0;j<fa.alignments.size();j++) {
					if (key.equals(fa.alignments.get(j))) {
						esps.get(fa.exp).add(key);
						numselected++;
						if(verbose)
							System.out.println("  adding esp " + key);
					}
				}
			}
		}

		// if there are no selected ESPs, set support to empty and return.
		if(numselected==0) {
			if(!supports.containsKey(diagram.cIDroot)) {
				System.out.println(diagram.cIDroot + " is missing");
				throw new Error();
			}
			supports.put(diagram.cIDroot,new HashMap<Experiment,ArrayList<Integer>>());
			exp_iter = MultiBreakSV.experiments.keySet().iterator();
			while(exp_iter.hasNext()) {
				explabel = exp_iter.next();
				supports.get(diagram.cIDroot).put(MultiBreakSV.experiments.get(explabel),new ArrayList<Integer>());
			}
			return;
		}

		////////////////////////////
		// Step 2: Iterative Set Cover Approximation
		////////////////////////////

		// DO FOR EACH LAMBDAD!
		supports.put(diagram.cIDroot,new HashMap<Experiment,ArrayList<Integer>>());
		exp_iter = MultiBreakSV.experiments.keySet().iterator();
		Experiment exp;
		while(exp_iter.hasNext()) {
			exp = MultiBreakSV.experiments.get(exp_iter.next());
			boolean found[] = new boolean[esps.get(exp).size()];
			Arrays.fill(found,false);
			ArrayList<Integer> diagramsup = new ArrayList<Integer>();
			ArrayList<DiagramNode> maximalclusters = new ArrayList<DiagramNode>();
			for(int i=0;i<diagram.roots.size();i++)
				maximalclusters.add(diagram.roots.get(i));

			boolean allfound = false;
			int maxcount,curcount;
			int maxind = -1;
			while(!allfound) {
				// get maximal cluster with largest number of undiscovered esps.
				maxcount = 0;
				maxind = -1;
				for(int i=0;i<maximalclusters.size();i++) {
					curcount = 0;
					for(int j=0;j<maximalclusters.get(i).esps.length;j++) {
						String str = maximalclusters.get(i).esps[j];
						if(esps.get(exp).contains(str) && !found[esps.get(exp).indexOf(str)])
							curcount++;
					}
					if(curcount>maxcount) {
						maxcount = curcount;
						maxind = i;
					}
				}

				// set the corresponding ESPs to 'found'
				if(maxind == -1) { 
					allfound = true;
					for(int i=0;i<found.length;i++)
						if(!found[i])
							allfound = false;
					if(!allfound)
						throw new Error();
					continue;
				}
				for(int j=0;j<maximalclusters.get(maxind).esps.length;j++) {
					String str = maximalclusters.get(maxind).esps[j];
					if(esps.get(exp).contains(str)) // it's ok if found[] is already true - just overwrite.
						found[esps.get(exp).indexOf(str)] = true;
				}

				// add the count to diagramsup
				diagramsup.add(maxcount);
				if(verbose)
					System.out.println("added support of " + maxcount);

				// remove the maximal cluster from consideration.
				maximalclusters.remove(maxind);

				// see if all have been found.
				allfound = true;
				for(int i=0;i<found.length;i++)
					if(!found[i])
						allfound = false;
			}

			////////////////////////////
			// Step 3: replace supports with diagramsup
			////////////////////////////
			Collections.sort(diagramsup);
			if(!supports.containsKey(diagram.cIDroot)) {
				System.out.println(diagram.cIDroot + " is missing");
				throw new Error();
			}

			supports.get(diagram.cIDroot).put(exp,diagramsup);
		}

		if(verbose)
			System.out.println("done setting breakpoints.");
	}

	public void removeESP(String esp, DiagramNode n, MultiBreakSVDiagram diagram) {
		// find index that matches this esp.
		int ind = -1;
		for(int i=0;i<n.esps.length;i++) {
			if(esp.equals(n.esps[i])) {
				ind = i;
				break;
			}
		}

		// if index exists, make new esp array and reset the node.
		if (ind > -1) { 
			String[] newarray = new String[n.esps.length-1];
			for(int i=0;i<ind;i++)
				newarray[i] = n.esps[i];
			for(int i=ind+1;i<n.esps.length;i++)
				newarray[i-1] = n.esps[i];

			// reset
			n.esps = newarray;
		}

		// remove node if there are no more esps
		if(n.esps.length == 0) {
			// remove children link
			for(int i=0;i<n.parents.size();i++) {
				n.parents.get(i).children.remove(n);
			}

			// remove from nodes
			diagram.nodes.remove(n);

			// DO NOT remove from leaves if it's a leaf YET
			// (loop that calls this iterates through leaves)

			// remove from roots if it's a root
			if (diagram.roots.contains(n))
				diagram.roots.remove(n);
		}

		// recurse for parents.
		for(int i=0;i<n.parents.size();i++)
			removeESP(esp,n.parents.get(i),diagram);
	}

	public boolean sameAs(MultiBreakSVAssignmentMatrix A) {
		return (this.summaryString().equals(A.summaryString()));	
	}

	public String summaryString() {

		String str = ""+assignments.get(keys.get(0));
		for(int i=1;i<keys.size();i++)
			str+="\t"+assignments.get(keys.get(i));

		return str;
	}

	public void printSupport() {
		System.out.print("\nSUPPORT:\n");
		string_iter = supports.keySet().iterator();
		String cID;
		while(string_iter.hasNext()) {
			cID = string_iter.next();
			Iterator<Experiment> iter = supports.get(cID).keySet().iterator();
			Experiment exp;
			while(iter.hasNext()) {
				exp = iter.next();
				if(supports.get(cID).get(exp).size() > 0) {
				System.out.print("  "+cID+": "+exp.label+" lambdaD="+exp.lambdaD+":("+supports.get(cID).get(exp).get(0));
				for(int j=1;j<supports.get(cID).get(exp).size();j++)
					System.out.print(","+supports.get(cID).get(exp).get(j));
				System.out.print(")\n");
				
				}
			}
		}
	}
	
	public void printAssignments() {
		System.out.println("ASSIGNMENTS:");
		string_iter = assignments.keySet().iterator();
		String fragmentID;
		FragmentAlignment fa;
		while(string_iter.hasNext()) {
			fragmentID = string_iter.next();
			System.out.print("  "+fragmentID+": "+assignments.get(fragmentID)+" --> ");
			if (assignments.get(fragmentID)==-1)
				System.out.println("error");
			else {
				fa = MultiBreakSV.fragments.get(fragmentID).fragmentalignments.get(assignments.get(fragmentID));
				System.out.print(fa.alignments.get(0));
				for (int i=1;i<fa.alignments.size();i++) 
					System.out.print(", "+fa.alignments.get(i));
				System.out.println(": " + fa.numerrors + " errors and " + fa.len + " length.");
			}
		}
	}
}
