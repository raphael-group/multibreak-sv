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

public class Fragment {

	public String fragmentID;
	public ArrayList<FragmentAlignment> fragmentalignments;
	public ArrayList<Double> avgAssignmentsshort,avgAssignmentslong,avgAssignmentsthinned,avgAssignmentlast,avgAssignmentfirst;
	public double avgNoassignmentshort,avgNoassignmentlong,avgNoassignmentthinned,avgNoassignmentlast,avgNoassignmentfirst;

	public Fragment(String id) {
		fragmentID = id;
		fragmentalignments = new ArrayList<FragmentAlignment>();
		initLists();
	}
	
	public void initLists() {
		avgAssignmentsshort = new ArrayList<Double>();
		avgNoassignmentshort = 0;
		avgAssignmentslong = new ArrayList<Double>();
		avgNoassignmentlong = 0;
		avgAssignmentsthinned = new ArrayList<Double>();
		avgNoassignmentthinned = 0;
		avgAssignmentlast = new ArrayList<Double>();
		avgNoassignmentlast = 0;
		avgAssignmentfirst = new ArrayList<Double>();
		avgNoassignmentfirst = 0;
	}
	
	public void print() {
		System.out.println("fragment id = "+ fragmentID + ", "+fragmentalignments.size() + " strobe alignments.");
		for(int i=0;i<fragmentalignments.size();i++) {
			System.out.println(" #"+i+": "+fragmentalignments.get(i).alignments.size() + " alignments");
			for(int j=0;j<fragmentalignments.get(i).alignments.size();j++) {
				System.out.println("  #"+j+" "+fragmentalignments.get(i).alignments.get(j));
			}
		}
	}
	
	public void addFragmentAlignment(FragmentAlignment a) {
		fragmentalignments.add(a);
		avgAssignmentsshort.add((double)0);
		avgAssignmentslong.add((double)0);
		avgAssignmentsthinned.add((double)0);
		avgAssignmentlast.add((double)0);
		avgAssignmentfirst.add((double)0);
	}

	public String concatenateLabels(int i) {
		String str = "";
		FragmentAlignment fa = fragmentalignments.get(i);
		
		if(fa.alignments.size() > 0) 
			str= str+""+fa.alignments.get(0);
		for(int j=1;j<fa.alignments.size();j++) {
			str= str+","+fa.alignments.get(j);
		}
		
		return str;
	}
	
	public ArrayList<String> getAllAlignments() {
		ArrayList<String> alignments = new ArrayList<String>();
		for(int i=0;i<fragmentalignments.size();i++) {
			alignments.addAll(fragmentalignments.get(i).alignments);
		}
		return alignments;
	}

	
}
