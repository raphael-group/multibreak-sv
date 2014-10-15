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

public class FragmentAlignment implements Comparable<FragmentAlignment>{

	public ArrayList<String> alignments;
	public int numerrors;
	public int len;
	public Experiment exp;
	
	public FragmentAlignment(Experiment e, String[] a, int ne, int l) {
		alignments = new ArrayList<String>();
		for (int i=0;i<a.length;i++)
			alignments.add(a[i]);
		exp = e;
		numerrors = ne;
		len = l;
	}
		
	public void print() {
		System.out.print("Fragment Alignment has " + alignments.size() + " alignments: ");
		for(int j=0;j<alignments.size();j++) {
			System.out.print(alignments.get(j)+" ");
		}
		System.out.println();
	}

	public boolean equals(FragmentAlignment A) {
		return contains(A) && A.contains(this);
	}
	
	public boolean contains(FragmentAlignment A) {
		for(int i=0;i<A.alignments.size();i++) {
			boolean match = false;
			for(int j=0;j<alignments.size();j++) {
				if (A.alignments.get(i).equals(alignments.get(j))) {
					match = true;
					break;
				}
			}
			if (!match)
				return false;
		}
		return true;
	}

	public boolean contains(String esp) {
		for(int i=0;i<alignments.size();i++) 
			if(alignments.get(i).equals(esp))
				return true;
		return false;
	}

	@Override
	public int compareTo(FragmentAlignment o) {
		if(numerrors < o.numerrors)
			return 1;
		if(numerrors > o.numerrors)
			return -1;
		return 0;
		
	}
}
