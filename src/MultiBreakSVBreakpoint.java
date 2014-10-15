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
