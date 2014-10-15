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
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;

/**
 * Node for traversing diagram. Not used to build cluster diagram.
 * @author aritz
 *
 */
public class DiagramNode {
	public String cID;
	public String[] esps;
	public ArrayList<DiagramNode> parents;
	public ArrayList<DiagramNode> children;
	public DiagramNode selectedparent;

	public DiagramNode(String c,String[] e) {
		cID = c;
		esps = e;
		parents = new ArrayList<DiagramNode>();
		children = new ArrayList<DiagramNode>();
		selectedparent = null;
	}

	public DiagramNode(String c,String e) {
		cID = c;
		esps = new String[1];
		esps[0] = e;
		parents = new ArrayList<DiagramNode>();
		children = new ArrayList<DiagramNode>();
		selectedparent = null;
	}

	public void print() {
		System.out.print(cID+" ("+parents.size()+" parents and " + children.size() + " children): ");
		for(int i=0;i<esps.length;i++)
			System.out.print(esps[i]+" ");
		System.out.println();
	}

	public void printLong() {
		System.out.print(cID+": ");
		for(int i=0;i<esps.length;i++)
			System.out.print(esps[i]+" ");
		System.out.println();
		System.out.print("  parents: ");
		for(int i=0;i<parents.size();i++) {
			System.out.print(parents.get(i).cID+" ");
		}
		System.out.println();
		System.out.print("  children: ");
		for(int i=0;i<children.size();i++) {
			System.out.print(children.get(i).cID+" ");
		}
		System.out.println();
	}

	public boolean contains(String s) {
		for(int i=0;i<esps.length;i++)
			if(esps[i].equals(s))
				return true;
		return false;
	}

	public String espNames() {
		String str = esps[0];
		for(int i=1;i<esps.length;i++)
			str=str+","+esps[i];
		return str;
	}

}
