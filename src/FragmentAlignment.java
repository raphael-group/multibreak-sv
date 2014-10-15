

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
