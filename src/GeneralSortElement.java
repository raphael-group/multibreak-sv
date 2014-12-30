public class GeneralSortElement implements Comparable<GeneralSortElement>{

	public String read;
	public double perc;
	public String row;
	
	public GeneralSortElement(String r,String s,double a) {
		read = s;
		perc = a;
		row = r;
	}
	public int compareTo(GeneralSortElement se) {
		if (this.read.compareTo(se.read) < 0) 
			return -1;
		if (this.read.compareTo(se.read) > 0)
			return 1;

		// strobes are the same; check percent ID
		if(!GeneralExternalSort.reverse) {
			if (this.perc > se.perc)
				return -1;
			if (this.perc < se.perc)
				return 1;
		}
		else {
			if (this.perc > se.perc)
				return 1;
			if (this.perc < se.perc)
				return -1;	
		}

		// strobes and scores are the same; they are equal
		return 0;
	}
}
