import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 * class for traversing diagram. Not used to build cluster diagram.
 * @author aritz
 *	// http://swjscmail1.java.sun.com/cgi-bin/wa?A2=ind9909&L=rmi-users&F=&S=&P=34494
 * for serializable info.
 */
public class MultiBreakSVDiagram {
	public ArrayList<DiagramNode> nodes;
	public HashMap<String,DiagramNode> leaves;
	public ArrayList<DiagramNode> roots;
	public HashMap<String,ArrayList<Integer>> faInds;
	public Boolean isCluster;
	public Integer maxindegree;
	public String cIDroot;

	public MultiBreakSVDiagram() {
		nodes = new ArrayList<DiagramNode>();
		leaves = new HashMap<String,DiagramNode>();
		roots = new  ArrayList<DiagramNode>();
		faInds = new HashMap<String,ArrayList<Integer>>();
		maxindegree = 0;
	}

	public void constructFromCluster(String cID, String[] names) {
		isCluster = true;
		cIDroot = cID;

		DiagramNode parent = new DiagramNode(cID,names);
		nodes.add(parent);
		roots.add(parent);

		DiagramNode n;
		int rCounter = 0;
		for(int i=0;i<names.length;i++) { 
			n = new DiagramNode(cID+".R"+rCounter,names[i]);
			rCounter++;
			n.parents.add(parent);
			parent.children.add(n);

			nodes.add(n);
			leaves.put(names[i],n);

		}

		maxindegree = parent.children.size(); // max indegree is at root.
		setSAinds();
	}

	public void constructFromConnComp(String cIDr,ArrayList<String[]> maximals,ArrayList<String> cIDs, boolean setSAs, boolean verbose) throws FileNotFoundException {
		isCluster = false;
		cIDroot = cIDr;

		HashMap<String,Boolean> seen = new HashMap<String,Boolean>();
		DiagramNode n;
		
		if (verbose) {
			System.out.println(cIDs.size() + " cIDs");
			for (String c : cIDs)
				System.out.println("CID " + c);
				
		}
		
		// make leaves.
		for(int i=0;i<maximals.size();i++) {
			for(int j=0;j<maximals.get(i).length;j++) {
				if(!seen.containsKey(maximals.get(i)[j])) {
					n = new DiagramNode(cIDroot,maximals.get(i)[j]);
					leaves.put(maximals.get(i)[j],n);
					nodes.add(n);
				}
			}
		}
		if (verbose)
			System.out.println(leaves.size() + " leaves");

		// make roots.
		DiagramNode child;
		for(int i=0;i<maximals.size();i++) {
			n = new DiagramNode(cIDroot,maximals.get(i));
			// add edges
			for(int j=0;j<maximals.get(i).length;j++) {
				child = leaves.get(maximals.get(i)[j]);
				n.children.add(child);
				child.parents.add(n);
			}
			roots.add(n);
			nodes.add(n);
		}
		if (verbose)
			System.out.println(roots.size() + " roots");
		if(verbose)
			System.out.println("  getting max indegree...");
		for(int i=0;i<nodes.size();i++) {

			// find max indegree.
			if(nodes.get(i).children.size() > maxindegree) 
				maxindegree = nodes.get(i).children.size();
		}

	

		if(setSAs) {
			if(verbose)
				System.out.println("  setting SAinds...");
			setSAinds();
		}
	}

	// sets the strobe alignment indices for each leaf.
	public void setSAinds() {
		Iterator<String> iter  = leaves.keySet().iterator();
		String fragmentID;
		Fragment f;
		String esp;
		while(iter.hasNext()) {
			esp = iter.next();
			faInds.put(esp,new ArrayList<Integer>());

			fragmentID = MultiBreakSV.getFragmentID(esp);

			// if the diagram has an ESP not in the clusters file, ignore.
			if(!MultiBreakSV.fragments.containsKey(fragmentID))
				continue;

			f = MultiBreakSV.fragments.get(fragmentID);

			for(int i=0;i<f.fragmentalignments.size();i++)
				if(f.fragmentalignments.get(i).contains(esp))
					faInds.get(esp).add(i);	
		}
	}

	public void print() {
		System.out.println("Diagram cIDroot = " + cIDroot);
		System.out.println(leaves.size() + " leaves and " + roots.size() + " roots");
		
	}


}
