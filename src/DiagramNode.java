

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
