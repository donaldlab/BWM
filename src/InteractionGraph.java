import java.util.LinkedHashSet;


public class InteractionGraph {

	private int v[] = null; //the vertex set (molecule-relative residue numbering, not pdb numbering)
	private boolean e[][] = null; //the edge set
	
	private int numAddedV = 0; //the number of vertices currently added to v[]
	
	InteractionGraph(int numV){
		
		v = new int[numV];
		e = new boolean[numV][numV];
		numAddedV = 0;
	}
	
	//Adds the vertex with molecule-relative number molResNum to the set of vertices for this graph
	public void addV(int molResNum){
		
		if (findVind(molResNum)>=0)
			return;
		
		v[numAddedV] = molResNum;
		numAddedV++;
	}
	
	//Adds an edge between vertices x and y (molecule-relative residue numbering)
	public void addEdge (int x, int y){
		
		int xInd = findVind(x);
		int yInd = findVind(y);
		
		e[xInd][yInd] = true;
		e[yInd][xInd] = true;
	}
	
	//Determines is an edge exists between vertices x and y (molecule-relative residue numbering)
	public boolean edgeExists(int x, int y){
		
		return e[findVind(x)][findVind(y)];
	}
	
	//Finds the index into v[] for the vertex with molecule-relative number molResNum; returns -1 if vertex not found
	private int findVind(int molResNum){
		
		for (int i=0; i<numAddedV; i++){
			if (molResNum==v[i]) //found vertex
				return i;
		}
		
		return -1;
	}
	
	//Returns all vertices of this graph
	public LinkedHashSet<Integer> getAllV(){
		
		LinkedHashSet<Integer> gv = new LinkedHashSet<Integer>();
		for (int i=0; i<numAddedV; i++)
			gv.add(v[i]);
		
		return gv;
	}
}
