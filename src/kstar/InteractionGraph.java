package kstar;
import java.util.LinkedHashSet;


public class InteractionGraph {

	private int v[] = null; //the vertex set (molecule-relative residue numbering, not pdb numbering)
	private boolean e[][] = null; //the edge set
	private float maxEnergy[][] = null; // SJ, added for keeping the max abs energy value between pair of residues
	private float minDistance[][] = null; // SJ, added for keeping track of the min distance between pair of residues
	
	private int numAddedV = 0; //the number of vertices currently added to v[]
	
	InteractionGraph(int numV){
		
		v = new int[numV];
		e = new boolean[numV][numV];
		maxEnergy = new float[numV][numV]; // SJ
		minDistance = new float[numV][numV];
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
	
	//SJ, Deletes the edge between the vertices x and y
	public void deleteEdge (int x, int y){
		
		int xInd = findVind(x);
		int yInd = findVind(y);
		
		e[xInd][yInd] = false;
		e[yInd][xInd] = false;
	}
	
	// SJ, adds the minDistance and maxEnergy for vertices x and y
	public void addDistEner (int x, int y, float dist, float ener){
		
		int xInd = findVind(x);
		int yInd = findVind(y);
		
		minDistance[xInd][yInd] = dist;
		minDistance[yInd][xInd] = dist;
		maxEnergy[xInd][yInd] = ener;
		maxEnergy[yInd][xInd] = ener;
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
	
	public String toString()
	{ 
		String out = "";
		for(int i = 0; i < e.length; i++)
		{
			for(int j = 0; j < e[i].length; j++)
			{
				if(e[i][j])
					out+= "("+i+","+j+")\n";
			}
		}
		return out;
	}
}
