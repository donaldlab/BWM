package kstar;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.StringTokenizer;
import BranchDecomposition.BranchNode;
import BranchDecomposition.BranchTree;
import BranchDecomposition.GraphVertices;


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
	
	public static InteractionGraph loadFromFile(String fileName)
	{
		ArrayList<Integer> leftColumn = new ArrayList<>();
		ArrayList<Integer> rightColumn = new ArrayList<>();
		BufferedReader bufread = null;
		try {
			File file = new File(fileName);		
			FileReader fr = new FileReader(file);  
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			System.out.println("ERROR: Residue interaction graph file not found");
			System.exit(1);
		}
		
		String str = null;
		boolean done = false;
		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
				System.out.println("ERROR: An error occurred while reading input");
				System.exit(1);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			
			else if ( (!getToken(str,1).equalsIgnoreCase("PIG:0")) && (!str.equalsIgnoreCase("0 0")) ) { //not the first/last lines in file
				
				Integer v1 = Integer.valueOf(getToken(str,1));
				Integer v2 = Integer.valueOf(getToken(str,2));
				leftColumn.add(v1);
				rightColumn.add(v2);
				
			}
		}
		

		InteractionGraph graph = new InteractionGraph(leftColumn.size());
		for(int i = 0; i < leftColumn.size(); i++)
		{
			Integer leftVertex = leftColumn.get(i);
			graph.addV(leftVertex);
			Integer rightVertex = rightColumn.get(i);
			graph.addV(rightVertex);
			
			graph.addEdge(leftVertex, rightVertex);
		}
		
		try { bufread.close(); } catch(Exception e){} //done, so close file for reading
		
		return graph;
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
	
	// This function returns the xth token in string s
	private static String getToken(String s, int x) {
	
		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");
		
		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
			  st.nextToken();
			else {
				System.out.println("ERROR: Unable to access argument " + x + " from input string");
				System.exit(1);
			}
		}
		
		if (st.hasMoreTokens())		
			return(st.nextToken());
		
		else {
			System.out.println("ERROR: Unable to access argument " + x + " from input string");
			System.exit(1);
			return null;
		}

	} // end getToken
}
