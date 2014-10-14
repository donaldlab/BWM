package kstar;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.Serializable;
import java.util.LinkedHashSet;
import java.util.StringTokenizer;

import BDAStar.BWMAStarNode;
import BDAStar.BWMSolutionSpace;

public class BranchTree implements Serializable {

	private TreeNode root = null; //the root of the tree
	
	private InteractionGraph G = null; //the residue interaction graph
	
	/* Enumeration objects */
	private static BWMAStarNode AStarRoot;
	private static BWMSolutionSpace solutionSpace;
	
	public BranchTree(String fName, Molecule m, int numUnprunedRot[], int molResidueMap[], int invResidueMap[], int sysStrandNum, int numResInAS, boolean ligPresent){
		
		int numV = numResInAS;
		if (ligPresent) //ligand is present
			numV++;
		
		G = new InteractionGraph(numV);
		
		readBranchDecomposition(fName, m, numUnprunedRot, molResidueMap, invResidueMap, sysStrandNum);
		
		System.out.print("Num unpruned rot: ");
		for (int i=0; i<numUnprunedRot.length; i++)
			System.out.print(numUnprunedRot[i]+" ");
		System.out.println();
	}
	
	//Traverses the current tree starting at the root
	public void traverseTree(StrandRotamers sysLR, StrandRotamers ligRot, Molecule m, RotamerLibrary rl,
			RotamerLibrary grl, PrunedRotamers<Boolean> prunedRot, int numTotalRot, int rotIndOffset[], PairwiseEnergyMatrix eMatrix){
		
		traverseTreeHelper(root, sysLR, ligRot, m, rl, grl, prunedRot, numTotalRot, rotIndOffset, eMatrix);
	}

	//Performs post-order traversal of the tree
	private void traverseTreeHelper(TreeNode n, StrandRotamers sysLR, StrandRotamers ligRot, Molecule m, RotamerLibrary rl,
			RotamerLibrary grl, PrunedRotamers<Boolean> prunedRot, int numTotalRot, int rotIndOffset[], PairwiseEnergyMatrix eMatrix){
		
		TreeNode lc = n.getlc();
		if (lc!=null) //traverse left subtree first
			traverseTreeHelper(lc, sysLR, ligRot, m, rl, grl, prunedRot, numTotalRot, rotIndOffset, eMatrix);
		
		TreeNode rc = n.getrc();
		if (rc!=null) //then traverse right subtree
			traverseTreeHelper(rc, sysLR, ligRot, m, rl, grl, prunedRot, numTotalRot, rotIndOffset, eMatrix);
		
		if (n==root) //done
			return;
		
		else { //not at the root, so do the computation for the edge for which node n is the child
			TreeEdge curEdge = n.getCofEdge();
			if (curEdge.getIsRootEdge()){ //this is the root edge, so finish computation and output results
				
				System.out.println("Starting A matrix computation for root edge..");
			
				curEdge.compLlambda(); //compute the L and lambda sets for the root edge
				curEdge.computeA(sysLR, ligRot, m, rl, grl, prunedRot, numTotalRot, rotIndOffset, eMatrix, G);
				
				String ligType = null;
				if (ligRot!=null)
					ligType = grl.getAAName(ligRot.getIndexOfNthAllowable(0,0));
					
				curEdge.outputBestStateE(m,rl,ligType); //output the results
			}
			
			else {
				curEdge.compLlambda(); //compute the L and lambda sets for the current edge
				if (curEdge.getIsLambdaEdge()){ //lambda edge, so perform the A matrix computation
					System.out.println("Starting A matrix computation for edge ("+curEdge.getNodeName1()+","+curEdge.getNodeName2()+")..");
					System.out.println("lambda: "+curEdge.getLambda().size()+", L: "+curEdge.getL().size()+", M: "+curEdge.getM().size());
					curEdge.computeA(sysLR, ligRot, m, rl, grl, prunedRot, numTotalRot, rotIndOffset, eMatrix, G);
				}
			}
		}
	}
	
	public void computeLambdaSets(TreeNode node)
	{
		computeLambdaSets(node, true);
	}
	public void computeLambdaSets(TreeNode node, boolean initMatrices)
	{
		TreeNode lc = node.getlc();
		if (lc!=null) //traverse left subtree first
			computeLambdaSets(lc, initMatrices);
		
		TreeNode rc = node.getrc();
		if (rc!=null) //then traverse right subtree
			computeLambdaSets(rc, initMatrices);
		
		if(node.getCofEdge() != null)
			node.getCofEdge().compLlambda(initMatrices);
	}
	
	public void outputBestStateE(Molecule m, StrandRotamers ligRot, RotamerLibrary rl, RotamerLibrary grl)
	{
		String ligType = null;
		if (ligRot!=null)
			ligType = grl.getAAName(ligRot.getIndexOfNthAllowable(0,0));
		root.getCofEdge().outputBestStateE(m, rl, ligType);
	}
	
	//Read the branch decomposition from file fName
	private void readBranchDecomposition(String fName, Molecule m, int numUnprunedRot[], int molResidueMap[], int invResidueMap[], int sysStrandNum){
		
		BufferedReader bufread = null;
		try {
			File file = new File(fName);		
			FileReader fr = new FileReader(file);  
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			System.out.println("ERROR: Residue interaction graph file not found");
			System.exit(1);
		}
		
		//Read the node information
		String str = readLine(bufread,fName);
		TreeNode tn[] = new TreeNode[new Integer(getToken(str,2)).intValue()];
		readLine(bufread,fName); //skip two comment lines
		readLine(bufread,fName);
		for (int i=0; i<tn.length; i++){ //read in the information for each node
			
			str = readLine(bufread,fName); //read node line
			
			int tnName = new Integer(getToken(str,1)).intValue();
			//int tnEdges = new Integer(getToken(str,2)).intValue();
			boolean tnIsLeaf = new Boolean(getToken(str,3)).booleanValue();
			
			if (tnIsLeaf) { //leaf node
				int tnV1 = m.mapPDBresNumToMolResNum(new Integer(getToken(str,4)).intValue()); //transform to molecule-relative numbering
				int tnV2 = m.mapPDBresNumToMolResNum(new Integer(getToken(str,5)).intValue());
				
				tn[i] = new TreeNode(tnName,tnIsLeaf,tnV1,tnV2); //add node
				
				//Add the two graph vertices and the corresponding edge to the residue interaction graph G
				G.addV(tnV1);
				G.addV(tnV2);
				G.addEdge(tnV1,tnV2);
			}
			
			else //internal node
				tn[i] = new TreeNode(tnName,tnIsLeaf,-1,-1);
		}
		
		//Read the edge information
		readLine(bufread,fName); //skip blank line
		str = readLine(bufread,fName);
		TreeEdge te[] = new TreeEdge[new Integer(getToken(str,2)).intValue()];
		readLine(bufread,fName); //skip three comment lines
		readLine(bufread,fName);
		readLine(bufread,fName);
		for (int i=0; i<te.length; i++){
			
			str = readLine(bufread,fName); //read edge line
			//int teName = new Integer(getToken(str,1)).intValue();
			int teNode1 = new Integer(getToken(str,2)).intValue();
			int teNode2 = new Integer(getToken(str,3)).intValue();
			int teWidth = new Integer(getToken(str,4)).intValue();
			
			str = readLine(bufread,fName); //read the M set
			LinkedHashSet<Integer> teM = new LinkedHashSet<Integer>();
			for (int j=0; j<teWidth; j++){
				teM.add(m.mapPDBresNumToMolResNum(new Integer(getToken(str,j+1)).intValue()));
			}
			
			te[i] = new TreeEdge(teNode1, teNode2, teM, numUnprunedRot, molResidueMap, invResidueMap, sysStrandNum, false);
		}

		//Transform the branch decomposition into a rooted tree
		transformRootedTree(tn, te, numUnprunedRot, molResidueMap, invResidueMap, sysStrandNum);
		
		try { bufread.close(); } catch(Exception e){} //done, so close file for reading
	}
	
	//Called by readBranchDecomposition();
	//Transforms the branch decomposition (nodes in tn[], edges in te[]) into a rooted tree:
	//		Select an arbitrary edge and add node s between the vertices incident with this edge; add the root and connect it to node s
	private void transformRootedTree(TreeNode tn[], TreeEdge te[], int numUnprunedRot[], int molResidueMap[], int invResidueMap[], int sysStrandNum){
		
		//First, add the root and a node s, and the edge between root and s
		root = new TreeNode(-3,false,-1,-1);
		TreeNode s = new TreeNode(-2,false,-1,-1);
		TreeEdge rs = new TreeEdge(root.getName(),s.getName(),(new LinkedHashSet<Integer>()), numUnprunedRot, molResidueMap, invResidueMap, sysStrandNum, true);
		setPCvar(root,s,rs,true); //set s to be the left (and only) child of the root
		
		
		//Next, select the first edge from te[], and use it to root the tree; delete the edge from te[]
		final int useEdgeInd = 0;
		TreeEdge se = te[useEdgeInd].deepCopy();
		TreeEdge tte[] = new TreeEdge[te.length-1];
		int cnt = 0;
		for (int i=0; i<te.length; i++) {
			if (i!=useEdgeInd) {
				tte[cnt] = te[i];
				cnt++;
			}
		}
		te = tte;
		
		int lcInd = getInd(tn,se.getNodeName1()); //set the left child of s 
		setPCvar(s,tn[lcInd],new TreeEdge(s.getName(),tn[lcInd].getName(),se.getM(),numUnprunedRot, molResidueMap, invResidueMap, sysStrandNum, false),true);
		
		int rcInd = getInd(tn,se.getNodeName2()); //set the right child of s
		setPCvar(s,tn[rcInd],new TreeEdge(s.getName(),tn[rcInd].getName(),se.getM(),numUnprunedRot, molResidueMap, invResidueMap, sysStrandNum, false),false);
		
		//Finally, build the tree by following the edges in te[]
		buildRootedTreeHelper(tn[lcInd],te,tn);
		buildRootedTreeHelper(tn[rcInd],te,tn);
	}
	
	//Find the edges incident with node r, update the corresponding variables, and follow the edges to build the entire tree
	private void buildRootedTreeHelper(TreeNode r, TreeEdge te[], TreeNode tn[]){
		
		if (r.getIsLeaf()) { //current node is leaf node, so only one incident edge (leading to this node)
			return;
		}
		
		else { //current node is internal, so exactly three incident edges (one leading to this node, and two starting from this node)
			
			//get the edges starting from node r, and the corresponding left and right children of r
			int cNames[] = new int[2];
			TreeEdge ce[] = new TreeEdge[2];
			getChildrenOfNode(r,te,cNames,ce);
			
			int lcInd = getInd(tn,cNames[0]); //set the left child of r
			setPCvar(r,tn[lcInd],ce[0],true);
			
			int rcInd = getInd(tn,cNames[1]); //set the right child of r
			setPCvar(r,tn[rcInd],ce[1],false);
			
			//recursively build the tree
			buildRootedTreeHelper(tn[lcInd],te,tn);
			buildRootedTreeHelper(tn[rcInd],te,tn);
		}
	}
	
	//Find the two children edges of the internal node r and the corresponding node names of the two children;
	//		do not call this when r is a leaf node
	private void getChildrenOfNode(TreeNode r, TreeEdge te[], int cNames[], TreeEdge ce[]){
		
		int cn = 0;
		for (int i=0; i<te.length; i++){
			
			if ( te[i].getNodeName1()==r.getName() ) { //found an edge incident with node r
				
				if ((te[i].getNodeName2()!=r.getp().getName())) { //the other node incident with this edge is not the parent of node r
					
					cNames[cn] = te[i].getNodeName2();
					ce[cn] = te[i];
					cn++;
				}
			}
			
			else if ( te[i].getNodeName2()==r.getName() ) { //found an edge incident with node r
				
				if ((te[i].getNodeName1()!=r.getp().getName())) { //the other node incident with this edge is not the parent of node r
					
					cNames[cn] = te[i].getNodeName1();
					ce[cn] = te[i];
					cn++;
				}
			}
			
			if (cn==cNames.length) //found both children
				break;
		}
		
		if (cn!=2) {
			System.out.println("ERROR: could not find all nodes adjacent with node "+r.getName());
			System.exit(1);
		}
	}
	
	//Sets the necessary variables for the parent node p, child node c (left if isL==true, right otherwise), and the edge pc between p and c
	private void setPCvar(TreeNode p, TreeNode c, TreeEdge pc, boolean isL){
		
		pc.setP(p);
		pc.setC(c);
		
		c.setP(p);
		c.setCofEdge(pc);
		
		if (isL)
			p.setLc(c);		
		else
			p.setRc(c);
	}
	
	//Returns the index of the element of a[] with name s; returns -1 if s is not found in a[]
	private int getInd(TreeNode a[], int s){
		for (int i=0; i<a.length; i++){
			if (a[i].getName()==s)
				return i;
		}
		return -1;
	}
	
	//Reads the next line from bufread
	private String readLine(BufferedReader bufread, String fName){
		try {
			return bufread.readLine();
		}
		catch ( Exception e ){
			System.out.println("ERROR: An error occurred while reading input from file "+fName);
			System.exit(1);
		}
		return null;
	}
	
	// This function returns the xth token in string s
	private String getToken(String s, int x) {
	
		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");
		
		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
			  st.nextToken();
			else {
				// System.out.println("ERROR: Unable to access parameter " + x);
				return(new String(""));
			}
		}
		
		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken

	public void setEnumerationObjects(BWMAStarNode asroot,
			BWMSolutionSpace space) {
		AStarRoot = asroot;
		solutionSpace = space;
		root.setEnumerationObjects(asroot,space);
	}

	public TreeNode getRoot() {
		return root;
	}

	public InteractionGraph getGraph() {
		return G;
	}
}
