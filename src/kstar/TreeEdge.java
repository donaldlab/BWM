package kstar;

import java.util.*;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import BDAStar.BWMAStarNode;
import BDAStar.BWMSolutionSpace;
import BDAStar.Conformation;
import BDAStar.Position;
import BDAStar.ProteinConformation;
import BDAStar.ProteinPosition;
import BDAStar.SolutionSpace;

public class TreeEdge implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private int nodeName1 = -1; //the name of the first incident node
	private int nodeName2 = -1; //the name of the second incident node

	@SuppressWarnings("unused")
	private TreeNode p = null; //the parent node incident to the edge in the rooted tree
	private TreeNode c = null; //the child node incident to the edge in the rooted tree
	
	private boolean isLambdaEdge = false; //determines if this is a lambda edge
	
	private LinkedHashSet<Integer> M = null; //the M set (vertices have molecule index-relative numbering)
	private LinkedHashSet<Integer> L = null; //the L set (vertices have molecule index-relative numbering)
	private LinkedHashSet<Integer> lambda = null; //the lambda set (vertices have molecule index-relative numbering)
	
	private LinkedHashSet<TreeEdge> Fset = null; // list to contain the F set of the edges
	
	private int molResMap[] = null; //maps redesign-order indices to molecule-relative residue numbers
	private int invResMap[] = null; //maps the indices of molecule-relative residue numbers to redesign-order indexing (used in the energy matrices); ligand is mapped last
	
	private int A[][] = null; //the A matrix for the current tree edge
	private float energy[] = null; //the computed energy for each state combination in A
	
	private int sysStrNum = -1; //the strand number in the molecule of the system strand
	
	private int numStates[] = null; //the number of allowed states (unpruned rotamers) for each redesign position
	
	private RotTypeMap rtm[][] = null; //the mapping between state indices and residue/aa/rot tuples
	
	private boolean isRootEdge = false; //determines if this is the root edge
	
	/* Enumeration objects */
	private static BWMAStarNode root;
	private static BWMSolutionSpace solutionSpace;
	private PriorityQueue<Conf>[] A2;
	private Set<TreeEdge> rightFSet;
	private Map<Integer[],Integer> rightSolutionOffset;
	private List<Integer[]> rightSolutions;
	
	public TreeEdge(int eNodeName1, int eNodeName2, LinkedHashSet<Integer> teM,
			int numUnprunedRot[], int molResidueMap[], int invResidueMap[], int sysStrandNum, boolean rootEdge){
	
		
		nodeName1 = eNodeName1;
		nodeName2 = eNodeName2;
		
		M = new LinkedHashSet<Integer>(teM);
		
		numStates = numUnprunedRot;
		molResMap = molResidueMap;
		invResMap = invResidueMap;
		sysStrNum = sysStrandNum;
		
		isRootEdge = rootEdge;
	}
	
	public void setBWMAStarObjects( BWMAStarNode node, BWMSolutionSpace space)
	{
		root = node;
		solutionSpace = space;
	}
	
	//Computes the L and lambda sets for this edge; must be called only after p and c have been assigned for all edges
	public void compLlambda(){
		
		TreeNode clc = c.getlc();
		if(clc==null)//child is leaf tree node, so the L set is the difference between the M set and the two graph vertices for the child
		{
			L = new LinkedHashSet<Integer>();
			L.add(c.getv1());
			L.add(c.getv2());
			L.removeAll(M);
			
			lambda = new LinkedHashSet<Integer>(L); // as the lambda set and L set for a leaf node would be the same
		}
		
		else
		{
			LinkedHashSet<Integer> uMc = new LinkedHashSet<Integer>(clc.getCofEdge().getM());
			uMc.addAll(c.getrc().getCofEdge().getM()); //the union of the M sets for the two incident edges with the two children
			lambda = new LinkedHashSet<Integer>(uMc);
			lambda.removeAll(M); //the difference between the M set for this edge and the uMc set is equal to the lambda set
			
			L = new LinkedHashSet<Integer>(lambda);
			L.addAll(clc.getCofEdge().getL()); //add the L set of the left edge
			L.addAll(c.getrc().getCofEdge().getL()); //add the L set of the right edge
		}
		
		if(!lambda.isEmpty())
		{
			isLambdaEdge=true; // initialising the matrices and calculating the Fset since it is a lambda edge
			initializeMatrices();
			rtm = new RotTypeMap[M.size()+lambda.size()][];
			Fset = new LinkedHashSet<TreeEdge>();
			computeFset(c);
		}
		else
			isLambdaEdge=false;
	}
	
	//Initialize the A[] and energy[] matrices;
	//		Each entry in the first dimension of A corresponds to a unique state assignment in M; for each such entry,
	//			the second dimension gives the best state for each graph vertex in lambda
	private void initializeMatrices(){
		
		int size = 1;
		for (int i=0; i<numStates.length; i++){
			if (M.contains(molResMap[i])){
				size *= numStates[i];
			}
		}
		A = new int[size][lambda.size()];
		energy = new float[size];
		
		for (int i=0; i<size; i++){
			energy[i] = Integer.MAX_VALUE;
			for (int j=0; j<A[i].length; j++){
				A[i][j] = -1;
			}
		}
	}
	
	public LinkedHashSet<Integer> getL(){
		return L;
	}
	
	public LinkedHashSet<Integer> getM(){
		return M;
	}
	
	public boolean getIsLambdaEdge(){
		return isLambdaEdge;
	}
	
	public LinkedHashSet<Integer> getLambda(){
		return lambda;
	}
	
	public int [][] getA(){
		return A;
	}
	
	public TreeNode getc(){
		return c;
	}
	
	public float [] getEnergy(){
		return energy;
	}
	
	//Computes and stores the A matrix for the current edge; must be called after the L and lambda sets for the current edge have already been computed (using compLlambda())
	public void computeA(StrandRotamers sysLR, StrandRotamers ligRot, Molecule m, RotamerLibrary rl, RotamerLibrary grl, 
			PrunedRotamers<Boolean> prunedRot, int numTotalRot, int rotIndOffset[], PairwiseEnergyMatrix eMatrix, InteractionGraph G){
		
		int maxDepth = M.size() + lambda.size();
		
		Object arrayM[] = M.toArray();
		Object arrayLambda[] = lambda.toArray();
		
		int curState[] = new int[maxDepth];
		int bestState[] = new int[maxDepth]; //This becomes a heap if we port directly.
		for (int i=0; i<maxDepth; i++){
			curState[i] = -1;
			bestState[i] = -1;
		}
		
		float bestEnergy[] = new float[]{Float.MAX_VALUE};
		float energy_store[] = new float[]{Float.MAX_VALUE};
		
		computeAhelper(0, maxDepth, arrayM, arrayLambda, sysLR, ligRot, m, rl, grl, prunedRot, numTotalRot, rotIndOffset,
				curState, eMatrix, G, bestState, bestEnergy,energy_store);
	}
	
	//Called by computeA(.)
	private void computeAhelper(int depth, int maxDepth, Object arrayM[], Object arrayLambda[], StrandRotamers sysLR, StrandRotamers ligRot,
			Molecule m, RotamerLibrary rl, RotamerLibrary grl, PrunedRotamers<Boolean> prunedRot, int numTotalRot, int rotIndOffset[], 
			int curState[], PairwiseEnergyMatrix eMatrix, InteractionGraph G, int bestState[], float bestEnergy[], float energy_store[]){		
		if (depth >= maxDepth){ //end level of recursive calls; call the backtracking procedure to look-up the optimal states for (L-lambda)
			
			
			float energy_ll=0; //energy got from the Fset edges - energy corresponding to L- lambda set 
			energy_ll=bTrack(curState);
			
			float en[] = new float[2];
			en[0]=0;
			en[1]=0;
			
			computeEforState(curState,eMatrix,m,G,en); //en[0] will contain the energy to compare, en[1] will contain the energy to store
			float total_energy=0;
			total_energy=en[0]+energy_ll;
                      
			PriorityQueue<Conf> conformationHeap = null;// A2[computeIndexInA(curState)];
			conformationHeap.add(new Conf(curState, en[0], computeIndexInA(curState)));
			
			if ( (total_energy<bestEnergy[0]) || (bestEnergy[0]==Float.MAX_VALUE) ) { //new best energy, so update to the current state assignment
				
				bestEnergy[0] = total_energy;
				energy_store[0]=energy_ll+en[1];
				System.arraycopy(curState, 0, bestState, 0, curState.length);
			}
		}
		
		else { //setup current level: Do not need to touch this part 
			
			Object vArray[] = null;
			int vInd = -1;
			if (depth < M.size()){ //work on the M set
				vArray = arrayM;
				vInd = depth;
			}
			else { //work on the lambda set
				vArray = arrayLambda;
				vInd = depth - M.size();
			}
			
			rtm[depth] = new RotTypeMap[numStates[invResMap[(Integer)vArray[vInd]]]];
			
			int curPos = invResMap[(Integer)vArray[vInd]];
			
			int curStrandResNum = m.residue[(Integer)vArray[vInd]].strandResidueNumber;
			
			StrandRotamers str = null;
			RotamerLibrary rotLib = null;
			if (m.residue[(Integer)vArray[vInd]].strandNumber==sysStrNum){ //this residue is in the system strand
				str = sysLR;
				rotLib = rl;
			}
			else { //this residue is in the ligand strand
				str = ligRot;
				rotLib = grl;
			}
			
			for(int q=0;q<str.getNumAllowable(curStrandResNum);q++) { //for all allowed amino acid types
				
				int AAindex = str.getIndexOfNthAllowable(curStrandResNum,q);
				
				int numRot = rotLib.getNumRotForAAtype(AAindex);
				if (numRot==0)
					numRot = 1;
				
				for(int w=0;w<numRot;w++) { //for all rotamers
					
					int rotInd = -1;
					if (m.residue[(Integer)vArray[vInd]].strandNumber==sysStrNum){ //this residue is in the system strand
						rotInd = curPos*numTotalRot + rotIndOffset[AAindex] + w;
					}
					else { //this residue is in the ligand strand
						rotInd = curPos*numTotalRot + w;
					}
					
					if (!prunedRot.get(curPos, AAindex, w)){ //rotamer not pruned, so check
						
						curState[depth]++;
						
						if (rtm[depth][curState[depth]]==null) //store the mapping from state index to residue/aa/rot tuple
							rtm[depth][curState[depth]] = new RotTypeMap(curPos,AAindex,w);
						
						computeAhelper(depth+1, maxDepth, arrayM, arrayLambda, sysLR, ligRot, m, rl, grl, prunedRot, numTotalRot, rotIndOffset,
								curState, eMatrix, G, bestState, bestEnergy,energy_store);
					}	
				}
			}
			
			curState[depth] = -1;
			
			if ( depth==M.size() ){//done with the lambda states for the current state assignment in M, so update A[] and energy[]
			    /** The self-balancing happens automatically here so long as we work bottom-up, which is the case! */
				if(isRootEdge) // as you will never look up the A matrix for the root, i.e it will never be in the F set of any edge, we can store the actual energy
				{
					storeBestStateLambda(bestState, arrayM, bestEnergy[0]); //store the best state for each vertex in lambda, for the current state assignment in M
				}
				else
				{
					storeBestStateLambda(bestState, arrayM, energy_store[0]); //store the best state for each vertex in lambda, for the current state assignment in M
				}
	
				bestEnergy[0] = Float.MAX_VALUE;
				energy_store[0] = Float.MAX_VALUE;
				
				for (int i=0; i<bestState.length; i++)
					bestState[i] = -1;
				
			}
		}
	}
	
	
	private float bTrack(int curState[]){		
		
		
		Object arrayF[] = Fset.toArray();
		float energy_return=0;
		
		//Find the state assignments for the edges in the M sets of each of the edges in the set F
		int fkMstate[][] = new int[arrayF.length][];
		for (int i=0; i<fkMstate.length; i++)
			fkMstate[i] = getMstateForEdgeCurState(curState,(TreeEdge)arrayF[i]);
		
		 // code to sum the energies from the F set and return that value
		for(int i=0; i<arrayF.length;i++)

		{
			TreeEdge fk = (TreeEdge)arrayF[i];
			int index = fk.computeIndexInA(fkMstate[i]);
			energy_return+=fk.getEnergy()[index];
		}
		
		/* Handle right side... 
		TreeEdge[] arrayRightF = rightFSet.toArray(new TreeEdge[]{});
		int[][] rightMState = new int[arrayRightF.length][];
                for (int i=0; i<rightMState.length; i++)
                    rightMState[i] = getMstateForEdgeCurState(curState,(TreeEdge)arrayF[i]);
		for(int j = 0; j < arrayRightF.length; j++)
		{
		    int index = arrayRightF[j].computeIndexInA(rightMState[j]);
		    energy_return += arrayRightF[j].getEnergy()[index];
		}
		*/
		return energy_return;
	}
	
	//Called by bTrackHelper(.); finds the tree edges belonging to the set F starting at the sub-tree rooted at the tree node tn
	private void computeFset(TreeNode tn){
		
		TreeNode clc = tn.getlc();
		TreeNode crc = tn.getrc();
		
		if (clc!=null){
			TreeEdge clce = clc.getCofEdge();
			if (clce.getIsLambdaEdge()) //if clce is a lambda edge, add it to Fi, and do not travers the subtree rooted at clc
				Fset.add(clce);
			else //not a lambda edge, so traverse the subtree rooted at clc
				computeFset(clc);				
		}
		
		if (crc!=null){
			TreeEdge crce = crc.getCofEdge();
			if (crce.getIsLambdaEdge()) 
			{
				Fset.add(crce);
				rightFSet.add(crce);
			}
			else
				computeFset(crc);			
		}
	}
	
	//Computes the state of all vertices in e.M (returned in eMstate[]) corresponding to the state assignment curState[] for this edge;
	//NOTE: This must be called only when edge e is in the F set of this edge
	public int [] getMstateForEdgeCurState(int curState[],TreeEdge e){
		
		int eMstate[] = new int[e.getM().size()];
		
		Object curM[] = M.toArray();
		Object curLambda[] = lambda.toArray();
		Object MLambda[] = new Object[curM.length+curLambda.length];
		System.arraycopy(curM, 0, MLambda, 0, curM.length);
		System.arraycopy(curLambda,0,MLambda,curM.length,curLambda.length);
		curM = null; curLambda = null;
		
		Object eM[] = e.getM().toArray();
		
		for (int i=0; i<eM.length; i++){ //find the state for each vertex in e.M (it must be either in this.M or this.lambda)
			
			eMstate[i] = -1;
			
			for (int j=0; j<MLambda.length; j++){
				
				int eMi = ((Integer)eM[i]).intValue();
				int MLambdaj = ((Integer)MLambda[j]).intValue();
				
				if (eMi==MLambdaj){
					
					int p = rtm[j][curState[j]].pos;
					int a = rtm[j][curState[j]].aa;
					int r = rtm[j][curState[j]].rot;
					
					RotTypeMap ertm[] = e.getrtm()[i];
					for (int k=0; k<ertm.length; k++){
						if ( (ertm[k].pos==p) && (ertm[k].aa==a) && (ertm[k].rot==r) ){ //found state for vertex eM[i]
							eMstate[i] = k;
							break;
						}
					}
					break;
				}
			}
		}
		
		return eMstate;
	}
	
	
	
	
	private void computeEforState(int curState[],PairwiseEnergyMatrix eMatrix,Molecule m, InteractionGraph G, float en[]){
		
		int pi=0,ai=0,ri=0,pj=0,aj=0,rj=0;
		en[0] =eMatrix.getShellShellE(); //Add shell shell energy
		int numPos = M.size() + lambda.size();
		for(int i=0;i<numPos;i++){
			
			pi=rtm[i][curState[i]].pos;
			ai=rtm[i][curState[i]].aa;
			ri=rtm[i][curState[i]].rot;
			if(i>=M.size()){
				/*
				en[1]+=eMatrix[pi][ai][ri][pi][0][1]; //add the self energy of the rotamer of the lambda residue
				en[1]+=eMatrix[pi][ai][ri][pi][0][0];
				*/
				en[1]+=eMatrix.getShellRotE(pi, ai, ri);
			}
			/*
			en[0]+=eMatrix[pi][ai][ri][pi][0][1]; //add the self energy of the rotamer of the lambda residue
			en[0]+=eMatrix[pi][ai][ri][pi][0][0];
			*/
			en[0]+=eMatrix.getShellRotE(pi, ai, ri);

			for(int j=i+1;j<numPos;j++){
			
				pj=rtm[j][curState[j]].pos;
				aj=rtm[j][curState[j]].aa;
				rj=rtm[j][curState[j]].rot;

				if(G.edgeExists(molResMap[pi],molResMap[pj])){  // if edge exists between residues in the interaction graph

					if(j<M.size()) //intreacting between M set
						//en[0]+=eMatrix[pi][ai][ri][pj][aj][rj];
						en[0]+=eMatrix.getPairwiseE(pi, ai, ri, pj, aj, rj);
					else if(j>=M.size()){ //interaction between lambda set or M and lambda set
						//en[0]+=eMatrix[pi][ai][ri][pj][aj][rj];
						en[0]+=eMatrix.getPairwiseE(pi, ai, ri, pj, aj, rj);
						//en[1]+=eMatrix[pi][ai][ri][pj][aj][rj];
						en[1]+=eMatrix.getPairwiseE(pi, ai, ri, pj, aj, rj);
					}
				}
			}
			
		}	

	}
	
	//Compute the best state for each vertex in lambda, for the given state assignment in M
	private void storeBestStateLambda(int curState[], Object arrayM[], float curEnergy) {
		
		int curIndInA = computeIndexInA(curState); //get the index corresponding to the current state assignment in M
		
		for (int i=0; i<lambda.size(); i++){
			A[curIndInA][i] = curState[arrayM.length+i];
		}
		
		energy[curIndInA] = curEnergy;
		
	}
	
	//Compute the index into the A matrix for this tree edge, given the state assignments for the vertices in M in curState[]
	public int computeIndexInA(int curState[]){
		
		if (isRootEdge) //the M set is empty for the root edge
			return 0;
		
		Object array_M[] = M.toArray();
		int index = 0;
		int s = 1;
		
		index += curState[array_M.length-1];
		for (int i=(array_M.length-2); i>=0; i--){ //find the state assignment for the vertices in M
			
			if (curState[i]<0){
				System.out.println("ERROR: GD state not fully assigned");
				System.exit(1);
			}
			
			s *= numStates[invResMap[(Integer)array_M[i+1]]];
			index += (curState[i] * s);
		}
		
		return index;
	}
	
	public RotTypeMap[][] getrtm(){
		return rtm;
	}
	
	public void setP(TreeNode pn){
		p = pn;
	}
	
	public void setC(TreeNode cn){
		c = cn;
	}
	
	public int getNodeName1(){
		return nodeName1;
	}
	
	public int getNodeName2(){
		return nodeName2;
	}
	
	public boolean getIsRootEdge(){
		return isRootEdge;
	}
	
	//Finds the best energy and outputs the corresponding state, using the information for this edge;
	//	This must be called only for the root edge, since otherwise not all elements of aaRotPos[bestInd][][] are defined
	//Modify this function and add what to do for rootedge in each and every function
	/** This is the method you need to tweak to make things work.
	 * To wit:
	 * 1. Poll 
	 * 2. Rebalance
	 * 3. Process rightSide?
	 * */
	public void outputBestStateE(Molecule m, RotamerLibrary rl, String ligType){
		
		if (!isRootEdge) {
			System.out.println("ERROR: cannot output best state for non-root edges");
			System.exit(1);
		}
	
		// the energy matrix and the A matrix will only have one entry 	
		RotTypeMap bestPosAARot[] = new RotTypeMap[molResMap.length]; // creating a variable to store the best energy returned by the Btrack Procedure, molresMap
										// length is equal to the number of residues being designed
		/** TODO: Backtrack entry point. */
		bTrackBestConf(bestPosAARot,A[0]);
		System.out.print("GMEC: ");
		
		for (int i=0; i<bestPosAARot.length; i++) { //output the AAs
			
			if (m.residue[molResMap[i]].strandNumber==sysStrNum) //this residue is in the system strand
				System.out.print(rl.getAAName(bestPosAARot[i].aa)+" ");
			
			else //this residue is in the ligand strand
				System.out.print(ligType+" ");
		}
		
		for (int i=0; i<bestPosAARot.length; i++) //output the rotamers
			System.out.print(bestPosAARot[i].rot+" ");
		
		System.out.println(energy[0]); //output the energy
	}
	
	//Returns a deep copy of this TreeEdge (modified from http://javatechniques.com/blog/faster-deep-copies-of-java-objects/, accessed 10/30/2008)
	public TreeEdge deepCopy(){
		TreeEdge c = null;
		try {
			ByteArrayOutputStream b = new ByteArrayOutputStream();
			ObjectOutputStream fout = new ObjectOutputStream(b);
			fout.writeObject(this);
			fout.flush();
			fout.close();
			
			ObjectInputStream fin = new ObjectInputStream(new ByteArrayInputStream(b.toByteArray()));
			c = (TreeEdge)fin.readObject();
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}
		return c;
	}
	/**
	 * This method recursively populates bestPosAARot[] with the values from bestState[]
	 * @param bestPosAARot The full rotamer assignment array, I presume.
	 * @param bestState The assignment for the MSet of the current node.
	 */
	public void bTrackBestConf (RotTypeMap bestPosAARot[], int bestState[])
	{
		int position=-1;
		int index=0;
		if(!isRootEdge)
			index=M.size();
		/* populate the outcome with this edge's assignments */
		for(int i=index; i<index+lambda.size();i++)
		{
			position=rtm[i][bestState[i]].pos;
			bestPosAARot[position]= new RotTypeMap(rtm[i][bestState[i]].pos,rtm[i][bestState[i]].aa,rtm[i][bestState[i]].rot);
		}
		/* If we have no further lambda children, return. */
		if(Fset == null)
			return;
		/* Calculate the M Set for children. 
		 * No changes here. */
		Object array_fset[] = Fset.toArray();
		int fkMstate[][] = new int[array_fset.length][];
		for (int i=0; i<fkMstate.length; i++)
			fkMstate[i] = getMstateForEdgeCurState(bestState,(TreeEdge)array_fset[i]);
		
		/* Index into the child array 
		 * No changes here?*/
		int fkLambdaState[][] = new int[array_fset.length][];
		for (int i=0; i<fkMstate.length; i++){
			Object fkM[] = ((TreeEdge)array_fset[i]).getM().toArray();
			fkLambdaState[i] = ((TreeEdge)array_fset[i]).getA()[((TreeEdge)array_fset[i]).computeIndexInA(fkMstate[i])];
		}
		
		//For each edge in the set F, combine the state assignments for the M and lambda sets into a single array
		/* Store M + lambda as "best state" for each child 
		 * No changes here */
		int fkMLambda[][] = new int[array_fset.length][];
		for (int i=0; i<fkMstate.length; i++){
			fkMLambda[i] = new int[fkMstate[i].length+fkLambdaState[i].length];
			System.arraycopy(fkMstate[i], 0, fkMLambda[i], 0, fkMstate[i].length);
			System.arraycopy(fkLambdaState[i], 0, fkMLambda[i], fkMstate[i].length, fkLambdaState[i].length);
		}
		
		//call the recursive procedure 
		/* Recurse
		 * Needs to:
		 * 1. Recurse on left children
		 * 2. Manage right children solution list. 
		 * */
		for(int i=0;i<array_fset.length;i++)
		{
			TreeEdge fk = (TreeEdge)array_fset[i];
			fk.bTrackBestConf(bestPosAARot,fkMLambda[i]);
		}
	
		return;
	}
	
	/* Update energies method */
	private void updateEnergies()
	{
	    
	}
	
	public double nextBestEnergy(RotTypeMap[] partialConformation)
	{
	    return 0;
	}

	public void setPositions(LinkedHashSet<Integer> set)
	{
	    lambda = set;
	    if(lambda.size() > 0)
	        isLambdaEdge = true;
	}

	public Set<Position> getPositionSet () {
	    LinkedHashSet<Position> out = new LinkedHashSet<Position>();
	    if(lambda == null)
	    {
	        for(int i = 0; i < 2; i ++)
	            lambda.add(i);
	    }

	    for(Integer i : lambda)
	    {
	        out.add(new Position(i));
	    }
	    return out;
	}

	public List<? extends Position> getPositionList () {
	    LinkedList<Position> out = new LinkedList<Position>();
	    for(Integer i : lambda)
	    {
	        int store = i;
	        if(solutionSpace != null)
	            out.add(solutionSpace.positionFromPos(invResMap[i]));
	        else out.add(new Position(store));
	    }
	    return out;
	}
	
	/**
	 * Algorithm:
	 * 1. RootEdge gets node at the top of heap
	 * 2. Pass root's partial conformation to children
	 * 3. Children use TreeNode to get the next heap, to extend the partial conformation
	 * 4. Children pass the extended partial conformation to their children.
	 * 5. Poll the next best energy, and modify the conformation.
	 * 6. From the bottom up, children reinsert conformations into their heaps using the next best subconformation energy.
	 * @author Jon
	 *
	 */
	public RotTypeMap[] ConftoRotTypeMap(Conf c)
	{
	    int[] bestState = c.conformation;
	    RotTypeMap[] bestPosAARot = new RotTypeMap[bestState.length];
	    for(int i=0; i<bestState.length;i++)
	    {
	        int position=rtm[i][bestState[i]].pos;
	        bestPosAARot[position]= new RotTypeMap(rtm[i][bestState[i]].pos,rtm[i][bestState[i]].aa,rtm[i][bestState[i]].rot);
	    }
	    return bestPosAARot;
	}
	
	private class Conf
	{
	    int[] conformation;
	    double energy;
	    int indexInA;
	    public Conf(int[] c, double e, int i)
	    {
	        conformation = c;
	        energy = e;
	        indexInA = i;
	    }
	    
	    public int hashCode()
	    {
	        return indexInA;
	    }
	}
	
	public int[] getInvResMap()
	{
		return invResMap;
	}
	
	private class ConformationComparator implements Comparator<Conf>
	{

            @Override
            public int compare (Conf arg0, Conf arg1) {
                if(arg0.energy - arg1.energy < 0) 
                    return -1;
                if(arg0.energy - arg1.energy > 0)
                    return 1;
                return 0;
            }
	    
	}

	public void setM(LinkedHashSet<Integer> lambda2) 
	{
		M = lambda2;
	}

}
