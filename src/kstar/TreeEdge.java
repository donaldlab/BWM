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
    private LinkedHashSet<Integer> leftL = null;
    private LinkedHashSet<Integer> rightL = null;

    private LinkedHashSet<Integer> leftOnlyL = null;

    private LinkedHashSet<TreeEdge> Fset = null; // list to contain the F set of the edges

    public int molResMap[] = null; //maps redesign-order indices to molecule-relative residue numbers
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
    private ArrayList<PriorityQueue<Conf>> A2;
    private Map<String,Integer> rightSolutionOffset;
    private Map<String, List<RightConf>> rightSolutionMap;
    Map<String, PriorityQueue<Conf>> leftHeapMap;
    private double firstRightEnergy;

    static boolean printHeap = false;
    HashMap<String, Integer> confsEnumerated = new HashMap<String, Integer>();

    TreeNode leftChild;
    TreeNode rightChild;
    TreeEdge parent;
    private float shellShellEnergy;
    private Map<String, LazyHeap<Conf>> secondaryHeapMap;
    private static ConformationComparator comparator;

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
        if(comparator == null)
            comparator = new ConformationComparator();
        rightSolutionOffset = new HashMap<String, Integer>();
        rightSolutionMap = new HashMap<String, List<RightConf>>();
        secondaryHeapMap = new HashMap<String, LazyHeap<Conf>>();
        leftHeapMap = new HashMap<String, PriorityQueue<Conf>>();
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
            leftL = L;
            leftOnlyL = leftL;
            rightL = L;
            lambda = new LinkedHashSet<Integer>(L); // as the lambda set and L set for a leaf node would be the same
        }

        else
        {
            LinkedHashSet<Integer> uMc = new LinkedHashSet<Integer>(clc.getCofEdge().getM());
            uMc.addAll(c.getrc().getCofEdge().getM()); //the union of the M sets for the two incident edges with the two children
            lambda = new LinkedHashSet<Integer>(uMc);
            lambda.removeAll(M); //the difference between the M set for this edge and the uMc set is equal to the lambda set

            L = new LinkedHashSet<Integer>(lambda);
            leftL = new LinkedHashSet<Integer>();
            leftOnlyL = new LinkedHashSet<Integer>(lambda);
            rightL = new LinkedHashSet<Integer>();
            leftL.addAll(clc.getCofEdge().getleftOnlyL());
            rightL.addAll(c.getrc().getCofEdge().getL());
            L.addAll(clc.getCofEdge().getL()); //add the L set of the left edge
            L.addAll(c.getrc().getCofEdge().getL()); //add the L set of the right edge
        }

        Fset = new LinkedHashSet<TreeEdge>();
        computeFset(c);
        initializeMatrices();
        
        if(!lambda.isEmpty())
        {
            isLambdaEdge=true; // initialising the matrices and calculating the Fset since it is a lambda edge
            
            rtm = new RotTypeMap[M.size()+lambda.size()][];

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
        A2 = new ArrayList<PriorityQueue<Conf>>(size);
        while(A2.size() < size)
        {
            A2.add(new PriorityQueue<Conf>(1, comparator));
        }
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

    public LinkedHashSet<Integer> getleftOnlyL(){
        return leftOnlyL;
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

            if(A2.size() <= computeIndexInA(curState))
            {
                A2.add(new PriorityQueue<Conf>());
            }

            PriorityQueue<Conf> conformationHeap = A2.get(computeIndexInA(curState));
            Conf newConf = new Conf(curState.clone(), 0, rtm);
            newConf.selfEnergy = en[1];
            newConf.energy = en[1] + energy_ll;
            if(leftChild != null)
                newConf.leftEnergy = energy_ll;
            //if(conformationHeap.size() < 2) 
            conformationHeap.add(newConf);

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

            if(rtm[depth] == null)
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

        return energy_return;
    }

    private float bTrackLeft(int curState[]){	
        /* TODO: Handle empty lambda set */
        if(leftChild == null) return 0;
        TreeEdge leftEdge = leftChild.getCofEdge();
        while (leftEdge.lambda.size() < 1 && leftEdge.leftChild != null)
            leftEdge = leftEdge.leftChild.getCofEdge();
        if(leftEdge.leftChild == null && leftEdge.lambda.size() < 1)
        {
            System.err.println("The impossible has happened.");
            return 0;

        }
        int[] MLambda = getMstateForEdgeCurState(curState, leftEdge);
        return (float)leftEdge.A2.get(leftEdge.computeIndexInA(MLambda)).peek().energy;
    }


    private float bTrackRight(int curState[]){		
        /* TODO: Handle empty lambda set */
        if(rightChild == null) return 0;
        TreeEdge rightEdge = rightChild.getCofEdge();
        int[] rightMLambda = getMstateForEdgeCurState(curState, rightEdge);

        return (float)rightEdge.A2.get(rightEdge.computeIndexInA(rightMLambda)).peek().energy;
    }



    //Called by bTrackHelper(.); finds the tree edges belonging to the set F starting at the sub-tree rooted at the tree node tn
    public void computeFset(TreeNode tn){

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
            }
            else
                computeFset(crc);			
        }
    }

    public void compactTree()
    {
        leftChild = searchSubtree(c.getlc());
        rightChild = searchSubtree(c.getrc());
        if(leftChild != null)
        {
            leftChild.getCofEdge().compactTree();
            leftChild.getCofEdge().parent = this;
        }
        if(rightChild != null)
        {
            rightChild.getCofEdge().compactTree();
            rightChild.getCofEdge().parent = this;
        }
        if(leftChild == null && rightChild != null)
        {
            leftChild = rightChild;
            rightChild = null;
        }
    }

    public void printTree(String prefix)
    {
        String out = prefix + "[";
        for(int i : getLambda())
            out+=invResMap[i]+", ";
        out+="]";

        out += " - L Set:[";

        for(int i : getL())
            out+=invResMap[i]+", ";
        out+="]";

        out += " - M Set:[";

        for(int i : getM())
            out+=invResMap[i]+", ";
        out+="]";
        boolean showHeaps = false;
        if(showHeaps)
        {
            out+="heaps: \n";
            for(PriorityQueue<Conf> heap : A2)
            {
                if(heap.size() > 0)
                    out+=""+heap+"\n";
            }
        }
        System.out.println(out);
        if(leftChild != null)
            leftChild.getCofEdge().printTree(prefix+"+L--");
        if(rightChild != null)
            rightChild.getCofEdge().printTree(prefix+"+R--");
    }

    public void printTreeMol(String prefix)
    {
        String out = prefix + "[";
        for(int i : getLambda())
            out+=i+", ";
        out+="]";

        out += " - L Set:[";

        for(int i : getL())
            out+=i+", ";
        out+="]";

        out += " - M Set:[";

        for(int i : getM())
            out+=i+", ";
        out+="]";
        boolean showHeaps = false;
        if(showHeaps)
        {
            out+="heaps: \n";
            for(PriorityQueue<Conf> heap : A2)
            {
                if(heap.size() > 0)
                    out+=""+heap+"\n";
            }
        }
        System.out.println(out);
        if(leftChild != null)
            leftChild.getCofEdge().printTree(prefix+"+L--");
        if(rightChild != null)
            rightChild.getCofEdge().printTree(prefix+"+R--");
    }

    private TreeNode searchSubtree(TreeNode tn){
        if(tn == null) return null;
        if(tn.getCofEdge().isLambdaEdge)
            return tn;
        if(tn.getIsLeaf())
            return null;
        TreeNode clc = searchSubtree(tn.getlc());
        TreeNode crc = searchSubtree(tn.getrc());
        if(clc != null && crc != null)
            return tn;
        if(clc == null && crc != null)
            return crc;
        if(clc != null && crc == null)
            return clc;
        return null;
    }

    //Computes the state of all vertices in e.M (returned in eMstate[]) corresponding to the state assignment curState[] for this edge;
    //NOTE: This must be called only when edge e is in the F set of this edge
    public int [] getMstateForEdgeCurState(int curState[],TreeEdge e){
        if(e.getrtm() == null)
        {
            System.out.println("Whoa!");
            System.out.println(e.getc());
            return curState; 
        }

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
        shellShellEnergy = en[0];
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

        if (isRootEdge || curState == null) //the M set is empty for the root edge
            return 0;
        
        if(M.size() < 1)
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

    public String outputBestStateE2(Molecule m, RotamerLibrary rl, String ligType){

        if (!isRootEdge) {
            System.out.println("ERROR: cannot output best state for non-root edges");
            System.exit(1);
        }

        // the energy matrix and the A matrix will only have one entry  
        RotTypeMap bestPosAARot[] = new RotTypeMap[molResMap.length]; // creating a variable to store the best energy returned by the Btrack Procedure, molresMap
        // length is equal to the number of residues being designed
        int[] bestState = new int[]{};
        double resultEnergy = getHeap(bestPosAARot, bestState).peek().energy;
        if(printHeap)
            System.out.println("Peek heap: "+resultEnergy);
        resultEnergy += shellShellEnergy;
        bTrackBestConfRemoveEarlyNew(bestPosAARot, new int[]{});
        System.out.println("RTM result: "+RTMToString(bestPosAARot));
        System.out.print("GMEC: ");

        for (int i=0; i<bestPosAARot.length; i++) { //output the AAs
            if (m.residue[molResMap[i]].strandNumber==sysStrNum) //this residue is in the system strand
                System.out.print(molResMap[i]+": "+rl.getAAName(bestPosAARot[i].aa)+" ");

            else //this residue is in the ligand strand
                System.out.print(ligType+" ");
        }

        for (int i=0; i<bestPosAARot.length; i++){ //output the rotamers

            if(bestPosAARot[i] == null) continue;
            System.out.print(bestPosAARot[i].rot+" ");
        }

        System.out.println(resultEnergy); //output the energy
        return RTMToString(bestPosAARot);
    }

    public double nextBestEnergy()
    {
        RotTypeMap bestPosAARot[] = new RotTypeMap[molResMap.length]; // creating a variable to store the best energy returned by the Btrack Procedure, molresMap
        // length is equal to the number of residues being designed
        int[] bestState = new int[]{};
        double resultEnergy = getHeap(bestPosAARot, bestState).peek().energy;
        return resultEnergy;
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

        for(int i=0;i<array_fset.length;i++)
        {
            TreeEdge fk = (TreeEdge)array_fset[i];
            fk.bTrackBestConf(bestPosAARot,fkMLambda[i]);
        }

        return;
    }

    /**
     * New Algorithm should:
     * 
     * 1. Take in current M conformation
     * 2. Populate result with best lambda conformation as a RotTypeMap[]
     * 3. Recurse, passing in a new RotTypeMap for the next best conformation given M + lambda
     * 4. Reinsert the lambda conformation into the A2 heap
     * 5. Populate RotTypeMap[] nextBestConf with what's on top of the heap now. 
     * 
     * @return
     */
    public void bTrackBestConfRemoveEarlyNew(RotTypeMap bestPosAARot[], int[] bestState)
    {
        boolean reinsert = false;
        PriorityQueue<Conf> outHeap = getHeap(bestPosAARot, bestState);
        
        
        if(printHeap)
        {
            String curString = RTMToString(bestPosAARot);
            System.out.println("====================== start "+curString+"========================");
            PriorityQueue<Conf> copy = new PriorityQueue<Conf>();

            for(Conf c : outHeap)
            {
                copy.add(c);
            }

            while(!copy.isEmpty())
            {
                Conf c = copy.poll();
                System.out.println(c+"$"+c.energy+", left "+c.leftEnergy+", self "+c.selfEnergy+", right "+c.rightEnergy);
            }
        }
        
        Conf nextState = outHeap.poll();
        /* Add lambda to the solution */
        nextState.fillRotTypeMap(bestPosAARot);
        RotTypeMap[] bestPosAARotOld = Arrays.copyOf(bestPosAARot, bestPosAARot.length);
        
        double curBest = nextState.energy;


        if(rightChild != null)
        {
            TreeEdge leftEdge = leftChild.getCofEdge();
            int[] leftMLambda = nextState.conformation; 
            int[] leftM = getMstateForEdgeCurState(leftMLambda, leftEdge);
            TreeEdge rightEdge = rightChild.getCofEdge();
            
            List<RightConf> rightConfs = getRightSolutions(bestPosAARot);
            if(rightConfs.size() < 1)
            {
                int[] rightMLambda = nextState.conformation; 
                int[] rightM = getMstateForEdgeCurState(rightMLambda, rightEdge);
                getNewRightconf(bestPosAARotOld, rightEdge, rightConfs, rightM);
            }
            
            
            LazyHeap<Conf> secondaryHeap = getSecondaryHeap(bestPosAARotOld, leftM);
            if((secondaryHeap.dirty || secondaryHeap.size() < 1) && leftEdge.moreConformations(bestPosAARotOld, leftM)) 
            {
                /* Acquire new left conformation for the heap */
                RotTypeMap[] bestPosAARotLeft = Arrays.copyOf(bestPosAARotOld, bestPosAARot.length);
                
                double newLeftEnergy = leftEdge.peekEnergy(bestPosAARotLeft, leftM);
                leftEdge.bTrackBestConfRemoveEarlyNew(bestPosAARotLeft, leftM);
                Conf newLeftConf = new Conf(bestPosAARotLeft, newLeftEnergy);
                newLeftConf.updateLeftEnergy(rightConfs.get(0).energy);
                secondaryHeap.cleanNode = newLeftConf;
                secondaryHeap.add(newLeftConf);
                secondaryHeap.dirty = false;
            }
            
            /* Maintain cleanliness */
            Conf leftConf = secondaryHeap.poll();
            leftConf.fillConf(bestPosAARot);
            if(leftConf == secondaryHeap.cleanNode)
                secondaryHeap.dirty = true;
            int index = getRightOffset(leftConf.toString());
           

            int[] rightMLambda = nextState.conformation; 
            int[] rightM = getMstateForEdgeCurState(rightMLambda, rightEdge);
            while(rightConfs.size() <= index + 1 && rightEdge.moreConformations(bestPosAARotOld, rightM))
            {
                
                getNewRightconf(bestPosAARotOld, rightEdge, rightConfs, rightM);
            }
                        
            RightConf rightConf = rightConfs.get(index);
            rightConf.fillRotTypeMap(bestPosAARot);
            
            if(index + 1 < rightConfs.size() || rightEdge.moreConformations(bestPosAARotOld, rightM))
            {
                RightConf nextRightConf = rightConfs.get(index + 1);
                leftConf.updateLeftEnergy(nextRightConf.energy);
                rightSolutionOffset.put(leftConf.toString(), index + 1);
                secondaryHeap.add(leftConf);
                reinsert = true;
            }

            if(secondaryHeap.size() > 0)
            {
                nextState.updateLeftEnergy(secondaryHeap.peek().energy);
                reinsert = true;
            }
            else
                reinsert = false;
            
        }

        else if(leftChild != null)
        {
            TreeEdge leftEdge = leftChild.getCofEdge();
            int[] leftMLambda = nextState.conformation; 
            int[] leftM = getMstateForEdgeCurState(leftMLambda, leftEdge);
            leftEdge.bTrackBestConfRemoveEarlyNew(bestPosAARot, leftM);
            nextState.updateLeftEnergy(leftEdge.peekEnergy(bestPosAARot, leftM));
            reinsert = leftEdge.moreConformations(bestPosAARot, leftM);
        }
        if(reinsert)
            outHeap.add(nextState);
        
        checkHeap(outHeap);
        if(outHeap.size() > 1)
        {
        double curNewBest = outHeap.peek().energy;
        if(curNewBest < curBest)
        {
            System.err.println("Out of order. Terminating...");
            System.exit(-1);
        }
        }
        
        if(printHeap)
        {
            String curString = RTMToString(bestPosAARot);
            PriorityQueue<Conf> copy = new PriorityQueue<Conf>();

            for(Conf c : outHeap)
            {
                copy.add(c);
            }

            while(!copy.isEmpty())
            {
                Conf c = copy.poll();
                System.out.println(c+"$"+c.energy+", left "+c.leftEnergy+", self "+c.selfEnergy+", right "+c.rightEnergy);
            }
            System.out.println("====================== end "+curString+"========================");

        }


    }

	private void getNewRightconf(RotTypeMap[] bestPosAARotOld, TreeEdge rightEdge, List<RightConf> rightConfs, int[] rightM) 
	{
		RotTypeMap[] bestPosAARotRight = Arrays.copyOf(bestPosAARotOld, bestPosAARotOld.length);
		double newRightEnergy = rightEdge.peekEnergy(bestPosAARotRight, rightM);
		rightEdge.bTrackBestConfRemoveEarlyNew(bestPosAARotRight, rightM);
		RightConf newRight = new RightConf(bestPosAARotRight, newRightEnergy);
		System.out.println("New Right Conf: "+RTMToString(bestPosAARotRight)+", energy "+newRight.energy);
		rightConfs.add(newRight);
	}


    /**
     * New NEW algorithm:
     * 
     * Remove late option:
     * 1. Given M conformation, peek lambda heap
     * 2. Use results of peek to populate output
     * 3. Recurse on left child
     * 4. If a new right conformation is needed, get one.
     * 5. If there was a right conformation,
     *    recalculate energies along this path, by traversing down to the leaf, assigning 
     *    it the new right energy, and allowing the results to propagate back up.
     * 6. Poll the current lambda heap for the root node.
     * 7. If there are remaining left or right conformations, update its energy with 
     *    the new energy from its child, and reinsert into the heap.
     * */
    public void bTrackBestConfRemoveEarly(RotTypeMap bestPosAARot[], int[] bestState, Stack<Conf> polledConfs, Stack<Boolean> reinserts)
    {
        //Peek and populate
        PriorityQueue<Conf> outHeap = getHeap(bestPosAARot, bestState);
        checkHeap(outHeap);
        String curString = RTMToString(bestPosAARot);
        /* TODO: Handle empty lambda set */
        boolean sort = false;
        if(sort)
        {
            PriorityQueue<Conf> sortCopy = new PriorityQueue<Conf>();
            while(!outHeap.isEmpty())
            {
                sortCopy.add(outHeap.poll());
            }

            while(!sortCopy.isEmpty())
            {
                Conf c = sortCopy.poll();
                outHeap.add(c);
            }
        }
        if(printHeap)
        {
            System.out.println("====================== start "+curString+"========================");
            PriorityQueue<Conf> copy = new PriorityQueue<Conf>();

            for(Conf c : outHeap)
            {
                copy.add(c);
            }

            while(!copy.isEmpty())
            {
                Conf c = copy.poll();
                System.out.println(c+"$"+c.energy+", left "+c.leftEnergy+", self "+c.selfEnergy+", right "+c.rightEnergy);
            }
        }

        Conf nextState = outHeap.poll();
        nextState.fillRotTypeMap(bestPosAARot);

        RotTypeMap[] bestPosAARot3 = Arrays.copyOf(bestPosAARot, bestPosAARot.length);

        //If leaf, return
        if(leftChild == null)
        {
            polledConfs.push(nextState);
            reinserts.push(true);
            return;
        }

        //Recurse
        TreeEdge leftEdge = leftChild.getCofEdge();
        int[] leftMLambda = nextState.conformation;
        int[] leftM = getMstateForEdgeCurState(leftMLambda, leftEdge);
        String sanityString = getLeftConfString(bestPosAARot);
        RotTypeMap last3 = bestPosAARot[3];
        leftChild.getCofEdge().bTrackBestConfRemoveEarly(bestPosAARot, leftM, polledConfs, reinserts);

        //Handle right side
        boolean reinsert = bTrackRightSideRemoveLate(bestPosAARot, leftEdge.getL(), leftMLambda);

        //Reinsertion time!
        boolean outOfConformations = !leftEdge.moreConformations(bestPosAARot3, leftM);

        if(!polledConfs.isEmpty() && reinsert)
        {
            String leftConfString = getLeftConfString(bestPosAARot);
            int index = getRightOffset(leftConfString);
            List<RightConf> rightSolutions = getRightSolutions(bestPosAARot);
            if(printHeap)
            {
                for(RightConf r : rightSolutions)
                {
                    System.out.println(r+" energy is "+r.energy);
                }
            }
            RightConf newRightConf = rightSolutions.get(index);
            polledConfs.push(nextState);
            reinserts.push(true);
            //Must remove our right side result before reinserting...
            RotTypeMap[] bestPosAARotCopy = new RotTypeMap[bestPosAARot.length];
            for(int i = 0; i < bestPosAARotCopy.length; i++)
            {
                bestPosAARotCopy[i] = bestPosAARot[i];
            }
            reinsertLeftConformation(bestPosAARotCopy, leftMLambda, polledConfs, reinserts, newRightConf.energy);
        }

        if(!reinsert && !outOfConformations)
        {
            //TODO: consolidate RTM modifications. This should go elsewhere, like into peekEnergy,
            // which should be called at the current node, not at the child.
            RotTypeMap[] bestPosAARotCopy = Arrays.copyOf(bestPosAARot, bestPosAARot.length);
            for(Integer l : L)
            {
                if(!lambda.contains(l))
                    bestPosAARotCopy[invResMap[l]] = null;
            }
            double nextLeftEnergy = leftEdge.peekEnergy(bestPosAARotCopy, leftMLambda);
            nextState.updateLeftEnergy(nextLeftEnergy);
            outHeap.add(nextState);
        }
        if(!polledConfs.isEmpty())
        {
            polledConfs.push(nextState);
            reinserts.push(!leftEdge.moreConformations(bestPosAARot3, leftM));
        }

        checkHeap(outHeap);

        if(printHeap)
        {
            PriorityQueue<Conf> copy = new PriorityQueue<Conf>();

            for(Conf c : outHeap)
            {
                copy.add(c);
            }

            while(!copy.isEmpty())
            {
                Conf c = copy.poll();
                System.out.println(c+"$"+c.energy+", left "+c.leftEnergy+", self "+c.selfEnergy+", right "+c.rightEnergy);
            }
            System.out.println("====================== end "+curString+"========================");

        }

    }

    public void bTrackBestConfRemoveLate(RotTypeMap bestPosAARot[], int[] bestState, Stack<Conf> polledConfs, Stack<Boolean> reinserts)
    {
        //Peek and populate
        PriorityQueue<Conf> outHeap = getHeap(bestPosAARot, bestState);
        checkHeap(outHeap);
        String curString = RTMToString(bestPosAARot);
        boolean sort = false;
        if(sort)
        {
            PriorityQueue<Conf> sortCopy = new PriorityQueue<Conf>();
            while(!outHeap.isEmpty())
            {
                sortCopy.add(outHeap.poll());
            }

            while(!sortCopy.isEmpty())
            {
                Conf c = sortCopy.poll();
                outHeap.add(c);
            }
        }
        if(printHeap)
        {
            System.out.println("====================== start "+curString+"========================");
            PriorityQueue<Conf> copy = new PriorityQueue<Conf>();

            for(Conf c : outHeap)
            {
                copy.add(c);
            }

            while(!copy.isEmpty())
            {
                Conf c = copy.poll();
                System.out.println(c+"$"+c.energy+", left "+c.leftEnergy+", self "+c.selfEnergy+", right "+c.rightEnergy);
            }
        }

        Conf nextState = outHeap.peek();
        nextState.fillRotTypeMap(bestPosAARot);

        RotTypeMap[] bestPosAARot3 = Arrays.copyOf(bestPosAARot, bestPosAARot.length);

        //If leaf, return
        if(leftChild == null)
        {
            polledConfs.push(nextState);
            reinserts.push(true);
            return;
        }

        //Recurse
        TreeEdge leftEdge = leftChild.getCofEdge();
        int[] leftMLambda = nextState.conformation;
        int[] leftM = getMstateForEdgeCurState(leftMLambda, leftEdge);
        String sanityString = getLeftConfString(bestPosAARot);
        RotTypeMap last3 = bestPosAARot[3];
        leftChild.getCofEdge().bTrackBestConfRemoveLate(bestPosAARot, leftM, polledConfs, reinserts);

        //Handle right side
        boolean reinsert = bTrackRightSideRemoveLate(bestPosAARot, leftEdge.getL(), leftMLambda);

        //Reinsertion time!
        boolean outOfConformations = !leftEdge.moreConformations(bestPosAARot3, leftM);

        if(!polledConfs.isEmpty() && reinsert)
        {
            String leftConfString = getLeftConfString(bestPosAARot);
            int index = getRightOffset(leftConfString);
            List<RightConf> rightSolutions = getRightSolutions(bestPosAARot);
            if(printHeap)
            {
                for(RightConf r : rightSolutions)
                {
                    System.out.println(r+" energy is "+r.energy);
                }
            }
            RightConf newRightConf = rightSolutions.get(index);
            polledConfs.push(nextState);
            reinserts.push(true);
            //Must remove our right side result before reinserting...
            RotTypeMap[] bestPosAARotCopy = new RotTypeMap[bestPosAARot.length];
            for(int i = 0; i < bestPosAARotCopy.length; i++)
            {
                bestPosAARotCopy[i] = bestPosAARot[i];
            }
            reinsertLeftConformation(bestPosAARotCopy, leftMLambda, polledConfs, reinserts, newRightConf.energy);
        }

        if(!reinsert && !outOfConformations)
        {
            //TODO: consolidate RTM modifications. This should go elsewhere, like into peekEnergy,
            // which should be called at the current node, not at the child.
            RotTypeMap[] bestPosAARotCopy = Arrays.copyOf(bestPosAARot, bestPosAARot.length);
            for(Integer l : L)
            {
                if(!lambda.contains(l))
                    bestPosAARotCopy[invResMap[l]] = null;
            }
            double nextLeftEnergy = leftEdge.peekEnergy(bestPosAARotCopy, leftMLambda);
            nextState.updateLeftEnergy(nextLeftEnergy);
            outHeap.add(nextState);
        }
        if(!polledConfs.isEmpty())
        {
            polledConfs.push(nextState);
            reinserts.push(!leftEdge.moreConformations(bestPosAARot3, leftM));
        }

        checkHeap(outHeap);

        if(printHeap)
        {
            PriorityQueue<Conf> copy = new PriorityQueue<Conf>();

            for(Conf c : outHeap)
            {
                copy.add(c);
            }

            while(!copy.isEmpty())
            {
                Conf c = copy.poll();
                System.out.println(c+"$"+c.energy+", left "+c.leftEnergy+", self "+c.selfEnergy+", right "+c.rightEnergy);
            }
            System.out.println("====================== end "+curString+"========================");

        }

    }

    private void checkHeap(PriorityQueue<Conf> heap)
    {
        if(heap.isEmpty()) return;
        PriorityQueue<Conf> copy = new PriorityQueue<Conf>();
        double lastEnergy = heap.peek().energy;
        copy.add(heap.poll());
        while(!heap.isEmpty())
        {
            Conf c = heap.peek();
            if(c.energy < lastEnergy)
            {
                System.out.println("AHHHHH");
            }
            copy.add(c);
            lastEnergy = c.energy;
            heap.poll();
        }

        while(!copy.isEmpty())
        {
            heap.add(copy.poll());
        }
    }

    /** 
     * Algorithm to prepopulate right energies:
     * 1. If a right child exists, recurse right and get its best energy.
     * 2. Recurse left.
     * 3. If at a left leaf, assign the leaf that right energy.
     * **/

    public void populateLeftHeaps(RotTypeMap[] bestPosAARot, int[] MLambda, double rightEnergy)
    {
        //System.out.println(RTMToString(bestPosAARot));
        double energy;
        if(leftChild == null)
        {
            PriorityQueue<Conf> heap = A2.get(computeIndexInA(MLambda));
            Stack<Conf> readd = new Stack<Conf>();
            while(!heap.isEmpty())
            {
                Conf c = heap.poll();
                readd.push(c);
                c.updateRightEnergy(rightEnergy);
            }
            heap.addAll(readd);
            return;
        }
        if(rightChild != null)
        {
            TreeEdge rightEdge = rightChild.getCofEdge();
            PriorityQueue<Conf> heap = A2.get(computeIndexInA(MLambda));
            Stack<Conf> readd = new Stack<Conf>();
            while(!heap.isEmpty())
            {
                Conf c = heap.poll();
                readd.push(c);
                RotTypeMap[] outCopy = Arrays.copyOf(bestPosAARot, bestPosAARot.length);
                c.fillRotTypeMap(outCopy);
                int[] rightMLambda = getMstateForEdgeCurState(c.conformation, rightEdge);
                rightEdge.populateLeftHeaps(outCopy, rightMLambda, 0);
                energy = rightEdge.peekEnergy(outCopy, rightMLambda);
                TreeEdge leftEdge = leftChild.getCofEdge();
                // TODO: This has to split for right and left, and then update left as well.
                // TODO: investigate correctness.
                int[] leftMLambda = getMstateForEdgeCurState(c.conformation, leftEdge);
                leftEdge.populateLeftHeaps(outCopy, leftMLambda, rightEnergy + energy);
                double leftEnergy = leftEdge.peekEnergy(outCopy, leftMLambda);
                c.updateLeftEnergy(leftEnergy);

            }
            heap.addAll(readd);

            checkHeap(heap);
        }
        else
        {
            PriorityQueue<Conf> heap = A2.get(computeIndexInA(MLambda));
            Stack<Conf> readd = new Stack<Conf>();
            while(!heap.isEmpty())
            {
                Conf c = heap.poll();
                readd.push(c);
                RotTypeMap[] outCopy = Arrays.copyOf(bestPosAARot, bestPosAARot.length);
                c.fillRotTypeMap(outCopy);
                TreeEdge leftEdge = leftChild.getCofEdge();
                int[] leftMLambda = getMstateForEdgeCurState(c.conformation, leftEdge);
                leftEdge.populateLeftHeaps(outCopy, leftMLambda, rightEnergy);
                double leftEnergy = leftEdge.peekEnergy(outCopy, leftMLambda);
                c.updateLeftEnergy(leftEnergy);

            }
            heap.addAll(readd);

            checkHeap(heap);
        }
    }

    public void populateLeftHeaps()
    {
        populateLeftHeaps(new RotTypeMap[molResMap.length], new int[]{}, 0);
    }


    private int getRightOffset(String leftConfString) {
        if(!rightSolutionOffset.containsKey(leftConfString))
            rightSolutionOffset.put(leftConfString, 0);
        int index = rightSolutionOffset.get(leftConfString);
        return index;
    }

    private String getLeftConfString(RotTypeMap[] bestPosAARot) {
        RotTypeMap[] bestPosAARot2 = Arrays.copyOf(bestPosAARot, bestPosAARot.length);
        for(int l : rightL)
        {
            bestPosAARot2[invResMap[l]] = null;
        }
        String leftConfString = RTMToString(bestPosAARot2);
        return leftConfString;
    }

    private PriorityQueue<Conf> getHeap(RotTypeMap[] bestPosAARot, int[] bestState)
    {
        String curString = RTMToPrefix(bestPosAARot);
        if(!leftHeapMap.containsKey(curString))
        {
            PriorityQueue<Conf> newHeap = new PriorityQueue<Conf>();
            if(A2.size() < 1)
            {
                if(lambda.size() < 1)
                {
                    newHeap.add(new Conf(new int[]{}, 0, rtm));
                }
                else
                {
                    System.err.println("Unrecognized sequence with existing lambda set...");
                }
            }
            else
            {
                PriorityQueue<Conf> outHeap = A2.get(computeIndexInA(bestState));
                checkHeap(outHeap);

                for(Conf c : outHeap)
                {
                    Conf newConf = c.copy();
                    newHeap.add(newConf);
                }
                checkHeap(newHeap);
            }
            leftHeapMap.put(curString, newHeap);
        }
        PriorityQueue<Conf> out = leftHeapMap.get(curString);
        checkHeap(out);
        return out;
    }

    private LazyHeap<Conf> getSecondaryHeap(RotTypeMap[] bestPosAARot, int[] bestState)
    {
        //System.err.println("UNFINISHED CODE. TERMINATE.");
        //System.exit(-1);
        String curString = RTMToPrefixLambda(bestPosAARot);
        if(!secondaryHeapMap.containsKey(curString))
        {
            PriorityQueue<Conf> template = leftChild.getCofEdge().getHeap(bestPosAARot, bestState);
            LazyHeap<Conf> newHeap = new LazyHeap<Conf>(template);
            secondaryHeapMap.put(curString, newHeap);
        }
        LazyHeap<Conf> out = secondaryHeapMap.get(curString);
        checkHeap(out);
        return out;
    }


    private String RTMToString(RotTypeMap[] rtm)
    {
        String output = "";
        for(int i = 0; i < rtm.length; i++)
        {
            RotTypeMap current  = rtm[i];

            if(current!=null)
                output+= current.pos+":"+current.aa+"-"+current.rot+" ";
        }        
        return output;
    } 

    private boolean bTrackRightSideRemoveLate (RotTypeMap[] bestPosAARot, LinkedHashSet<Integer> leftL,
            int[] leftMLambda) {
        /* TODO: Handle empty lambda set */
        if(rightChild != null)
        {
            String leftConfString = getLeftConfString(bestPosAARot);
            int index = getRightOffset(leftConfString);
            List<RightConf> rightSolutions = getRightSolutions(bestPosAARot);
            while(index + 1 >= rightSolutions.size())
            {
                if(printHeap)
                    System.out.println("=================================CREATING NEW RIGHT CONFORMATION=============================================");
                //TODO: Create new maps on the fly, MAKE SURE ITS THE M SET!!
                RotTypeMap[] bestPosAARot3 = Arrays.copyOf(bestPosAARot, bestPosAARot.length);
                for(int l : leftL)
                {
                    bestPosAARot3[invResMap[l]] = null;
                }
                TreeEdge rightEdge = rightChild.getCofEdge();
                int[] rightM = getMstateForEdgeCurState(leftMLambda, rightEdge);

                if(!rightEdge.moreConformations(bestPosAARot3, rightM))
                {
                    System.out.println(leftConfString + "is depleted.");
                    break;
                }
                double rightEnergy = rightEdge.peekEnergy(bestPosAARot3, rightM);
                rightEdge.bTrackBestConfRemoveEarly(bestPosAARot3, rightM, new Stack<Conf>(), new Stack<Boolean>());
                RightConf conf = new RightConf(bestPosAARot3, rightEnergy);
                if(printHeap)
                    System.out.println("New energy is "+rightEnergy+", conf "+conf);
                rightSolutions.add(conf);
            }
            if(printHeap)
            {
                for(RightConf c : rightSolutions)
                {
                    System.out.println(c+", energy "+c.energy);
                }
            }
            if(index < rightSolutions.size())
            {
                RightConf rightConf = rightSolutions.get(index);
                rightConf.fillRotTypeMap(bestPosAARot);
                rightSolutionOffset.put(leftConfString, index + 1);
            }

            return index + 1 < rightSolutions.size();
        }
        return false;
    }

    private List<RightConf> getRightSolutions (RotTypeMap[] bestPosAARot) {
        String prefix = RTMToPrefixLambda(bestPosAARot);
        if(!rightSolutionMap.containsKey(prefix))
            rightSolutionMap.put(prefix, new LinkedList<RightConf>());
        List<RightConf> rightSolutions = rightSolutionMap.get(prefix);
        return rightSolutions;
    }

    private void reinsertLeftConformation(RotTypeMap bestPosAARot[], int[] bestState, Stack<Conf> removedConfs, Stack<Boolean> reinserts, double newRightEnergy)
    {
        /* TODO: Handle empty lambda set */
        if(rightChild != null)
        {
            TreeEdge rightEdge = rightChild.getCofEdge();
            for(int l: rightEdge.getL())
            {
                bestPosAARot[invResMap[l]] = null;
            }
        }
        //Peek
        PriorityQueue<Conf> outHeap = getHeap(bestPosAARot, bestState);

        Conf toAdd = removedConfs.pop();
        boolean reinsert = reinserts.pop();

        //If leaf, return
        if(leftChild == null)
        {
            toAdd.updateRightEnergy(newRightEnergy);
            if(printHeap)
                System.out.println("Updated leaf with "+newRightEnergy+" new energy is "+toAdd.energy);
            outHeap.add(toAdd);
            checkHeap(outHeap);

            if(printHeap)
            {
                System.out.println("====================== new Heap: ========================");
                PriorityQueue<Conf> copy = new PriorityQueue<Conf>();

                for(Conf c : outHeap)
                {
                    copy.add(c);
                }

                while(!copy.isEmpty())
                {
                    Conf c = copy.poll();
                    System.out.println(c+"$"+c.energy+", left "+c.leftEnergy+", self "+c.selfEnergy+", right "+c.rightEnergy);
                }
            }
            return;
        }

        //Recurse
        TreeEdge leftEdge = leftChild.getCofEdge();
        int[] leftMLambda = toAdd.conformation;
        int[] leftM = getMstateForEdgeCurState(leftMLambda, leftEdge);

        if(rightChild != null)
        {
            String leftConfString = getLeftConfString(bestPosAARot);
            int index = getRightOffset(leftConfString);

            List<RightConf> rightSolutions = getRightSolutions(bestPosAARot);
            if(index >= rightSolutions.size())
            {
                double resetRightEnergy = rightSolutions.get(0).energy;
                newRightEnergy += resetRightEnergy;
                rightSolutionOffset.put(bestPosAARot.toString(), 0);
            }
        }

        leftEdge.reinsertLeftConformation(bestPosAARot, leftM, removedConfs, reinserts, newRightEnergy);
        double newLeftEnergy = leftEdge.peekEnergy(bestPosAARot, leftM);
        if(toAdd.leftEnergy > newLeftEnergy)
        {
            outHeap.remove(toAdd);
            reinsert = true;
        }
        toAdd.updateLeftEnergy(newLeftEnergy);
        if(reinsert){
            if(!outHeap.isEmpty() && outHeap.peek().equals(toAdd))
                System.out.println("=========AHHH ABORT================\n"
                        + "do  not reinsert "+toAdd+" to "+outHeap
                        + "\n======================"+
                        "===========\n=================");
            outHeap.add(toAdd);
            checkHeap(outHeap);
        }

        if(printHeap)
        {
            System.out.println("====================== new Heap: ========================");
            PriorityQueue<Conf> copy = new PriorityQueue<Conf>();

            for(Conf c : outHeap)
            {
                copy.add(c);
            }

            while(!copy.isEmpty())
            {
                Conf c = copy.poll();
                System.out.println(c+"$"+c.energy+", left "+c.leftEnergy+", self "+c.selfEnergy+", right "+c.rightEnergy);
            }
        }
    }



    public boolean moreConformations(RotTypeMap[] bestPosAARot, int[] state)
    {
        PriorityQueue<Conf> outHeap = getHeap(bestPosAARot, state);
        return !outHeap.isEmpty();
    }


    public double peekEnergy(RotTypeMap[] bestPosAARot, int[] MLambda)
    {
        PriorityQueue<Conf> heap = getHeap(bestPosAARot, MLambda); 
        if(heap.size() < 1) return Double.MAX_VALUE;
        return heap.peek().energy;
    }

    private String RTMToLString (RotTypeMap[] bestPosAARot) {
        String output = "[";
        //TODO: STORE THESE!!!
        Integer[] MArray = M.toArray(new Integer[]{});
        Integer[] LArray = leftL.toArray(new Integer[]{});
        Integer[] outArray = new Integer[MArray.length + LArray.length];
        for(int i = 0; i < MArray.length; i++)
        {
            outArray[i] = MArray[i]; 
        }
        for(int i = 0; i < LArray.length; i++)
        {
            outArray[i+MArray.length] = LArray[i];
        }
        Arrays.sort(outArray);
        for(int i : outArray)
        {
            RotTypeMap current  = bestPosAARot[invResMap[i]];

            if(current!=null)
                output+= current.pos+":"+current.aa+"-"+current.rot+" ";
        }
        //System.out.println("Constructed "+output);
        return output +"]";
    }

    private String RTMToPrefix (RotTypeMap[] bestPosAARot) {
        String output = "[";
        String[] prefixElements = new String[bestPosAARot.length];

        for(RotTypeMap current: bestPosAARot)
        {
            if(current!=null)
                prefixElements[current.pos] = current.pos+":"+current.aa+"-"+current.rot+" ";
        }

        for(int l : L)
        {
            prefixElements[invResMap[l]] = "";
        }

        for(int i = 0; i < prefixElements.length; i++)
            if(prefixElements[i]!=null)
                output+=prefixElements[i];
        return output +"]";
    }

    private String RTMToPrefixLambda (RotTypeMap[] bestPosAARot) {
        String output = "[";
        String[] prefixElements = new String[bestPosAARot.length];

        for(RotTypeMap current: bestPosAARot)
        {
            if(current!=null)
                prefixElements[current.pos] = current.pos+":"+current.aa+"-"+current.rot+" ";
        }

        for(int l : L)
        {
            prefixElements[invResMap[l]] = "";
        }

        for(int l : lambda)
        {
            RotTypeMap current = bestPosAARot[invResMap[l]];
            prefixElements[invResMap[l]] = current.pos+":"+current.aa+"-"+current.rot+" ";
        }

        for(int i = 0; i < prefixElements.length; i++)
            if(prefixElements[i]!=null)
                output+=prefixElements[i];
        return output +"]";
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

    public boolean equals (Conf other)
    {
        return toString().equals(other.toString());
    }

    private class RightConf
    {
        RotTypeMap[] fullConformation;
        double energy;

        public RightConf(RotTypeMap[] fullRightConf, double e)
        {
            fullConformation = fullRightConf;
            energy = e;
        }

        public void fillRotTypeMap(RotTypeMap[] bestPosAARot)
        {
            for(int i=0; i<fullConformation.length;i++)
            {
                if(fullConformation[i] != null)
                {
                    RotTypeMap result = fullConformation[i];
                    RotTypeMap previous = bestPosAARot[i];
                    if(previous != null)
                    {
                        if (result.aa != previous.aa)
                            System.out.println("AAH!!! OVERWRITE "+previous.aa+" with "+result.aa);
                        if(result.rot != previous.rot)
                            System.out.println("AAH!!! OVERWRITE "+previous.rot+" with "+result.rot);
                    }
                    bestPosAARot[i] = fullConformation[i];

                }
            }
        }

        public String toString()
        {
            String out = "[";
            RotTypeMap[] conformation = fullConformation;
            for(int i=0; i<conformation.length;i++)
            {
                RotTypeMap current = conformation[i];
                if(current != null)
                    out+= current.pos+":"+current.aa+"-"+current.rot+" ";
            }
            out+="]";
            return out;
        }

    }

    private class Conf implements Comparable<Conf>
    {
        int[] conformation;
        double energy;
        double selfEnergy;
        double leftEnergy;
        double rightEnergy;
        RotTypeMap[] conf;
        int indexInA;
        RotTypeMap[][] rtm;
        public Conf(int[] c, double e, RotTypeMap[][] matrix)
        {
            conformation = c;
            selfEnergy = e;
            energy = e;
            rtm = matrix;
        }
        
        public void fillConf(RotTypeMap[] bestPosAARot)
        {
            for(int i=0; i<conf.length;i++)
            {
                if(conf[i] != null)
                {
                    RotTypeMap result = conf[i];
                    RotTypeMap previous = bestPosAARot[i];
                    if(previous != null)
                    {
                        if (result.aa != previous.aa)
                            System.out.println("AAH!!! OVERWRITE "+previous.aa+" with "+result.aa);
                        if(result.rot != previous.rot)
                            System.out.println("AAH!!! OVERWRITE "+previous.rot+" with "+result.rot);
                    }
                    bestPosAARot[i] = conf[i];

                }
            }
        }

        public Conf(RotTypeMap[] conformation, double energy)
        {
            conf = conformation;
            selfEnergy = energy;
        }

        public Conf copy()
        {
            Conf out = new Conf(conformation, selfEnergy, rtm);
            out.rightEnergy = rightEnergy;
            out.leftEnergy = leftEnergy;
            out.energy = energy;
            return out;
        }

        public void fillRotTypeMap(RotTypeMap[] bestPosAARot)
        {
            for(int i=0; i<conformation.length;i++)
            {   
                int position=rtm[i][conformation[i]].pos;
                RotTypeMap result =  new RotTypeMap(rtm[i][conformation[i]].pos,rtm[i][conformation[i]].aa,rtm[i][conformation[i]].rot);
                RotTypeMap previous = bestPosAARot[position];
                if(previous != null)
                {
                    if (result.aa != previous.aa)
                        System.out.println("AAH!!! OVERWRITE "+previous.aa+" with "+result.aa);
                    if(result.rot != previous.rot)
                        System.out.println("AAH!!! OVERWRITE "+previous.rot+" with "+result.rot);
                }
                bestPosAARot[position] = result;
            }
        }

        public void updateLeftEnergy(double newLeftEnergy)
        {
            leftEnergy = newLeftEnergy;
            energy = leftEnergy + selfEnergy + rightEnergy;
        }

        public void updateRightEnergy(double newRightEnergy)
        {
            rightEnergy = newRightEnergy;
            energy = leftEnergy + selfEnergy + rightEnergy;
        }

        public String toString()
        {
            if(conformation == null || conformation.length < 1)
            {
                return toConfString();
            }
            String out = "[";
            for(int i=0; i<conformation.length;i++)
            {
                RotTypeMap current = rtm[i][conformation[i]];
                out+= current.pos+":"+current.aa+"-"+current.rot+" ";
            }
            out+="]";
            return out;
        }
        public String toConfString()
        {
            String out = "[";
            RotTypeMap[] conformation = conf;
            if(conformation != null)
            for(int i=0; i<conformation.length;i++)
            {
                RotTypeMap current = conformation[i];
                if(current != null)
                    out+= current.pos+":"+current.aa+"-"+current.rot+" ";
            }
            out+="]";
            return out;
        }

        public int hashCode()
        {
            return indexInA;
        }

        @Override
        public int compareTo(Conf c) {
            if(energy - c.energy < 0) 
                return -1;
            if(energy - c.energy > 0)
                return 1;
            return 0;
        }
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
    
    private class LazyHeap<T> extends PriorityQueue<T>
    {
        private boolean dirty = false;
        private T cleanNode;
        private PriorityQueue<T> templateHeap;
        
        public LazyHeap (PriorityQueue<T> template)
        {
            templateHeap = template;
        }
        
        public T lazyPoll()
        {
            if(size() < 1 || dirty)
            {
                cleanNode = templateHeap.poll();
                add(cleanNode);
            }
            
            T next = super.poll();
            if(next == cleanNode)
            {
                dirty = true;
            }
            
            
            return next;
        }
    }

    public void setM(LinkedHashSet<Integer> lambda2) 
    {
        M = lambda2;
    }

}
