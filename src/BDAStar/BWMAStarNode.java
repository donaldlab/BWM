package BDAStar;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import kstar.TreeNode;

/*
 * TODO:
1. BranchingNode: Stores left subtree and right subtree, and creates a linked list of 
values for one subtree to point down.
2. Recalculating node: Do we want to count the top K all at once? that's kq at every
level but it never exceeds 3kq, I think
3. A* tree node: expands everything, just does it smarter because
everything is precomputed.

GetNextConformation function is different.

AbstractBDAStarNode
ListingBDAStarNode
TopKBDAStarNode
ExpandingBDAStarNode
 */
public class BWMAStarNode implements Comparable<BWMAStarNode> {

    private PriorityQueue<BWMAStarNode> children;
    private BWMAStarNode leftSubtree;
    private BWMAStarNode rightSubtree;
    private LinkedList<Conformation> solutionList;
    private Map<String, Integer> solutions;
    private Conformation partialConformation;
    boolean isLeaf = false;
    boolean branching = false;
    
    public BWMAStarNode(Conformation conf) 
    {
        children = new PriorityQueue<BWMAStarNode>();
        solutionList = new LinkedList<Conformation>();
        solutions = new HashMap<String, Integer>();
        partialConformation = conf;
    }
    
    
    public static BWMAStarNode CreateTree(TreeNode root, Conformation previous, SolutionSpace s)
    {
        List<Position> lambda = root.getCofEdge().getPositionList();
        BWMAStarNode AStarRoot = new BWMAStarNode(previous);
        AStarRoot.populateHeap(AStarRoot.children, lambda, 0, previous, s);
        if(!root.getIsLeaf()) //branch 
        {
        	for(BWMAStarNode child : AStarRoot.children)
        		child.addBranches(CreateTree(root.getlc(), AStarRoot.partialConformation, s), CreateTree(root.getrc(), AStarRoot.partialConformation, s));
        }
        else AStarRoot.isLeaf = true;
        return AStarRoot;
    }
    
    public void populateHeap(PriorityQueue<BWMAStarNode> heap, List<Position> positions, int index, Conformation currentConf, SolutionSpace s)
    {
		for(Choice c : s.getChoices(positions.get(index)))
		{
			Conformation nextConf = currentConf.copy();
			nextConf.append(positions.get(index), c);
			System.out.println("Creating conformation " + nextConf);
			if(index == positions.size() - 1)
				heap.add(new BWMAStarNode(nextConf));
			else 
				populateHeap(heap, positions, index+1, nextConf, s);
		}
    }
    
    public void resort()
    {
    	for(BWMAStarNode b : children)
    	{
    		b.resort();
    	}
    	if(branching)
    	{
    		leftSubtree.resort();
    		rightSubtree.resort();
    	}
    	PriorityQueue<BWMAStarNode> newHeap = new PriorityQueue<BWMAStarNode>();
    	newHeap.addAll(children);
    	children = newHeap;
    }
    
    private Conformation getConformation () {
        return partialConformation;
    }
    
    private void addBranches(BWMAStarNode leftChild, BWMAStarNode rightChild)
    {
        leftSubtree = leftChild;
        rightSubtree = rightChild;
        branching = true;
    }
    
    public double nextBestScore()
    {
        /* Redundant code paths 
        if(isLeaf) //This case could only occur at the leaf of a single tree.
            return partialConformation.score();
        if(branching)
            return leftSubtree.peekNextConformation().join(rightSubtree.peekNextConformation()).score();
            */
        return peekNextConformation().score();
    }

    private Conformation peekNextConformation () {
    	System.out.println(this+", conformation "+partialConformation+" called for heuristic...");
        if(branching){
        	Conformation peeked = leftSubtree.peekNextConformation().join(partialConformation);
        	if(peeked == null){
        		System.out.println("No left sub tree");
        		return partialConformation;
        	}
        	if(!solutions.containsKey(peeked.toString())) {
        		System.out.println("Returning default optima"+ peeked.join(rightSubtree.peekNextConformation()));
        		return peeked.join(rightSubtree.peekNextConformation()).join(partialConformation);
        	}
        	int offset = solutions.get(peeked.toString());
        	Conformation rightSide = null;
        	if(offset >= solutionList.size()){
        		rightSide = rightSubtree.peekNextConformation();
        		System.out.println("Getting next best optima:" + partialConformation.join(peeked.join(rightSide)));
        	}
        	else rightSide = solutionList.get(offset);
        	if(rightSide == null && rightSubtree.moreConformations()) 
        		rightSide = rightSubtree.peekNextConformation();
        	else {
        		System.out.println("Eliminate node from left subtree");
        		leftSubtree.getNextConformation();
        		return peekNextConformation();
        	}
        	if(rightSide == null){
        		System.out.println("No right subtree");
        		return partialConformation;
        	}
            return partialConformation.join(peeked.join(rightSide));
        }
        if(children.size() < 1)
        {
        	return partialConformation;
        }
        if(isLeaf) 
        {
        	return children.peek().partialConformation.join(partialConformation);
        }

        BWMAStarNode peek = children.peek();
        //System.out.println("Return child score: "+children.peek().peekNextConformation().join(partialConformation));
        return children.peek().peekNextConformation().join(partialConformation);
    }
    
    public boolean moreConformations()
    {
    	if(isLeaf) return true;
    	if(branching) return rightSubtree.moreConformations() || 
    			(solutionList.size() > solutions.get(leftSubtree.peekNextConformation().toString()));
    	boolean more = false;
    	for(BWMAStarNode b : children)
    	{
    		if(b.moreConformations())
    			more = true;
    	}
    	return more;
    }
    
    public boolean branchHasNext()
    {
        return leftSubtree.children.size() > 1;
    }
     
    public Conformation getNextConformation()
    {
        return getNextConformation(false);
    }
    
    public Conformation getNextConformation (boolean reinsert) 
    {
        System.out.println("PROCESS NODE: "+partialConformation+", score: "+partialConformation.score()+", heuristic: "+peekNextConformation()+" = "+nextBestScore());
        //System.out.println("Children size: "+children.size());
        

    	if(branching){
            //System.out.println("Process branching.");
            /*
             * We need to keep track of children and not remove them until they are all had....
             */
            BWMAStarNode leftChild = leftSubtree;
            if(leftChild == null)
            {
                //System.out.println("...wha?");
            }
            Conformation leftConformation = leftChild.getNextConformation(true);
            if(solutions.get(leftConformation.toString())==null){
                //System.out.println("Initialize "+leftConformation);
                solutions.put(leftConformation.toString(), 0);
            }
            //System.out.println("Starting "+leftConformation+" at "+solutions.get(leftConformation.toString()));
            Iterator<Conformation> pointer = solutionList.listIterator(solutions.get(leftConformation.toString()));
            //System.out.println ("Branch!");
            
            
            
            Conformation rightConformation = null;
            if(pointer.hasNext()){
                rightConformation = pointer.next();
            }
            else
            {
                //System.out.println("polling for new right conformation...");
                rightConformation = rightSubtree.getNextConformation(false);
                if(rightConformation != null){ 
                    //put child back in, it has more to go!
                    //System.out.println("Reinsert "+leftChild.partialConformation);
                    solutionList.add(rightConformation);
                }
                else // no more conformations
                {
                    //System.out.println("Right subtree is depleted, removing "+leftConformation);
                    leftSubtree.getNextConformation();
                    solutions.put(leftConformation.toString(), (solutions.get(leftConformation.toString())+1));
                    
                    return getNextConformation(reinsert);
                }
            }
            solutions.put(leftConformation.toString(), (solutions.get(leftConformation.toString())+1));
            
            return partialConformation.join(leftConformation.join(rightConformation));
        }
        if(isLeaf){
            System.out.println("Process leaf!.");
            if(children.size()<1) {
                //System.out.println("ITS THE END OF THE WORLD!!!");
                return null;
            }
        	BWMAStarNode out = children.poll();
        	if(reinsert) 
        		children.add(out);
            return out.partialConformation.join(partialConformation);
        }
        if(children.size() < 1)
        {
        	return null;
        }
        //System.out.println("Process intermediate node...");
        BWMAStarNode next = children.poll();
        if(next.children.size() <= 1 && !next.branching && !next.isLeaf)
        {
            //System.out.println("HOLY WHAT BALLS??!?!?!===");
        }
        //System.out.println("Polled: "+next.partialConformation);
        ////System.out.println("Next conformation is: "+next.partialConformation+", score "+next.partialConformation.score());
        //update next conformation

        Conformation nextConf = next.getNextConformation(reinsert).join(partialConformation);
        if(reinsert || next.children.size() > 0 || (next.branching && next.branchHasNext()))
        {
            //System.out.println("Reinsert"+next.partialConformation+": "+next.children.size()+" children, branching is "+next.branching);
            if(next.children.size()<=1 && !next.isLeaf && !next.branching)
            {
                //System.out.println("EMPTY CHILD INSERTION");
            }
            children.add(next);
        }
        //else System.out.println("Not reinserting "+next.partialConformation+", "+children.size()+" children remaining.");
        //System.out.println("intermediate node complete");
        
        return nextConf;
    }
    
    public void printTree(String prefix, Conformation c)
    {
        Conformation joined = partialConformation.join(c);
        System.out.println(prefix+joined);
        for(BWMAStarNode child: children)
        {
                child.printTree(prefix+"+--",joined);
        }
        if(branching)
        {
                leftSubtree.printTree(prefix+"L--",joined);
                rightSubtree.printTree(prefix+"R--",joined);
        }
    }
    
    public int remainingConformations()
    {
        if(isLeaf)
            return children.size();
        if(branching)
            return leftSubtree.children.size()*rightSubtree.children.size();
        int sum = 0;
        for(BWMAStarNode b: children)
        {
            sum += b.remainingConformations();
        }
        return sum;
    }
    
    public void printTree(String prefix)
    {
        printTree(prefix, partialConformation);
    }
    
    public int compareTo(BWMAStarNode node)
    {
        if(nextBestScore() - node.nextBestScore() < 0) return -1;
        if(nextBestScore() - node.nextBestScore() > 0) return 1;
        return 0;
    }
   
}
