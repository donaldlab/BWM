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
    
    boolean dirty = true;
    private Conformation cachedConf;
    
    private Conformation rightSideConformation;
    
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
        	for(BWMAStarNode child : AStarRoot.children){
        		child.addBranches(CreateTree(root.getlc(), child.partialConformation, s), CreateTree(root.getrc(), child.partialConformation, s));
        		child.solutionList.add(child.rightSubtree.getNextConformation());
        		for(BWMAStarNode leftChild : child.leftSubtree.children)
        		{
        		    leftChild.rightSideConformation = child.solutionList.getFirst();
        		}
        	}
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
    	cachedConf = peekNextConformation();
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
    
    private int getRightConformation(Conformation leftConformation)
    {
        if(!solutions.containsKey(leftConformation.toString()))
                solutions.put(leftConformation.toString(), 0);
        return solutions.get(leftConformation.toString());
    }

    private Conformation peekNextConformation () {
    	//System.out.println(this+", conformation "+partialConformation+" called for heuristic...");

        if(branching){
        	Conformation peeked = leftSubtree.peekNextConformation();

                if(!solutions.containsKey(peeked.toString()))
                {
                    solutions.put(peeked.toString(), 0);
                }
                if(solutionList.size() <= solutions.get(peeked.toString()) && rightSubtree.children.size() > 0) {
        		System.out.println("Returning default optima"+ peeked.join(rightSubtree.peekNextConformation()));
        		cachedConf = peeked.join(rightSubtree.peekNextConformation());
        		dirty = false;
        		return peeked.join(rightSubtree.peekNextConformation());
        	}
        	int offset = solutions.get(peeked.toString());
        	Conformation rightSide = null;
        	if(offset >= solutionList.size()){
        		rightSide = rightSubtree.peekNextConformation();
        		System.out.println("Getting next best optima:" + peeked.join(rightSide));
        	}
        	else rightSide = solutionList.get(offset);
        	//System.out.println("Calculating heuristic for "+partialConformation+", left side is "+peeked+", offset is "+offset+", solution is "+peeked+", score "+peeked.score());
                
        	if(rightSide == null){
        		System.out.println("No right subtree");
        	            Conformation out = partialConformation;
        	            if(rightSideConformation != null)
        	            {
        	                out = out.join(rightSideConformation);
        	            }
        	            cachedConf = out;
        	            dirty = false;
        	            return out;
        	}
            return peeked;
        }
        if(children.size() < 1)
        {
            Conformation out = partialConformation;
            if(rightSideConformation != null)
            {
                out = out.join(rightSideConformation);
            }

            return out;
        }
        if(isLeaf) 
        {
            Conformation out = children.peek().peekNextConformation();
            if(rightSideConformation != null)
            {
                out = out.join(rightSideConformation);
            }            
            return out;
        }

        BWMAStarNode peek = children.peek();
        //System.out.println("Return child score: "+children.peek().peekNextConformation().join(partialConformation));
        Conformation out = children.peek().peekNextConformation().join(partialConformation);
        if(rightSideConformation != null)
        {
            out = out.join(rightSideConformation);
        }

        return out;
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
    
    public Conformation peekPartial()
    {
        if(branching){
            Conformation peeked = leftSubtree.peekNextConformation();

            if(!solutions.containsKey(peeked.toString()))
            {
                solutions.put(peeked.toString(), 0);
            }
            if(solutionList.size() <= solutions.get(peeked.toString()) && rightSubtree.children.size() > 0) {
                    System.out.println("Returning default optima"+ peeked);
                    return peeked;
            }
            int offset = solutions.get(peeked.toString());
            Conformation rightSide = null;
            if(offset >= solutionList.size()){
                    rightSide = rightSubtree.peekNextConformation();
                    System.out.println("Getting next best optima:" + peeked);
            }
            else rightSide = solutionList.get(offset);
            //System.out.println("Calculating heuristic for "+partialConformation+", left side is "+peeked+", offset is "+offset+", solution is "+peeked+", score "+peeked.score());
            
            if(rightSide == null){
                    System.out.println("No right subtree");
                        Conformation out = partialConformation;

                        return out;
            }
        return peeked;
    }
    if(children.size() < 1)
    {
        Conformation out = partialConformation;

        return out;
    }
    if(isLeaf) 
    {
        Conformation out = children.peek().peekPartial();

        return out;
    }

    BWMAStarNode peek = children.peek();
    //System.out.println("Return child score: "+children.peek().peekNextConformation().join(partialConformation));
    Conformation out = children.peek().peekNextConformation();
    return out;
    }
    
    public Conformation getNextConformation (boolean reinsert) 
    {
        //System.out.println("PROCESS NODE: "+partialConformation+", score: "+partialConformation.score()+", heuristic: "+peekNextConformation()+" = "+nextBestScore());
        //System.out.println("Children size: "+children.size());

    	if(branching){
            //System.out.println("Process branching.");
            /*
             * We need to keep track of children and not remove them until they are all had....
             */
            Conformation leftConformation = leftSubtree.peekPartial();
            if(!solutions.containsKey(leftConformation.toString()))
            {
                //System.out.println("Initialize "+leftConformation);
                solutions.put(leftConformation.toString(), 0);
            }
            //System.out.println("Starting "+leftConformation+" at "+solutions.get(leftConformation.toString()));
            

            Conformation rightConformation = null;
            int offset = solutions.get(leftConformation.toString());
            //System.out.println("Offset for "+leftConformation+" is "+offset);
            rightConformation = solutionList.get(offset);
            
            solutions.put(leftConformation.toString(), offset+1);

            //System.out.println("Offset for "+leftConformation+" incremented to "+solutions.get(leftConformation.toString()));
            if(solutionList.size() -1 <= offset)
            {
                if(rightSubtree.children.size() < 1){
                    System.out.println("Remove left child "+leftConformation +", offset "+offset+", total solutions "+solutionList.size());
                    //printTree("BEFORE");
                    leftSubtree.getNextConformation(false);
                    //printTree("AFTER");
                }
                if(rightSubtree.children.size() > 0)
                {
                    //System.out.println("Removing a right child...");
                    solutionList.add(rightSubtree.getNextConformation(false));
                }

                rightConformation = solutionList.get(offset);
            }
            
            if(leftSubtree.children.size() > 0)
            {
                BWMAStarNode leftChild = leftSubtree.children.poll();
                
                if(offset + 1 < solutionList.size())
                {
                    Conformation nextConformation = solutionList.get(offset + 1);
                    leftChild.rightSideConformation = nextConformation;
                }
                leftSubtree.children.add(leftChild);
            }

            return leftConformation.join(rightConformation);
        }
        if(children.size()<1) {
            //System.out.println("ITS THE END OF THE WORLD!!!");
            return null;
        }
        if(isLeaf){
            //System.out.println("Process leaf!.");
        	BWMAStarNode out = children.poll();
        	if(reinsert) 
        		children.add(out);
        	//else System.out.println("Remove leaf "+out.partialConformation);
            return out.partialConformation;
        }
        BWMAStarNode next = children.poll();
        if(next.children.size() <= 1 && !next.branching && !next.isLeaf)
        {
            //System.out.println("HOLY WHAT BALLS??!?!?!===");
        }
        //System.out.println("Polled: "+next.partialConformation);
        ////System.out.println("Next conformation is: "+next.partialConformation+", score "+next.partialConformation.score());
        //update next conformation

        Conformation nextConf = next.getNextConformation(reinsert).join(partialConformation);
        if(reinsert || next.children.size() > 0 || (next.branching && next.leftSubtree.children.size() > 0))
        {
            //System.out.println("Reinsert"+next.partialConformation+": "+next.children.size()+" children, branching is "+next.branching);
            if(next.children.size()<=1 && !next.isLeaf && !next.branching)
            {
                //System.out.println("EMPTY CHILD INSERTION");
            }
            children.add(next);
        }
        //else System.out.println("Not reinserting "+next.partialConformation+", it has "+next.children.size()+", we have "+children.size()+" children remaining.");
        //System.out.println("intermediate node complete");
        
        return nextConf;
    }

    
    public void printTree(String prefix, Conformation c)
    {
        Conformation joined = partialConformation;
        if(rightSideConformation != null)
            joined = joined.join(rightSideConformation);
        String output = prefix+joined+", "+peekNextConformation().score();

        System.out.println(output);
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
