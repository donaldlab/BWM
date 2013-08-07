package BDAStar;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
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
public class BDAStarNode implements Comparable<BDAStarNode> {

    private PriorityQueue<BDAStarNode> children;
    private BDAStarNode leftSubtree;
    private BDAStarNode rightSubtree;
    private LinkedList<Conformation> solutionList;
    private Map<String, Integer> solutions;
    private Conformation partialConformation;
    boolean isLeaf = false;
    boolean branching = false;
    
    public BDAStarNode(Conformation conf) 
    {
        children = new PriorityQueue<BDAStarNode>();
        solutionList = new LinkedList<Conformation>();
        solutions = new HashMap<String, Integer>();
        partialConformation = conf;
    }
    
    public static BDAStarNode CreateTree(TreeNode root, Conformation previous, SolutionSpace s)
    {
        Set<Position> lambda = root.getCofEdge().getPositionSet();
        BDAStarNode AStarRoot = new BDAStarNode(new TestConformation(previous));
        System.out.println("Root conformation "+AStarRoot.partialConformation+", score "+AStarRoot.partialConformation.score());
        
        Queue<BDAStarNode> remaining = new LinkedList<BDAStarNode>();
        remaining.add(AStarRoot);
        int cycles = 0;
        while(remaining.size() > 0 && cycles <20)
        {
            cycles ++;
            //System.out.println("Remaining: "+remaining.size());
            BDAStarNode current = remaining.poll();
            //System.out.println("Remaining after poll: "+remaining.size());
            Conformation currentConf = new TestConformation(current.getConformation());
            //System.out.println("Parent conformation positions: "+currentConf.getPositions().size());
            Set<Position> lambdaCopy = difference(lambda, currentConf.getPositions());//copy(lambda);
            /*
            for(Position p : currentConf.getPositions())
            {
            	//System.out.println("Position "+p+": "+p.pos);
            }
            lambdaCopy.removeAll(currentConf.getPositions());*/
            boolean isLeaf = lambdaCopy.size() == 1;
            System.out.println("Current conformation "+currentConf+", score "+currentConf.score());
            //System.out.println("Lambda Size: "+lambdaCopy.size());
           // //System.out.println("Leaf? "+isLeaf);
            if(!lambdaCopy.iterator().hasNext()) break;
            Position p = lambdaCopy.iterator().next();
                //System.out.println("Processing next position: " + p.pos);
                for (Choice c : s.getChoices(p))
                {
                    //System.out.println("Position "+p.pos+", choice "+c.choice);
                    currentConf.append(p, c);
                    ////System.out.println("New conformation positions: "+currentConf.getPositions().size());
                    for(Position p2 : currentConf.getPositions())
                    {
                    	////System.out.println("Position "+p2+": "+p2.pos);
                    }
                    
                    BDAStarNode newNode = new BDAStarNode(new TestConformation(currentConf));
                    //System.out.println("New node:"+newNode.partialConformation+" score: "+newNode.partialConformation.score());
                    if(isLeaf)
                    {
                        if(root.getIsLeaf())
                        {
                            newNode.isLeaf = true;
                        }
                        else 
                        {
                            //System.out.println("BRANCHING!! ");
                            newNode.addBranches(CreateTree(root.getlc(), currentConf, s), CreateTree(root.getrc(), currentConf, s));
                            //System.out.println("Branch complete.");
                            newNode.branching = true;
                        }
                    }
                    current.addChild(newNode);
                    ////System.out.println("New Node conformation size check: "+newNode.getConformation().getPositions().size());
                    remaining.add(newNode);
                    ////System.out.println("Next one is "+remaining.peek()+" size: "+remaining.peek().getConformation().getPositions().size());
                    currentConf.delete(p);
                }
            
        }
        return AStarRoot;
    }
    
    public void resort()
    {
    	for(BDAStarNode b : children)
    	{
    		b.resort();
    	}
    	if(branching)
    	{
    		leftSubtree.resort();
    		rightSubtree.resort();
    	}
    	PriorityQueue<BDAStarNode> newHeap = new PriorityQueue<BDAStarNode>();
    	newHeap.addAll(children);
    	children = newHeap;
    }
    
    private static Set<Position> copy(Set<Position> set)
    {
        Set<Position> newset = new LinkedHashSet<Position>();
        for(Position p : set)
        {
            newset.add(p);
        }
        return newset;
    }
    
    private static Set<Position> difference(Set<Position> set, Collection<Position> removed)
    {
        Set<Position> newset = new LinkedHashSet<Position>();
        for(Position p : set)
        {
        	boolean skip = false;
        	//System.out.println("Should we add "+p.pos+"?");
            for(Position p2 : removed)
            {
            	//System.out.println("Check against "+p2.pos);
            	if(p.equals(p2)){ skip = true;
            	//System.out.println("Skip is true");
            	}
            }
            if(!skip){
            	//System.out.println("Adding "+p.pos);
            	newset.add(p);
            }
        }

        return newset;
    }
    
    private Conformation getConformation () {
        return partialConformation;
    }

    private void addChild(BDAStarNode node)
    {
        children.add(node);
    }
    
    private void addBranches(BDAStarNode leftChild, BDAStarNode rightChild)
    {
        leftSubtree = leftChild;
        rightSubtree = rightChild;
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
    	//System.out.println(partialConformation+" called for heuristic...");
        if(branching){
        	Conformation peeked = leftSubtree.peekNextConformation().join(partialConformation);
        	if(peeked == null){
        		System.out.println("No left sub tree");
        		return partialConformation;
        	}
        	if(!solutions.containsKey(peeked.toString())) {
        		//System.out.println("Returning default optima"+ peeked.join(rightSubtree.peekNextConformation()));
        		return peeked.join(rightSubtree.peekNextConformation()).join(partialConformation);
        	}
        	int offset = solutions.get(peeked.toString());
        	Conformation rightSide = null;
        	if(offset >= solutionList.size()){
        		rightSide = rightSubtree.peekNextConformation();
        		//System.out.println("Getting next best optima:" + partialConformation.join(peeked.join(rightSide)));
        	}
        	else rightSide = solutionList.get(offset);
        	if(rightSide == null && rightSubtree.moreConformations()) 
        		rightSide = rightSubtree.peekNextConformation();
        	else {
        		leftSubtree.getNextConformation();
        		return peekNextConformation();
        	}
        	if(rightSide == null){
        		System.out.println("No right subtree");
        		return partialConformation;
        	}
            return partialConformation.join(peeked.join(rightSide));
        }
        if(isLeaf) 
        	return partialConformation;
        if(children.size() < 1)
        {
        	return partialConformation;
        }
        BDAStarNode peek = children.peek();
        //System.out.println("Return child score: "+children.peek().peekNextConformation().join(partialConformation));
        return children.peek().peekNextConformation().join(partialConformation);
    }
    
    public boolean moreConformations()
    {
    	if(isLeaf) return true;
    	if(branching) return rightSubtree.moreConformations() || 
    			(solutionList.size() > solutions.get(leftSubtree.peekNextConformation().toString()));
    	boolean more = false;
    	for(BDAStarNode b : children)
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
            BDAStarNode leftChild = leftSubtree;
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
            
            
            
            Conformation rightChild = null;
            if(pointer.hasNext()){
                rightChild = pointer.next();
                //leftSubtree.children.add(leftChild);
            }
            if(rightChild == null)
            {
                //System.out.println("polling for new right conformation...");
                rightChild = rightSubtree.getNextConformation(false);
                if(rightChild != null){ 
                    //put child back in, it has more to go!
                    //System.out.println("Reinsert "+leftChild.partialConformation);
                    solutionList.add(rightChild);
                }
                else // no more conformations
                {
                    //System.out.println("Right subtree is depleted, removing "+leftConformation);
                    branching = leftSubtree.children.size()>0;
                    leftSubtree.getNextConformation();
                    solutions.put(leftConformation.toString(), (solutions.get(leftConformation.toString())+1));
                    
                    return getNextConformation(reinsert);
                }
            }
            solutions.put(leftConformation.toString(), (solutions.get(leftConformation.toString())+1));
            
            return partialConformation.join(leftConformation.join(rightChild));
        }
        if(isLeaf){
            //System.out.println("Process leaf!.");
            return partialConformation;
        }
        if(children.size()<1) {
            //System.out.println("ITS THE END OF THE WORLD!!!");
            return null;
        }
        //System.out.println("Process intermediate node...");
        BDAStarNode next = children.poll();
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
        for(BDAStarNode child: children)
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
        for(BDAStarNode b: children)
        {
            sum += b.remainingConformations();
        }
        return sum;
    }
    
    public void printTree(String prefix)
    {
        printTree(prefix, partialConformation);
    }
    
    public int compareTo(BDAStarNode node)
    {
        if(nextBestScore() - node.nextBestScore() < 0) return -1;
        if(nextBestScore() - node.nextBestScore() > 0) return 1;
        return 0;
    }
   
}
