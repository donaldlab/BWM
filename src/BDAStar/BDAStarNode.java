package BDAStar;

import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import kstar.TreeNode;

public class BDAStarNode implements Comparable<BDAStarNode> {

    private PriorityQueue<BDAStarNode> children;
    private BDAStarNode leftSubtree;
    private BDAStarNode rightSubtree;
    private Conformation partialConformation;
    boolean isLeaf = false;
    boolean branching = false;
    
    public BDAStarNode(Conformation conf) 
    {
        children = new PriorityQueue<BDAStarNode>();
        partialConformation = conf;
    }
    
    public static BDAStarNode CreateTree(TreeNode root, SolutionSpace s)
    {
        Set<Position> lambda = root.getCofEdge().getPositionSet();
        BDAStarNode AStarRoot = new BDAStarNode(new TestConformation());
        Queue<BDAStarNode> remaining = new LinkedList<BDAStarNode>();
        remaining.add(AStarRoot);
        int cycles = 0;
        while(remaining.size() > 0 && cycles <20)
        {
            cycles ++;
            System.out.println("Remaining: "+remaining.size());
            BDAStarNode current = remaining.poll();
            System.out.println("Remaining after poll: "+remaining.size());
            Conformation currentConf = new TestConformation(current.getConformation());
            System.out.println("Parent conformation positions: "+currentConf.getPositions().size());
            Set<Position> lambdaCopy = difference(lambda, currentConf.getPositions());//copy(lambda);
            /*
            for(Position p : currentConf.getPositions())
            {
            	System.out.println("Position "+p+": "+p.pos);
            }
            lambdaCopy.removeAll(currentConf.getPositions());*/
            boolean isLeaf = lambdaCopy.size() == 1;
            System.out.println("Current conformation positions: "+currentConf.getPositions().size());
            System.out.println("Lambda Size: "+lambdaCopy.size());
           // System.out.println("Leaf? "+isLeaf);
            if(!lambdaCopy.iterator().hasNext()) break;
            Position p = lambdaCopy.iterator().next();
                System.out.println("Processing next position: " + p.pos);
                for (Choice c : s.getChoices(p))
                {
                    System.out.println("Position "+p.pos+", choice "+c.choice);
                    currentConf.append(p, c);
                    //System.out.println("New conformation positions: "+currentConf.getPositions().size());
                    for(Position p2 : currentConf.getPositions())
                    {
                    	//System.out.println("Position "+p2+": "+p2.pos);
                    }
                    
                    BDAStarNode newNode = new BDAStarNode(new TestConformation(currentConf));
                    System.out.println("New node:"+newNode.partialConformation+" score: "+newNode.partialConformation.score());
                    if(isLeaf)
                    {
                        if(root.getIsLeaf())
                        {
                            newNode.isLeaf = true;
                        }
                        else 
                        {
                            System.out.println("BRANCHING!! ");
                            newNode.addBranches(CreateTree(root.getlc(), s), CreateTree(root.getrc(), s));
                            System.out.println("Branch complete.");
                            newNode.branching = true;
                        }
                    }
                    current.addChild(newNode);
                    //System.out.println("New Node conformation size check: "+newNode.getConformation().getPositions().size());
                    remaining.add(newNode);
                    //System.out.println("Next one is "+remaining.peek()+" size: "+remaining.peek().getConformation().getPositions().size());
                    currentConf.delete(p);
                }
            
        }
        return AStarRoot;
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
        	System.out.println("Should we add "+p.pos+"?");
            for(Position p2 : removed)
            {
            	System.out.println("Check against "+p2.pos);
            	if(p.equals(p2)){ skip = true;
            	System.out.println("Skip is true");
            	}
            }
            if(!skip){
            	System.out.println("Adding "+p.pos);
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
        if(branching)
            return partialConformation.join(leftSubtree.peekNextConformation().join(rightSubtree.peekNextConformation()));
        if(isLeaf || /*no children yet...*/ children.size() < 1)
            return partialConformation;
        BDAStarNode peek = children.peek();
        return children.peek().peekNextConformation().join(partialConformation);
    }

    public Conformation getNextConformation () 
    {
        System.out.println("Children size: "+children.size());
        System.out.println("My conformation is: "+partialConformation+", score: "+partialConformation.score());
        
        if(branching){
        	/*
        	 * What actually has to happen here:
        	 * 1. Arbitrarily assign one tree to be the major subtree. Is bigger or smaller better? We need to count conformations.
        	 * 2. Poll the minor tree for the next conformation, and store it in a list of any sort.
        	 * 3. Give all leaves of the major subtree a pointer to the head of that list. From this point forward, all queries traverse the list.
        	 * 4. If the end of the list is hit, poll the minor tree again, and append to the end of the list.
        	 * 5. If the minor tree is depleted, do not reinsert. Otherwise, reinsert.
        	 */
            System.out.println ("Branch!");
            return partialConformation.join(leftSubtree.getNextConformation().join(rightSubtree.getNextConformation()));
        }
        if(isLeaf)
            return partialConformation;
        //if(children.size()<1) return null;
        BDAStarNode next = children.poll();
        //System.out.println("Polled: "+next);
        //System.out.println("Next conformation is: "+next.partialConformation+", score "+next.partialConformation.score());
        Conformation nextConf = next.getNextConformation().join(partialConformation);
        //update next conformation
        if(next.children.size()>0 || next.branching)
        {
            System.out.println("Reinsert"+next.partialConformation+": "+next.children.size()+" children, branching is "+next.branching);
            children.add(next);
        }
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
