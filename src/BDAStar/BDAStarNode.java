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
        while(remaining.size() > 0 && cycles <15)
        {
            cycles ++;
            System.out.println("Remaining: "+remaining.size());
            BDAStarNode current = remaining.poll();
            Conformation currentConf = new TestConformation(current.getConformation());
            Set<Position> lambdaCopy = copy(lambda);
            lambdaCopy.removeAll(currentConf.getPositions());
            boolean isLeaf = lambdaCopy.size() == 0;
            System.out.println("Total conformation positions: "+currentConf.getPositions().size());
            System.out.println("Lambda Size: "+lambdaCopy.size());
            Position p = lambdaCopy.iterator().next();
                System.out.println("WHEEE" + p.pos);
                for (Choice c : s.getChoices(p))
                {
                    System.out.println("Position "+p.pos+", choice "+c.choice);
                    currentConf.append(c);
                    System.out.println("Current conformation positions: "+currentConf.getPositions().size());
                    
                    BDAStarNode newNode = new BDAStarNode(currentConf);
                    if(isLeaf)
                    {
                        if(root.getIsLeaf())
                        {
                            newNode.isLeaf = true;
                        }
                        else 
                        {
                            newNode.addBranches(CreateTree(root.getlc(), s), CreateTree(root.getrc(), s));
                            newNode.branching = true;
                        }
                    }
                    current.addChild(newNode);
                    remaining.add(current);
                    currentConf.deleteLast();
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
            return leftSubtree.peekNextConformation().join(rightSubtree.peekNextConformation());
        if(isLeaf || /*no children yet...*/ children.size() < 1)
            return partialConformation;
        return children.peek().peekNextConformation();
    }

    public Conformation getNextConformation () {
        if(branching)
            return leftSubtree.getNextConformation().join(rightSubtree.getNextConformation());
        if(isLeaf)
            return partialConformation;
        BDAStarNode next = children.poll();
        Conformation nextConf = next.getNextConformation();
        //update next conformation
        children.add(next);
        return nextConf;
    }
    
    public int compareTo(BDAStarNode node)
    {
        if(nextBestScore() - node.nextBestScore() < 0) return -1;
        if(nextBestScore() - node.nextBestScore() > 0) return 1;
        return 0;
    }
   
}
