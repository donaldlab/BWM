package BDAStar;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
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
    boolean isSubtreeRoot = false;
    boolean branching = false;


    private Conformation rightSideConformation;

    public BWMAStarNode(Conformation conf) 
    {
        children = new PriorityQueue<BWMAStarNode>();
        solutionList = new LinkedList<Conformation>();
        solutions = new HashMap<String, Integer>();
        partialConformation = conf;
    }

    public void addChild(Conformation conf)
    {
        children.add(new BWMAStarNode(conf));
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

    private void addBranches(BWMAStarNode leftChild, BWMAStarNode rightChild)
    {
        leftSubtree = leftChild;
        rightSubtree = rightChild;
        branching = true;
    }

    public double nextBestScore()
    {
        return peekNextConformation().score();
    }

    private int getRightConformation(Conformation leftConformation)
    {
        if(!solutions.containsKey(leftConformation.toString()))
            solutions.put(leftConformation.toString(), 0);
        return solutions.get(leftConformation.toString());
    }

    private Conformation peekNextConformation () 
    {

        if(isSubtreeRoot) 
        {
            return fullConformation(children.peek().peekNextConformation());
        }

        if(branching)
        {
            Conformation peeked = leftSubtree.peekPartial();
            int offset = getRightConformation(peeked);
            Conformation rightSide = null;
            if(offset >= solutionList.size()){
                rightSide = rightSubtree.peekNextConformation();
            }
            else rightSide = solutionList.get(offset);

            if(rightSide == null){
                return fullConformation(partialConformation);
            }
            return peeked.join(rightSide);
        }
        
        if(children.size() < 1)
        {
            return fullConformation(partialConformation);
        }


        return fullConformation(children.peek().peekNextConformation());
    }

    private Conformation fullConformation(Conformation out)
    {
        if(rightSideConformation != null)
        {
            out = out.join(rightSideConformation);
        }
        return out;
    }


    public Conformation peekPartial()
    {
        if(isSubtreeRoot) 
        {
            return children.peek().peekPartial();
        }

        if(branching)
        {
            return leftSubtree.peekNextConformation();
        }

        if(children.size() < 1)
        {
            return partialConformation;
        }

        return children.peek().peekNextConformation();
    }

    public boolean moreConformations()
    {
        return children.size() > 0 || (branching && leftSubtree.moreConformations());
    }

    public Conformation getNextConformation () 
    {
    	if(isSubtreeRoot){
            return children.poll().partialConformation;
        }

        if(branching){
            Conformation leftConformation = leftSubtree.peekPartial();
            Conformation rightConformation = null;

            int offset = getRightConformation(leftConformation);
            rightConformation = updateConformationList(leftConformation, offset);
            removeFinishedConformation(offset);

            return leftConformation.join(rightConformation);
        }


        BWMAStarNode next = children.poll();
        Conformation nextConf = next.getNextConformation();
        if(next.moreConformations())
        {
            children.add(next);
        }

        return nextConf;
    }


    private void removeFinishedConformation (int offset) 
    {
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
    }


    private Conformation updateConformationList (Conformation leftConformation, int offset) 
    {
        Conformation rightConformation;
        rightConformation = solutionList.get(offset);
        solutions.put(leftConformation.toString(), offset+1);

        if(solutionList.size() -1 <= offset)
        {
            if(rightSubtree.children.size() < 1){
                leftSubtree.getNextConformation();
            }
            if(rightSubtree.children.size() > 0)
            {
                solutionList.add(rightSubtree.getNextConformation());
            }
            rightConformation = solutionList.get(offset);
        }
        return rightConformation;
    }


    public void printTree(String prefix, Conformation c)
    {
        Conformation joined = partialConformation;
        if(rightSideConformation != null)
            joined = joined.join(rightSideConformation);
        String output = prefix+joined+", "+peekNextConformation().score()+" leaf? "+isSubtreeRoot;

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
    
    public void insertChild(Set<Position> MSet, BWMAStarNode node)
    {
        if(partialConformation.getPositions().containsAll(MSet) || MSet.isEmpty())
            children.add(node);
        else
        {
        	HashSet<Position> nextSet = new HashSet<Position>();
        	nextSet.addAll(MSet);
        	MSet.removeAll(partialConformation.getPositions());
            if(shareElements(leftSubtree.partialConformation.getPositions(), MSet))
            {
                leftSubtree.insertChild(MSet, node);
            }
            else 
                rightSubtree.insertChild(MSet, node);
        }
    }
    
    private boolean shareElements(Iterable<Position> set1, Set<Position> set2)
    {
        for(Position p1: set1)
            if(set2.contains(p1))
                return true;
        return false;   
    }
    

    public static BWMAStarNode CreateTree(TreeNode root, Conformation previous, SolutionSpace s)
    {
        List<Position> lambda = root.getCofEdge().getPositionList();
        BWMAStarNode AStarRoot = new BWMAStarNode(previous);
        populateHeap(AStarRoot.children, lambda, 0, previous, s);
        if(!root.getIsLeaf())
        {
            for(BWMAStarNode child : AStarRoot.children)
            {
                child.addBranches(CreateTree(root.getlc(), child.partialConformation, s), CreateTree(root.getrc(), child.partialConformation, s));
                Conformation rightHandSide = child.rightSubtree.getNextConformation();
                child.solutionList.add(rightHandSide);
                for(BWMAStarNode leftChild : child.leftSubtree.children)
                {
                    leftChild.rightSideConformation = rightHandSide;
                }
            }
        }
        else {
        	AStarRoot.isSubtreeRoot = true;
        }
        return AStarRoot;
    }
    
    public static void populateHeap(PriorityQueue<BWMAStarNode> heap, List<Position> positions, int index, Conformation currentConf, SolutionSpace s)
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
    
}
