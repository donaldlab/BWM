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

TODO:
Compress leftSubtree into leftChildren, rightSubtree into rightChildren

GetNextConformation function is different.

AbstractBDAStarNode
ListingBDAStarNode
TopKBDAStarNode
ExpandingBDAStarNode
 */
public class BWMAStarNode implements Comparable<BWMAStarNode> {

    private PriorityQueue<BWMAStarNode> children;
    private PriorityQueue<BWMAStarNode> leftChildren;
    private PriorityQueue<BWMAStarNode> rightChildren;
    private LinkedList<Conformation> solutionList;
    private Map<String, BWMAStarNode> childrenMap;
    private Map<String, Integer> solutions;
    private Conformation partialConformation;
    boolean isSubtreeRoot = false;
    boolean branching = false;


    private Conformation rightSideConformation;

    public BWMAStarNode(Conformation conf) 
    {
    	childrenMap = new HashMap<String, BWMAStarNode>();
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
            for(BWMAStarNode b : leftChildren)
            {
                b.resort();
            }
            for(BWMAStarNode b : rightChildren)
            {
                b.resort();
            }
            PriorityQueue<BWMAStarNode> newLeftHeap = new PriorityQueue<BWMAStarNode>();
            newLeftHeap.addAll(leftChildren);
            PriorityQueue<BWMAStarNode> newRightHeap = new PriorityQueue<BWMAStarNode>();
            newRightHeap.addAll(rightChildren);
        }
        PriorityQueue<BWMAStarNode> newHeap = new PriorityQueue<BWMAStarNode>();
        newHeap.addAll(children);
        children = newHeap;
    }

    private void addBranches(PriorityQueue<BWMAStarNode> leftChild,  PriorityQueue<BWMAStarNode> rightChild)
    {
        leftChildren = leftChild;
        rightChildren = rightChild;
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
            Conformation peeked = leftChildren.peek().partialConformation;
            int offset = getRightConformation(peeked);
            Conformation rightSide = null;
            if(offset >= solutionList.size()){
                rightSide = rightChildren.peek().peekNextConformation();
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
            return leftChildren.peek().peekNextConformation();
        }

        if(children.size() < 1)
        {
            return partialConformation;
        }

        return children.peek().peekNextConformation();
    }

    public boolean moreConformations()
    {
        return children.size() > 0 || (branching && leftChildren.size() > 0);
    }

    public Conformation getNextConformation () 
    {
        if(branching){
            Conformation leftConformation = leftChildren.peek().partialConformation;
            Conformation rightConformation = null;

            int offset = getRightConformation(leftConformation);
            rightConformation = updateConformationList(leftConformation, offset);
            removeFinishedConformation(offset);

            return leftConformation.join(rightConformation);
        }
        
        if(children.size() < 1){
            return fullConformation(partialConformation);
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
        if(leftChildren.size() > 0)
        {
            BWMAStarNode leftChild = leftChildren.poll();
            if(offset + 1 < solutionList.size())
            {
                Conformation nextConformation = solutionList.get(offset + 1);
                leftChild.rightSideConformation = nextConformation;
            }
            leftChildren.add(leftChild);
        }
    }


    private Conformation updateConformationList (Conformation leftConformation, int offset) 
    {
        Conformation rightConformation;
        rightConformation = solutionList.get(offset);
        solutions.put(leftConformation.toString(), offset+1);

        if(solutionList.size() -1 <= offset)
        {
            if(rightChildren.size() < 1){
                leftChildren.poll();
            }
            if(rightChildren.size() > 0)
            {
                solutionList.add(rightChildren.poll().getNextConformation());
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
            for(BWMAStarNode leftChild : leftChildren)
            {
                leftChild.printTree(prefix+"L--",joined);
            }
            for(BWMAStarNode rightChild : rightChildren)
            {
                rightChild.printTree(prefix+"R--",joined);
            }
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
    
    public void insertChild(Set<ProteinPosition> MSet, BWMAStarNode node)
    {
        if(MSet.isEmpty())
            children.add(node);
        else
        {
        	HashSet<ProteinPosition> nextSet = new HashSet<ProteinPosition>();
        	nextSet.addAll(MSet);
        	nextSet.removeAll(partialConformation.getPositions());
            childrenMap.get(MSet).insertChild(nextSet, node);
        }
    }
    
    private boolean shareElements(Iterable<Position> set1, Set<Position> set2)
    {
        for(Position p1: set1)
            if(set2.contains(p1))
                return true;
        return false;   
    }
    
    private void addChildMap(BWMAStarNode child)
    {
    	childrenMap.put(child.partialConformation.toString(), child);
    }
    
    

    public static BWMAStarNode CreateTree(TreeNode root, Conformation previous, SolutionSpace s, int index)
    {
        List<? extends Position> lambda = root.getCofEdge().getPositionList();
        BWMAStarNode AStarRoot = new BWMAStarNode(previous);
        if(index < lambda.size())
            populateHeap(root, AStarRoot, AStarRoot.children, lambda, index, previous, s);
        if(!root.getIsLeaf())
        {
            for(BWMAStarNode child : AStarRoot.children)
            {
            	child.leftChildren = new  PriorityQueue<BWMAStarNode>();
            	child.rightChildren = new  PriorityQueue<BWMAStarNode>();
                populateHeap(root.getlc(), child, child.leftChildren, root.getlc().getCofEdge().getPositionList(), 0, child.partialConformation, s);
                populateHeap(root.getrc(), child, child.rightChildren, root.getrc().getCofEdge().getPositionList(), 0, child.partialConformation, s);
            	Conformation rightHandSide = child.rightChildren.poll().partialConformation;
            	child.branching = true;
                child.solutionList.add(rightHandSide);
                for(BWMAStarNode leftChild : child.leftChildren)
                {
                    leftChild.rightSideConformation = rightHandSide;
                }
            }
        }
        return AStarRoot;
    }
    
    public static void populateHeap(TreeNode root, BWMAStarNode node, PriorityQueue<BWMAStarNode> heap, List<? extends Position> positions, int index, Conformation currentConf, SolutionSpace s)
    {
    	if(index == 0 && positions.size() < 1) return;
        for(Choice c : s.getChoices(positions.get(index)))
        {
            Conformation nextConf = currentConf.copy();
            nextConf.append(positions.get(index), c);
            System.out.println("Creating conformation " + nextConf);
            if(index >= positions.size() - 1)
            {
            	BWMAStarNode next = CreateTree(root, nextConf, s, index+1);
            	node.addChildMap(next);
            	heap.add(next);
            }
            else 
                populateHeap(root, node, heap, positions, index+1, nextConf, s);
        }
    }
    
}
