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

    /* TODO: REMOVE THIS */
    private BWMAStarNode parent;


    private Conformation rightSideConformation;
    private int totalPossible;
    private int remaining;

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
    
    public void setParent(BWMAStarNode p)
    {
        parent = p;
    }

    public void resort()
    {
        for(BWMAStarNode b : children)
        {
            b.resort();
        }
        if(branching)
        {
            for(BWMAStarNode b : children)
            {
                b.resort();
            }
            for(BWMAStarNode b : rightChildren)
            {
                b.resort();
            }
            PriorityQueue<BWMAStarNode> newLeftHeap = new PriorityQueue<BWMAStarNode>();
            newLeftHeap.addAll(children);
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

    /* TODO: DELETE catch function */
    private boolean allZeros(Conformation c)
    {
    	for(Position p : c.getPositions())
    	{
    		//System.out.println(c.getChoiceAt(p).choice);
    		if(c.getChoiceAt(p).choice != 0)
    			return false;
    	}
    	return true;
    }
    
    private Conformation peekNextConformation () 
    {
        if(allZeros(partialConformation))
        {
        	System.out.println("Catch!");
        }

        if(children.size() < 1)
        {
            //System.out.println("Returning on "+partialConformation+", my right side is "+rightSideConformation);
            return fullConformation(partialConformation);
        }

        if(branching)
        {
            Conformation peeked = children.peek().peekPartial();
            int offset = getRightConformation(peeked);
            Conformation rightSide = null;
            //System.out.println("PEEK "+peeked + " offset is "+offset+", should append "+solutionList.get(offset));
            if(offset >= solutionList.size()){
                if(rightChildren.size() > 0)
                    rightSide = rightChildren.peek().peekNextConformation();
            }
            else rightSide = solutionList.get(offset);
            //System.out.println("Solution list for "+peeked+" : "+solutionList);
            children.peek().rightSideConformation = rightSide;

            if(rightSide == null){
                return fullConformation(partialConformation);
            }
           //System.out.println("Returning "+peeked.join(rightSide));
            return fullConformation(peeked.join(rightSide));
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

        if(children.size() < 1)
        {
            return partialConformation;
        }
        return children.peek().peekNextConformation();
    }

    public boolean moreConformations()
    {
        return children.size() > 0;
    }

    public Conformation getNextConformation () 
    {
        remaining --;
        if(branching){
        	BWMAStarNode leftChild = children.peek();
            Conformation leftConformation = leftChild.peekPartial();
            Conformation rightConformation = null;
            if(rightChildren != null)
            {
                int offset = getRightConformation(leftConformation);
                rightConformation = updateConformationList(leftConformation, offset);
                removeFinishedConformation(offset);
                /* TODO: Recalculation of node position is necessary here. */
            }
            children.add(leftChild);
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
        if(children.size() > 0)
        {
            BWMAStarNode leftChild = children.poll();
            if(offset + 1 < solutionList.size())
            {
                Conformation nextConformation = solutionList.get(offset + 1);
                leftChild.rightSideConformation = nextConformation;
            }
            children.add(leftChild);
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
                BWMAStarNode nextChild = children.poll();
                nextChild.getNextConformation();
                if(nextChild.moreConformations())
                {
                    children.add(nextChild);
                }
            }
            if(rightChildren.size() > 0)
            {
                BWMAStarNode rightChild = rightChildren.poll();
                solutionList.add(rightChild.getNextConformation());
                if(rightChild.moreConformations())
                    rightChildren.add(rightChild);
            }
            rightConformation = solutionList.get(offset);
        }
        return rightConformation;
    }

    public int remainingConformations()
    {
        if(children.size() < 1)
        {
            return remainingRightConformations(partialConformation);
        }
        int sum = 0;
        if(branching)
        {
            for(BWMAStarNode leftChild: children)
            {
                int add = leftChild.remainingConformations();
                sum += add;
            }
        }
        else
        for(BWMAStarNode child: children)
        {
            int add = child.remainingConformations();
            sum += add;
        }
        //printTree("");
        return sum;
    }
    
    public int remainingRightConformations(Conformation c)
    {
        if(!solutions.isEmpty() && !solutions.containsKey(c.toString()))
        {
            if(c.getPositions().containsAll(children.peek().peekPartial().getPositions()))
                getRightConformation(c);
        }
        if(solutions != null && solutions.containsKey(c.toString()))
        {
            return totalRightCombinations() - solutions.get(c.toString());
        }
        if(parent == null)
        {
            return 0;
        }
        return parent.remainingRightConformations(c);

    }

    
    public int totalPossibleCombinations()
    {
        if(children.size() < 1)
        {
            return 1;
        }
        int sum = 0;
        for(BWMAStarNode child : children)
        {
            sum += child.totalPossibleCombinations();
        }
        totalPossible = sum;
        return sum*totalRightCombinations();
    }
    
    private int totalRightCombinations()
    {
        int rightCombinations = 1;
        
        if(branching)
        {
            rightCombinations = 0;           
            for(BWMAStarNode rightChild : rightChildren)
            {
                rightCombinations += rightChild.totalPossibleCombinations();
            }
        }
        return rightCombinations + solutionList.size();
    }
    
    
    public void printTree(String prefix, Conformation c)
    {
        if(prefix.length() < 1)
            System.out.println("BEGIN PRINT TREE==================================================");
        Conformation joined = partialConformation;
        Conformation peeked = peekNextConformation();
        if(rightSideConformation != null)
            joined = joined.join(rightSideConformation);
        String output = prefix+joined+", peek to "+peeked+", "+peeked.score()+
                " totalConformations: "+totalPossibleCombinations()+", rightConformations : "+totalRightCombinations();
        if(parent!= null)
            output += " remaining right conformations "+remainingRightConformations(peekPartial());

        if(parent != null)
            output += ", parent solution list size "+parent.solutionList.size();
        System.out.println(output);
        if(branching)
        {
            for(BWMAStarNode leftChild : children)
            {
                leftChild.printTree(prefix+"L--",joined);
            }
            for(BWMAStarNode rightChild : rightChildren)
            {
                rightChild.printTree(prefix+"R--",joined);
            }
        }
        else
            for(BWMAStarNode child: children)
            {
                child.printTree(prefix+"+--",joined);
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

    public PriorityQueue<BWMAStarNode> getChildren()
    {
        return children;
    }

    public static void CreateTree2(TreeNode root, BWMAStarNode parent, Conformation previous, PriorityQueue<BWMAStarNode> heap, SolutionSpace s, int index, List<? extends Position> positions)
    {
        //System.out.println("Recurse: parent is "+parent.partialConformation+", previous is "+previous+", index "+index);
         if(index < positions.size())
        {
            for(Choice c : s.getChoices(positions.get(index)))
            {
                Conformation nextConf = previous.copy();
                nextConf.append(positions.get(index), c);
                CreateTree2(root, parent, nextConf, heap, s, index + 1, positions);
            }
        }
        else
        {
            BWMAStarNode newNode = new BWMAStarNode(previous);
            newNode.setParent(parent);
            TreeNode nextLeftEdge = searchSubtree(root.getlc());
            TreeNode nextRightEdge = searchSubtree(root.getrc());
           
            if(nextRightEdge != null)
            {
                
                PriorityQueue<BWMAStarNode> nextHeap = newNode.children;
                if(nextLeftEdge != null)
                {
                    newNode.rightChildren = new PriorityQueue<BWMAStarNode>();
                    nextHeap = newNode.rightChildren;
                    newNode.branching = true;
                }
                CreateTree2(nextRightEdge, newNode, newNode.partialConformation, nextHeap, s, 0, nextRightEdge.getCofEdge().getPositionList());
            }
            
            if(nextLeftEdge != null)
            {
                CreateTree2(nextLeftEdge, newNode, newNode.partialConformation, newNode.children, s, 0, nextLeftEdge.getCofEdge().getPositionList());
            }
            
            if(nextLeftEdge != null && nextRightEdge != null)
            {
                /* assign partial conformation */
                BWMAStarNode rightChild = newNode.rightChildren.poll();
                Conformation rightHandSide = rightChild.getNextConformation();
                if(rightChild.moreConformations())
                    newNode.rightChildren.add(rightChild);
                newNode.solutionList.add(rightHandSide);
                for(BWMAStarNode leftChild : newNode.children)
                {
                    leftChild.rightSideConformation = rightHandSide;
                    newNode.getRightConformation(leftChild.peekPartial());
                }
            }

            newNode.remaining = newNode.totalPossibleCombinations();
            if(heap != null)
                heap.add(newNode);
        }
    }
   

    private static TreeNode searchSubtree(TreeNode node) {
        if(node == null) return null;
        boolean isLeaf = node.getIsLeaf();
        boolean isLambdaEdge = node.getCofEdge().getIsLambdaEdge();
        if(isLambdaEdge)
        {
            return node;
        }
        if(isLeaf) 
            return null;
        TreeNode leftChild = node.getlc();
        TreeNode rightChild = node.getrc();
        if(leftChild.getIsLeaf() && !leftChild.getCofEdge().getIsLambdaEdge())
            return searchSubtree(rightChild);
        if(rightChild.getIsLeaf() && !rightChild.getCofEdge().getIsLambdaEdge())
            return searchSubtree(leftChild);
        /* TODO:  incomplete algorithm */
        return null;
    }

}
