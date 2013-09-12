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
                if(rightChildren.size() > 0)
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
            if(rightChildren != null)
            {
	            int offset = getRightConformation(leftConformation);
	            rightConformation = updateConformationList(leftConformation, offset);
	            removeFinishedConformation(offset);
	        }

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
        String output = prefix+joined+", "+peekNextConformation().score()+
        		" children: "+children+children.hashCode();
        if(branching)
        	output+= ", left : "+leftChildren+leftChildren.hashCode()+", right "+
        		rightChildren+rightChildren.hashCode();

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
    
    private String conformationStringFromMSet(Set<ProteinPosition> MSet)
    {
        return "";
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
    
    public static BWMAStarNode CreateTree2(TreeNode root, Conformation previous, SolutionSpace s, int index, PriorityQueue<BWMAStarNode> heap,
    		List<? extends Position> positions)
    {
    	/*
    	 * 1. If not complete lambda set, recurse, incrementing index and appending to the conformation.
    	 * 2. (tree compacting section) If:
    	 *      a) the current node has one leaf child with no lambda set, merge it into the current node and don't recurse.
    	 *      b) the current node (now) has no lambda set and has only one child, merge it into its parent and delete it from the structure.
    	 * 3. Create new tree node with the current conformation.
    	 * 4. Add it to the provided heap.
    	 * 5. find next treenode by applying recursive compacting algorithm on whatever tree exists, preferentially left
    	 * 6. run the compacting search on the right treenode
    	 * 7. if both left and right exist, recurse on the left, passing in the first heap, and then recurse on the right, passing in the right heap,
    	 *     then assign the partial conformation for the right subtree to the left subtree.
    	 * 8. if only one exist, recurse on that child.
    	 * 
    	 */
    	return null;
    }
    

    public static BWMAStarNode CreateTree(TreeNode root, Conformation previous, SolutionSpace s, int index)
    {
        /*
        System.out.println("PRINTING TREE");
        root.printTree("");
        System.out.println("END PRINT TREE");
        */
        List<? extends Position> lambda = root.getCofEdge().getPositionList();
        System.out.println(root+", Index "+index+", "+previous+", "+lambda.size());
        BWMAStarNode AStarRoot = new BWMAStarNode(previous);
        if(index < lambda.size())
            populateHeap(root, AStarRoot, AStarRoot.children, lambda, index, previous, s);
        System.out.println("Population of "+root+" complete. Resulting heap size: "+AStarRoot.children.size());
        if(!root.getIsLeaf())
        {
            System.out.println("Processing left and right subtrees...."+AStarRoot.children.size());
            if(AStarRoot.children.size() < 1)
                AStarRoot.children.add(new BWMAStarNode(previous));
            for(BWMAStarNode child : AStarRoot.children)
            {
                
                //TODO: Update algorithm to account for missing left or right subtrees.
            	child.branching = true;
                if(root.getlc() != null)
                {
                    System.out.println("Left Child of "+root);
                    TreeNode leftTreeNode = root.getlc();
                    System.out.println("is: "+leftTreeNode);
                    child.leftChildren = new  PriorityQueue<BWMAStarNode>();
                    populateHeap(leftTreeNode, child, child.leftChildren, leftTreeNode.getCofEdge().getPositionList(), 0, child.partialConformation, s);

                }
                if(root.getrc() != null)
                {
                    System.out.println("Right Child of "+root);
                    System.out.println(root.getrc());
                    child.rightChildren = new  PriorityQueue<BWMAStarNode>();
                    populateHeap(root.getrc(), child, child.rightChildren, root.getrc().getCofEdge().getPositionList(), 0, child.partialConformation, s);
                    
                    if(child.rightChildren.size() > 0)
                    {
                        Conformation rightHandSide = child.rightChildren.peek().partialConformation;
                        child.branching = true;
                        child.solutionList.add(rightHandSide);
                        for(BWMAStarNode leftChild : child.leftChildren)
                        {
                            leftChild.rightSideConformation = rightHandSide;
                        }
                    }
                    
                }
            }
        }
        else System.out.println(root+" has no children.");
        return AStarRoot;
    }
    
    public static void populateHeap(TreeNode root, BWMAStarNode node, PriorityQueue<BWMAStarNode> heap, List<? extends Position> positions, int index, Conformation currentConf, SolutionSpace s)
    {
        System.out.println("Populating "+root+", index: "+index +", number of positions: "+positions.size());
    	if(index == 0 && positions.size() < 1)
    	{
    	    System.out.println("Emptynode: "+index+" "+positions.size()+" leaf: "+root.getIsLeaf());
    	    System.out.println("Parent: "+root.getp());
            System.out.println("Current: "+root);
            //root.getp().printTree("");
    	    System.out.println("BLANK");
        
            BWMAStarNode next = CreateTree(root, currentConf, s, index + 1);
            System.out.println("Adding "+next+" to "+heap+heap.hashCode());
            node.addChildMap(next);
            heap.add(next);
            return;
    	}
    	System.out.println("Creating new conformation... number of choices: "+s.getChoices(positions.get(index)).size());
           
        for(Choice c : s.getChoices(positions.get(index)))
        {
             Conformation nextConf = currentConf.copy();
            nextConf.append(positions.get(index), c);
            System.out.println("Creating conformation " + nextConf);
            if(index >= positions.size() - 1)
            {
            	BWMAStarNode next = CreateTree(root, nextConf, s, index+1);
            	node.addChildMap(next);

                System.out.println("Adding "+next+" to "+heap+heap.hashCode());
            	heap.add(next);
            	System.out.println("Heap is now: "+heap.size());
            	if(heap.size() == 3)
               	{
            	    System.out.println("Catch!");
            	}
            }
            else 
                populateHeap(root, node, heap, positions, index+1, nextConf, s);
        }
    }
    
}
