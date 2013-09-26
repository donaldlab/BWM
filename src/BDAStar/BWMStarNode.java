package BDAStar;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
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
public class BWMStarNode 
{

    /* TODO: REMOVE THIS */
    private BWMStarNode parent;

    private ProteinConformationTrie conformationTrees;
    private LinkedList<Conformation> solutionList;
    private Map<String, Integer> solutions;
    private BWMStarNode child;
    private BWMStarNode rightChild;

    private Conformation rightSideConformation;
    private int totalPossible;
    private int remaining;

    public BWMStarNode() 
    {
        solutionList = new LinkedList<Conformation>();
        solutions = new HashMap<String, Integer>();
    }


    
    public void setParent(BWMStarNode p)
    {
        parent = p;
    }

    public void resort()
    {
        child.resort();
        rightChild.resort();
        conformationTrees.resort();
        
    }

    public double nextBestScore(Conformation partialConformation)
    {
        /* TODO: Do we want this or do we want a cheaper peek? 
         * We could probably do it just from the trie's score. */
        return peekNextConformation(partialConformation).score();
    }

    private int getRightConformationOffset(Conformation leftConformation)
    {
        if(!solutions.containsKey(leftConformation.toString()))
            solutions.put(leftConformation.toString(), 0);
        return solutions.get(leftConformation.toString());
    }
    
    private Conformation peekNextConformation (Conformation partial) 
    {
        if(child == null)
        {
            return conformationTrees.getAStarRoot(partial, 0).peekNextConformation();
        }
        Conformation soFar = conformationTrees.getAStarRoot(partial, 0).peekNextConformation();
        Conformation conformation = child.peekNextConformation(soFar);
        Conformation rightConformation = getRightConformation(conformation);
        return conformation.join(rightConformation);
    }
    
    private Conformation getRightConformation (Conformation conformation) {
        if(!solutions.containsKey(conformation.toString()))
        {
            solutions.put(conformation.toString(), 0);
        }
        return solutionList.get(solutions.get(conformation.toString()));
    }



    private double peekNextScore()
    {
        return conformationTrees.getAStarRoot(null, 0).peekNextConformation().score();
    }

    public Conformation getNextConformation (Conformation partial) 
    {
        //If we don't actually have a lambda set, recurse
        
        BDAStarNode AStarRoot = conformationTrees.getAStarRoot(partial, 0);
        Conformation rootPartial = AStarRoot.getNextConformation();
        // If we're a leaf, return whatever partialConformation fits.
        if(child == null)
        {
            return rootPartial;
        }
        //If we have one child, return the joint result.
        Conformation jointPartial = rootPartial.join(partial);
        Conformation childPartial = child.getNextConformation(jointPartial);
        if(child.moreConformations(rootPartial) && rightChild != null)
        {
            AStarRoot.insertConformation(rootPartial);
        }
        Conformation out = jointPartial.join(childPartial);
        if(rightChild == null)
            return out;
        //If we have a right child as well, get the right conformation and join it.
        Conformation rightConformation = getRightConformation(out);
        updateRightConformations(partial, rootPartial, out, AStarRoot);
        return out.join(rightConformation);
    }
    
    
    private boolean moreConformations (Conformation out) 
    {
        return conformationTrees.getAStarRoot(out, 0).moreConformations();
    }



    private void insertConformation (Conformation c) {
        //remove this node's lambda set
        //insert the part that applies to this trie
        conformationTrees.insertConformation(c);
        //pass on the rest
        if(child != null)
            child.insertConformation(c);
        if(rightChild != null)
            rightChild.insertConformation(c);
    }



    private void updateRightConformations (Conformation MSet, Conformation lambda, Conformation lastUsed, BDAStarNode lastNode) {
        int offset = getRightConformationOffset(lastUsed);
        if(offset + 1 == solutionList.size())
        {
            BDAStarNode rightTree = rightChild.conformationTrees.getAStarRoot(MSet, 0);
            if(rightTree.moreConformations())
            {
                Conformation nextRightConf = rightTree.getNextConformation();
                solutionList.add(nextRightConf);
            }
            else 
            {
                conformationTrees.getAStarRoot(MSet, 0).deleteConformation(lastUsed);
            }
        }
    }



    public static BWMStarNode CreateTree(TreeNode node, BWMSolutionSpace space)
    {
        BWMStarNode root = new BWMStarNode();;
        if(node.getCofEdge().getIsLambdaEdge())
        {
            createConformationTrie(node, root, space);
        }
        TreeNode leftChild = searchSubtree(node.getlc());
        TreeNode rightChild = searchSubtree(node.getrc());
        if(leftChild != null)
        {
            root.child = CreateTree(leftChild, space);
        }
        if(rightChild != null)
        {
            BWMStarNode newChild = CreateTree(rightChild, space);
            if(leftChild != null)
            {
                root.rightChild = newChild;
            }
            else root.child = newChild;
        }
        return root;
    }
    
    private static void createConformationTrie(TreeNode node, BWMStarNode root, BWMSolutionSpace space)
    {
        LinkedHashSet<Integer> MSet = node.getCofEdge().getM();
        LinkedHashSet<Integer> lambda = node.getCofEdge().getLambda();
        Integer[] dummy = new Integer[]{};
        Position[] MArray = space.MSetFromArray(MSet, null).toArray(new Position[]{});
        root.conformationTrees = ProteinConformationTrie.createTrie(space, MArray, 0);
        root.generateConformations(MSet.toArray(dummy), lambda.toArray(dummy), 0, space, space.getEmptyConformation());
        
    }
    
    private void generateConformations(Integer[] MSet, Integer[] lambda, int index, BWMSolutionSpace space, Conformation current)
    {
        if(index == MSet.length + lambda.length - 2)
        {
            conformationTrees.insertConformation(current);
        }
        if(index >= MSet.length)
        {
            Position p = space.positionFromPos(lambda[index]);
            for(Choice c : space.getChoices(p))
            {
                Conformation nextConf = current.copy();
                nextConf.append(p, c);
                generateConformations(MSet, lambda, index + 1, space, nextConf);
            }
            
        }
        if(index < MSet.length)
        {
            Position p = space.positionFromPos(MSet[index]);            
            for(Choice c : space.getChoices(p))
            {
                Conformation nextConf = current.copy();
                nextConf.append(p, c);
                generateConformations(MSet, lambda, index + 1, space, nextConf);
            }
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
        /* TODO:  incomplete algorithm, there are other cases */
        return null;
    }

}
