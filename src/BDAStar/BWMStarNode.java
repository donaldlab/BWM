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

    private Set<? extends Position> MSet;
    private Set<? extends Position> lambdaSet;

    private Conformation rightSideConformation;
    private Conformation emptyConformation;
    private int totalPossible;
    private int remaining;

    public BWMStarNode(Set<? extends Position> M, Set<? extends Position> lambda, 
            Conformation empty)
    {
        solutionList = new LinkedList<Conformation>();
        solutions = new HashMap<String, Integer>();
        MSet = M;
        lambdaSet = lambda;
        emptyConformation = empty;
    }



    public void setParent(BWMStarNode p)
    {
        parent = p;
    }

    public void resort()
    {
        if(child != null)
            child.resort();
        if(rightChild != null)
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
            return conformationTrees.getAStarRoot(partial).peekNextConformation();
        }
        Conformation MConf = partial.extract(MSet);
        Conformation soFar = conformationTrees.getAStarRoot(MConf).peekNextConformation();
        Conformation conformation = child.peekNextConformation(soFar);
        Conformation rightConformation = getRightConformation(conformation);
        return conformation.join(rightConformation);
    }

    private double peekNextScore()
    {
        return conformationTrees.getAStarRoot(null).peekNextConformation().score();
    }

    public Conformation getNextConformation ()
    {
        return getNextConformation(emptyConformation);
    }

    public Conformation getNextConformation (Conformation parentConformation) 
    {
        /*
         * Algorithm overview:
         * 1. Get AStarRoot from MConf (can be empty conformation)
         * 2. Get the partial conformation (lambdaConf) from AStarRoot
         * 3. If we're a leaf, return lambdaConf.
         * 4. If we have a child, add MConf to lambdaConf (childMSet) and recurse on child.
         * 5. If the child believes childMSet has more conformations, reinsert lambdaConf.
         */
        Conformation MConformation = parentConformation.extract(MSet);
        BDAStarNode AStarRoot = conformationTrees.getAStarRoot(MConformation);
        Conformation lambdaConformation = AStarRoot.getNextConformation();
        // If we're a leaf, return whatever partialConformation fits.
        if(child == null)
        {
            return lambdaConformation;
        }
        //If we have one child, return the joint result.
        Conformation childMConformation = lambdaConformation.join(parentConformation);
        Conformation childNextConformation = child.getNextConformation(childMConformation);
        if(child.moreConformations(childMConformation))
        {
            System.out.println("Reinserting "+childMConformation);
            AStarRoot.insertConformation(lambdaConformation, parentConformation);
        }
        Conformation out = childMConformation.join(childNextConformation);
        if(rightChild == null)
            return out;
        //If we have a right child as well, get the right conformation and join it.
        Conformation rightConformation = getRightConformation(out);
        //updateRightConformations(parentConformation, lambdaConformation, out);
        return out.join(rightConformation);
    }


    private boolean moreConformations (Conformation out) 
    {
        return conformationTrees.getAStarRoot(out).moreConformations();
    }



    private void insertConformation (Conformation c) {
        //remove this node's lambda set
        Conformation MConformation = c.extract(MSet);
        Conformation lambdaConformation = c.extract(lambdaSet);
        //insert the part that applies to this trie
        conformationTrees.insertConformation(MConformation, lambdaConformation);
        //pass on the rest
        if(child != null)
            child.insertConformation(c);
        if(rightChild != null)
            rightChild.insertConformation(c);
    }

    private Conformation getRightConformation (Conformation conformation) {
        int offset = getRightConformationOffset(conformation);

        if(offset < solutionList.size())
            return solutionList.get(offset);
        Conformation MConf = conformation.extract(MSet);
        Conformation lambdaConf = conformation.extract(lambdaSet);
        updateRightConformations(MConf, lambdaConf, conformation);
        solutions.put(conformation.toString(), offset+1);

        return solutionList.get(offset);
    }

    private void updateRightConformations (Conformation MSet, Conformation lambda, Conformation lastUsed) {
        int offset = getRightConformationOffset(lastUsed);
        if(offset >= solutionList.size())
        {
            BDAStarNode rightTree = rightChild.conformationTrees.getAStarRoot(MSet);
            if(rightTree.moreConformations())
            {
                rightTree.printTree();
                
                Conformation nextRightConf = rightTree.getNextConformation();
                solutionList.add(nextRightConf);
                boolean more = rightTree.moreConformations(); 
                System.out.println("Reinserting1 "+lastUsed+" index is now "+offset+", are there more conformations? "+more);
                if(rightTree.moreConformations())
                {
                    child.insertConformation(lastUsed);
                }
                else 
                {
                    System.out.println(lastUsed + " depleted at index" + (offset+1));
                    conformationTrees.getAStarRoot(MSet).deleteConformation(lastUsed);
                }
            }
            else 
            {
                System.out.println(lastUsed + " depleted at index" + offset);
                conformationTrees.getAStarRoot(MSet).deleteConformation(lastUsed);
            }
        }
        else 
        {
            System.out.println("Reinserting2 "+lastUsed);
            child.insertConformation(lastUsed);
        }
    }



    public static BWMStarNode CreateTree(TreeNode node, SolutionSpace space)
    {
        Set<? extends Position> MSet = space.MSetFromArray(node.getCofEdge().getM());
        Set<? extends Position> lambdaSet = space.MSetFromArray(node.getCofEdge().getLambda());
        BWMStarNode root = new BWMStarNode(MSet, lambdaSet, space.getEmptyConformation());
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
                /* TODO: handle right conformation lists */;
            }
            else root.child = newChild;
        }
        return root;
    }

    private static void createConformationTrie(TreeNode node, BWMStarNode root, SolutionSpace space)
    {
        LinkedHashSet<Integer> MSet = node.getCofEdge().getM();
        LinkedHashSet<Integer> lambda = node.getCofEdge().getLambda();
        Integer[] dummy = new Integer[]{};
        Position[] MArray = space.MSetFromArray(MSet).toArray(new Position[]{});
        root.conformationTrees = ProteinConformationTrie.createTrie(space, MArray, 0);
        root.generateConformations(MSet.toArray(dummy), lambda.toArray(dummy), 0, space, 
                space.getEmptyConformation(), space.getEmptyConformation());

    }

    private void generateConformations(Integer[] MSet, Integer[] lambda, int index, SolutionSpace space,
            Conformation current, Conformation tree)
    {
        System.out.println("Recurse! "+current+", "+tree+", "+index);
        if(index >= MSet.length + lambda.length)
        {
            System.out.println("Inserting "+tree+" to "+current);
            conformationTrees.insertConformation(current, tree);
            return;
        }
        if(index >= MSet.length)
        {
            Position p = space.positionFromPos(lambda[index - MSet.length]);
            for(Choice c : space.getChoices(p))
            {
                Conformation nextConf = tree.copy();
                nextConf.append(p, c);
                generateConformations(MSet, lambda, index + 1, space, current, nextConf);
            }

        }
        if(index < MSet.length)
        {
            Position p = space.positionFromPos(MSet[index]);            
            for(Choice c : space.getChoices(p))
            {
                Conformation nextConf = current.copy();
                nextConf.append(p, c);
                generateConformations(MSet, lambda, index + 1, space, nextConf, tree);
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

    public boolean moreConformations() {
        return moreConformations(emptyConformation);
    }

}
