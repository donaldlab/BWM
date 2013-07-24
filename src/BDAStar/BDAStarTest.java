package BDAStar;

import java.util.LinkedHashSet;
import kstar.TreeEdge;
import kstar.TreeNode;

public class BDAStarTest {
    
    public static void main (String[] args)
    {
        System.out.println("HOORAY");
        TreeNode root = new TreeNode(0, false, 0, 0);
        TreeNode lchild = new TreeNode(0, true, 0, 0);
        TreeNode rchild = new TreeNode(0, true, 0, 0);
    	LinkedHashSet<Integer> lambda = new LinkedHashSet<Integer>();
    	LinkedHashSet<Integer> lambda2 = new LinkedHashSet<Integer>();
    	LinkedHashSet<Integer> lambda3 = new LinkedHashSet<Integer>();
        for(int i = 0; i < 2; i ++)
        {
            lambda.add(i);
            lambda2.add(i+2);
            lambda3.add(i+4);
        }
        TreeEdge edge = new TreeEdge(0, 1, lambda, null, null, null, 0, false);
        edge.setPositions(lambda);
        TreeEdge edge2 = new TreeEdge(0, 1, lambda2, null, null, null, 0, false);
        edge2.setPositions(lambda2);
        TreeEdge edge3 = new TreeEdge(0, 1, lambda3, null, null, null, 0, false);
        edge3.setPositions(lambda3);
        
        
        root.setCofEdge(edge);
        lchild.setCofEdge(edge2);
        rchild.setCofEdge(edge3);
        root.setLc(lchild);
        root.setRc(rchild);
        SolutionSpace space = new SolutionSpace();
        BDAStarNode rootNode = BDAStarNode.CreateTree(root, space);
        System.out.println("Done!");
        rootNode.printTree("");
        System.out.println("RESULT: "+rootNode.getNextConformation());
    }

}
