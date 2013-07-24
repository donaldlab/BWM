package BDAStar;

import java.util.LinkedHashSet;
import kstar.TreeEdge;
import kstar.TreeNode;

public class BDAStarTest {
    
    public static void main (String[] args)
    {
        System.out.println("HOORAY");
        TreeNode root = new TreeNode(0, true, 0, 0);
        TreeNode child = new TreeNode(0, false, 0, 0);
        LinkedHashSet<Integer> lambda = new LinkedHashSet<Integer>();
        TreeEdge edge = new TreeEdge(0, 1, lambda, null, null, null, 0, true);
        root.setCofEdge(edge);
        SolutionSpace space = new SolutionSpace();
        BDAStarNode rootNode = BDAStarNode.CreateTree(root, space);
        System.out.println("Done!");
        rootNode.getNextConformation();
    }

}
