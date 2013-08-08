package BDAStar;

import java.util.HashMap;
import java.util.HashSet;
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
        BWMAStarNode rootNode = BWMAStarNode.CreateTree(root, new TestConformation(), space);
        rootNode.resort();
        System.out.println("Done!");
        rootNode.printTree("");
        Conformation c = rootNode.getNextConformation();
        int rank = 0;
        HashSet<String> solutions = new HashSet<String>();
        while(c!= null && rank < 100)
        {
            rank++;
            if(rank == 4)
            {
            	System.out.println("Catch!");
            }
            System.out.println("RESULT "+rank+": "+c+", SCORE: "+c.score());

            if(solutions.contains(c.toString()))
            	System.out.println("DUPLICATE!!!");
            solutions.add(c.toString());
            c = rootNode.getNextConformation();
        }
        
        Conformation test = new TestConformation();
        for(int i = 0; i < 6; i++)
        {
            int whee = i;
            if(i == 5)
                whee = 4;
            test.append(new Position(i), new Choice(whee));
        }
        System.out.println(test+", "+test.score());

    }

}
