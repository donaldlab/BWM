package BDAStarTest;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;

import BDAStar.BWMAStarNode;
import BDAStar.Choice;
import BDAStar.Conformation;
import BDAStar.Position;
import BDAStar.SolutionSpace;
import kstar.PrunedRotamers;
import kstar.RotamerLibrary;
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
    	int choices = 2;
        for(int i = 0; i < choices; i ++)
        {
            lambda.add(i);
            lambda2.add(i+choices);
            lambda3.add(i+2*choices);
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
        SolutionSpace space = new TestSolutionSpace(2, 10);
        //BDAStarNode rootNode = BDAStarNode.CreateTree(root, new TestConformation(), space);
        //root.setEnumerationObjects(null, space);
        //BWMAStarNode rootNode = BWMAStarNode.CreateTree(root, space.getEmptyConformation(), space, 0);
        BWMAStarNode rootNode = new BWMAStarNode(space.getEmptyConformation());
        BWMAStarNode.CreateTree2(root, rootNode, space.getEmptyConformation(), rootNode.getChildren(), space, 0, root.getCofEdge().getPositionList());
        
        rootNode.resort();
        System.out.println("Done!");
        rootNode.printTree("");
        int rank = 0;
        HashSet<String> solutions = new HashSet<String>();
        double lastScore = -1000;
        System.out.println("Total solutions: "+rootNode.totalPossibleCombinations());
        while(rootNode.moreConformations() && rank < 500)
        {
            rank++;
            Conformation c = rootNode.getNextConformation();
            if(lastScore > c.score())
            	System.out.println("WRONG ORDER");
            System.out.println("RESULT "+rank+": "+c+", SCORE: "+c.score());
            lastScore = c.score();
            
            if(rank == 57)
            {
                System.out.println("");
            }

            if(solutions.contains(c.toString()))
            	System.out.println("DUPLICATE!!!");
            solutions.add(c.toString());
            rootNode.printTree("");
            if(rank + rootNode.remainingConformations() != 64)
            {
                System.out.println("Dur?");
            }
            System.out.println("Remaining conformations: "+rootNode.remainingConformations());

        }
        
        Conformation test = new TestConformation();
        for(int i = 0; i < 6; i++)
        {
            int whee = i;
            if(i == 5)
                whee = 4;
            test.append(new Position(i), new Choice(whee));
        }
        
        HashSet<String> fullSpace = new HashSet<String>(); 
        buildAllSolutions(fullSpace, new TestConformation(), 0, 5, 2);
        System.out.println("Total possible solutions: "+fullSpace.size());
        fullSpace.removeAll(solutions);
        System.out.println("Missing solutions: "+fullSpace.size());
        for(String s: fullSpace)
        {
            System.out.println(s+": "+scoreString(s));
        }
        

    }
    
    public static double scoreString(String conf)
    {
        conf = conf.substring(1, conf.length()-1);
        String[] pairs = conf.split(", ");
        double out = 0;
        for(String pair : pairs)
        {
            pair = pair.substring(1, pair.length() - 1);
            Integer position = Integer.valueOf(pair.substring(0,1));

            Integer choice = Integer.valueOf(pair.substring(pair.length()-1),pair.length());
            out+=1.0*Math.abs(1.0*choice-position)/(Math.max(1, position));
        }
        return out;
    }
    
    public static void buildAllSolutions(HashSet<String> out, Conformation current, int index, int positions, int choices)
    {
        //System.out.println(current);
        for(int choice = 0; choice < choices; choice ++)
        {    
            Conformation newConf = current.copy();
            newConf.append(new Position(index), new Choice(choice));
           if(index == positions)
           {
               out.add(newConf.toString());
           }
        
           else buildAllSolutions(out, newConf, index+1, positions, choices);
        }
    }

}