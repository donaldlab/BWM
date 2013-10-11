package BDAStarTest;

import java.util.ArrayList;
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
	
	public static final int CHOICES = 2;
	
	private static void BuildBalancedTree(TreeNode root, int height, int currentHeight, 
			int currentIndex)
	{
		boolean isLeaf = height == currentHeight;
		TreeNode leftChild = new TreeNode(0, isLeaf, 0, 0);
		TreeNode rightChild = new TreeNode(0, isLeaf, 0, 0);
		LinkedHashSet<Integer> leftLambda = new LinkedHashSet<Integer>();
		LinkedHashSet<Integer> rightLambda = new LinkedHashSet<Integer>();
		TreeEdge leftEdge = new TreeEdge(0,1, leftLambda, null, null, null, 0, false);
		TreeEdge rightEdge = new TreeEdge(0,1, leftLambda, null, null, null, 0, false);
		for(int i = 0; i < CHOICES; i++)
		{
			leftLambda.add(i + currentIndex*CHOICES);
			rightLambda.add(i + (currentIndex+1)*CHOICES);
		}
		leftEdge.setPositions(leftLambda);
		rightEdge.setPositions(rightLambda);
		leftChild.setCofEdge(leftEdge);
		rightChild.setCofEdge(rightEdge);
		root.setLc(leftChild);
		root.setRc(rightChild);
		
		if(currentHeight < height)
		{
			BuildBalancedTree(leftChild, height, currentHeight + 1, currentIndex + 2);
			BuildBalancedTree(leftChild, height, currentHeight + 1, currentIndex + 
					(int)Math.pow(2, height - currentHeight));
		}
	}
    
    public static void main (String[] args)
    {
    	TreeNode testRoot = new TreeNode(0, false, 0, 0);
    	LinkedHashSet<Integer> testLambda = new LinkedHashSet<Integer>();
    	for(int i = 0; i < CHOICES; i ++)
        {
            testLambda.add(i);
        }
    	TreeEdge testRootEdge = new TreeEdge(0,1,testLambda,null,null,null,0,false);
    	testRootEdge.setPositions(testLambda);
    	testRoot.setCofEdge(testRootEdge);
    	
    	BuildBalancedTree(testRoot, 3, 0, 1);
    	testRoot.printTree("");

        SolutionSpace space = new TestSolutionSpace(10, 2);
        BWMAStarNode rootNode = new BWMAStarNode(space.getEmptyConformation());
        BWMAStarNode.CreateTree2(testRoot, rootNode, space.getEmptyConformation(), rootNode.getChildren(), space, 0, testRoot.getCofEdge().getPositionList());
        
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
            //rootNode.printTree("");
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
