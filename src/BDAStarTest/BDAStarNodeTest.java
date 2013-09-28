package BDAStarTest;

import BDAStar.BDAStarNode;
import BDAStar.Choice;
import BDAStar.Conformation;
import BDAStar.Position;
import BDAStar.SolutionSpace;

public class BDAStarNodeTest {
    /**
     * Tests to write:
     * 1. Correct insertion
     * 2. Correct deletion
     * 3. Correct enumeration
     */

    public static int numPositions = 10;

    public static void main (String[] args)
    {
        System.out.println("Start!");
        TestSolutionSpace space = new TestSolutionSpace(10, 2);
        BDAStarNode root = new BDAStarNode(null, null, space.getEmptyConformation());
        insertConformations(root, new TestConformation(), space, 0);
        int rank = 0;
        while(root.moreConformations())
        {	
        	rank++;
            Conformation next = root.getNextConformation();
            System.out.println("Result "+rank+": "+next+" "+next.score());
        }
    }

    private static void insertConformations(BDAStarNode root, Conformation current, TestSolutionSpace space, int index)
    {
        if(index == numPositions)
        {
            //System.out.println("Inserting "+current);
            root.insertConformation(current);
            return;
        }
        Position p = new Position(index);
        for(Choice c : space.getChoices(p))
        {
            Conformation nextConf = current.copy();
            nextConf.append(p, c);
            insertConformations(root, nextConf, space, index+1);
        }

    }

}
