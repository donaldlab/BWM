package BDAStarTest;

import BDAStar.BDAStarNode;
import BDAStar.Choice;
import BDAStar.Conformation;
import BDAStar.Position;
import BDAStar.ProteinConformationTrie;

public class ProteinConformationTrieTest {
	/**
	 * Tests to write:
	 * 1. Correct insertion
	 * 2. Correct creation
	 * 3. Memory/Runtime efficiency?
	 */
	
	static int numPositions = 10;
	
	static int numRotamers = 2;
	
	public static void main(String[] args)
	{
		int[] rotamers = new int[numPositions];
		for(int i = 0; i < rotamers.length; i++)
		{
			rotamers[i] = numRotamers;
		}
		ProteinConformationTrie trie = new ProteinConformationTrie(numPositions, rotamers, new Position(-1));
		TestSolutionSpace space = new TestSolutionSpace(2);
		
		insertConformations(trie, space.getEmptyConformation(), space, 0);
	}
	
    private static void insertConformations(ProteinConformationTrie root, Conformation current, TestSolutionSpace space, int index)
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
