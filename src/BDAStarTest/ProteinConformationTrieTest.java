package BDAStarTest;

import java.util.HashSet;
import java.util.Set;

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
	
	static int numMPositions = 5;
	
	static int numRotamers = 2;
	
	public static void main(String[] args)
	{
		int[] rotamers = new int[numPositions];
		for(int i = 0; i < rotamers.length; i++)
		{
			rotamers[i] = numRotamers;
		}
		TestSolutionSpace space = new TestSolutionSpace(10, 2);
		
		ProteinConformationTrie trie = new ProteinConformationTrie(space, 
				space.createConformationMap(new Position(-1)), new Position(-1));
		
		insertConformations(trie, space.getEmptyConformation(), space.getEmptyConformation(), space, 0);
		HashSet<Conformation> MConfs = new HashSet<Conformation>();
		generateMSets(MConfs, space, space.getEmptyConformation(), 0);
		int MConfCount = 0;
		for(Conformation mConf : MConfs)
		{
			MConfCount++;
			if(MConfCount > 2) return;
			System.out.println("Results for "+mConf);
			BDAStarNode lambdaConfNode = trie.getAStarRoot(mConf, 0);
			int rank = 0;
			while(lambdaConfNode.moreConformations() && rank < 1024)
			{
				rank ++;
				System.out.println("Result "+rank+": "+lambdaConfNode.getNextConformation());
			}
		}
	}
	
	private static void generateMSets(Set<Conformation> set, TestSolutionSpace space, 
			Conformation current, int index)
	{
		if(index < numMPositions)
		{
			Position p = new Position(index);
			for(Choice c: space.getChoices(p))
			{
				Conformation nextConf = current.copy();
				nextConf.append(p, c);
				generateMSets(set, space, nextConf, index + 1);
			}
			return;
		}
		set.add(current);
		
	}
	
    private static void insertConformations(ProteinConformationTrie root, Conformation current, 
    		Conformation lambda, TestSolutionSpace space, int index)
    {
        if(index == numPositions)
        {
            //System.out.println("Inserting "+current);
            root.insertConformation(current, lambda);
            return;
        }
        Position p = new Position(index);
        Conformation target = current;
        if(index > numMPositions)
        	target = lambda;
        for(Choice c : space.getChoices(p))
        {
            Conformation nextConf = target.copy();
            nextConf.append(p, c);
            if(index < numMPositions)
            	insertConformations(root, nextConf, current, space, index+1);
            else
            	insertConformations(root, current, nextConf, space, index + 1);
        }

    }

}
