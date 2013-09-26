package BDAStarTest;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

import BDAStar.Choice;
import BDAStar.Conformation;
import BDAStar.Position;
import BDAStar.ProteinChoice;
import BDAStar.SolutionSpace;
import kstar.PrunedRotamers;
import kstar.RotInfo;

public class TestSolutionSpace implements SolutionSpace {
    private Map<Position, Collection<ProteinChoice>> choices;
    private Collection<Choice> out;
    public Collection<Choice> getChoices (Position p){
        return out;
    }
    
    public TestSolutionSpace(int numChoices)
    {
        /* We have to port over the rotamer library here, I think it's the RotamerSearch class. */
        /* TODO: 
         * 1. Convert library's contents into <Position, Collection<Choice>> Mapping.
         * 2. Optimize? Store?
         * 
         * NEW TODO:
         * 1. Get PrunedRotAtRes
         * 2. For everything in PruntedRotAtRes, produce a collection at each residue.
         */
    	out = new ArrayList<Choice>();
    	for(int i = 0; i < numChoices; i ++)
    	{
    		out.add(new Choice(i));
    	}
    }

    @Override
    public Conformation getEmptyConformation () {
        // TODO Auto-generated method stub
        return new TestConformation();
    }

}
