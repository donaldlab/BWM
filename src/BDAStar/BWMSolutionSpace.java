package BDAStar;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

import kstar.PrunedRotamers;
import kstar.RotInfo;

public class BWMSolutionSpace implements SolutionSpace {
    private Map<Position, Collection<Choice>> choices;
    public Collection<Choice> getChoices (Position p){
        return choices.get(p);
    }
    
    public BWMSolutionSpace(PrunedRotamers<Boolean> library)
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
    	for(RotInfo<Boolean> r: library)
    	{
    		if(!r.state)
    			continue;
    		Position p = new Position(r.curPos);
    		if(!choices.containsKey(p))
    			choices.put(p, new ArrayList<Choice>());
    		choices.get(p).add(new ProteinChoice(r.curAA, r.curRot));
    	}
    }

    @Override
    public Conformation getEmptyConformation () {
        // TODO Auto-generated method stub
        return new ProteinConformation();
    }

}
