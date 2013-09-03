package BDAStar;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import kstar.PrunedRotamers;
import kstar.RotInfo;
import kstar.RotTypeMap;

public class BWMSolutionSpace implements SolutionSpace {
	private EnergyFunction energyFunction;
    private Map<Position, Collection<ProteinChoice>> choices;
    public Collection<ProteinChoice> getChoices (Position p){
        return choices.get(p);
    }
    
    public BWMSolutionSpace(PrunedRotamers<Boolean> library, EnergyFunction e)
    {
    	energyFunction = e;
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
    			choices.put(p, new ArrayList<ProteinChoice>());
    		choices.get(p).add(new ProteinChoice(r.curAA, r.curRot));
    	}
    }

    @Override
    public Conformation getEmptyConformation () {
        // TODO Auto-generated method stub
        return new ProteinConformation();
    }
    
    public static Conformation createFromArray (int[] curState, RotTypeMap[][] rtm) {
        // TODO Auto-generated method stub
        ProteinConformation conf = new ProteinConformation();
    	for(int i = 0; i < curState.length; i++)
    	{
    		RotTypeMap rotamerMap = rtm[i][curState[i]];
			int position = rotamerMap.pos;
			int aminoAcid = rotamerMap.aa;
			int rotamer = rotamerMap.rot;
		ProteinPosition pos = new ProteinPosition(position);
		ProteinChoice choice = new ProteinChoice(aminoAcid, rotamer);
		conf.append(pos, choice);
    	}
    	/*
    	 * 1. Remap the indices to residue Positions
    	 * 2. Remap each array element (an int) into a Rotamer from rtm array
    	 * 3. Create a new Conformation, adding the Choice at each position.
    	 */
        return conf;
    }

    public static Set<Position> MSetFromArray (LinkedHashSet<Integer> m) {
        // TODO Auto-generated method stub
        return null;
    }

}
