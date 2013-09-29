package BDAStarTest;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import BDAStar.Choice;
import BDAStar.Conformation;
import BDAStar.ConformationMap;
import BDAStar.Position;
import BDAStar.ProteinChoice;
import BDAStar.ProteinPosition;
import BDAStar.SolutionSpace;
import kstar.PrunedRotamers;
import kstar.RotInfo;

public class TestSolutionSpace implements SolutionSpace {
    private Map<Position, Collection<ProteinChoice>> choices;
    private Collection<Choice> out;
    private Collection<Position> positions;
    public Collection<Choice> getChoices (Position p){
        return out;
    }
    
    public TestSolutionSpace(int numPositions, int numChoices)
    {

    	out = new ArrayList<Choice>();
    	for(int i = 0; i < numChoices; i ++)
    	{
    		out.add(new Choice(i));
    	}
    	positions = new ArrayList<Position>();
		for(int i = 0; i < numPositions; i++)
		{
			positions.add(new Position(i));
		}
    }

    @Override
    public Conformation getEmptyConformation () {
        // TODO Auto-generated method stub
        return new TestConformation();
    }

	@Override
	public Collection<? extends Position> getPositions() {
		return positions;
	}

	@Override
	public ConformationMap createConformationMap(Position p) {
		// TODO Auto-generated method stub
		return new TestConformationMap(out.size());
	}
	
    public Set<? extends Position> MSetFromArray (LinkedHashSet<Integer> m) {
        /*
         * 1. For each integer, get the strand and sequence numbers
         * 2. create the corresponding ProteinPosition
         */
        HashSet<Position> MSet = new HashSet<Position>();
        for(int position: m)
        {
            MSet.add(new Position(position));
        }
        return MSet;
    }

	@Override
	public Position positionFromPos(Integer integer) {
		// TODO Auto-generated method stub
		return new Position(integer);
	}

}
