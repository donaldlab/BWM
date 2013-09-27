package BDAStar;

import java.util.Collection;
import java.util.Map;
import java.util.PriorityQueue;

public class ProteinConformationTrie {
    
    private BDAStarNode subroot;
    private ConformationMap children;
    
    public ProteinConformationTrie(BWMSolutionSpace space, Position p)
    {
        children.initialize(space, p); 
    }
    
    public ProteinConformationTrie(int numAA, int[] numRotForAA, Position p)
    {
    	//if p is null, we're the root.

    }
    
    public BDAStarNode getAStarRoot (Conformation partial, int index)
    {
        Collection<Position> positions = partial.getPositions();
        Position[] positionArray = positions.toArray(new Position[]{});
        return getAStarRoot(partial, positionArray, index);
    }

    private BDAStarNode getAStarRoot (Conformation partial, Position[] positions,
            int index) {
        if(index < positions.length - 1)
        {
        	Choice currentChoice = partial.getChoiceAt(positions[index]);
            ProteinChoice choice = (ProteinChoice)currentChoice;
            return children.get(currentChoice).getAStarRoot(partial, positions, index+1);
        }
        return subroot;
    }
    
    public void insertConformation(Conformation partial)
    {
       Position[] positions = partial.getPositions().toArray(new Position[]{});
       insertConformation(partial, positions, 0);
    }

    private void insertConformation (Conformation partial,
            Position[] positions, int index) {
        if(index < positions.length - 1)
        {
        	Choice currentChoice = partial.getChoiceAt(positions[index]);
            ProteinChoice choice = (ProteinChoice)partial.getChoiceAt(positions[index]);
            children.get(currentChoice).insertConformation(partial, positions, index+1);            
        }
        
    }

    public void resort () {
    }
    
    public static ProteinConformationTrie createTrie(BWMSolutionSpace space, Position[] MSet, int index)
    {
        Position currentPosition = MSet[index];
        ProteinConformationTrie root = new ProteinConformationTrie(space, currentPosition);
        for(Choice c : space.getChoices(currentPosition))
        {
            ProteinChoice pc = (ProteinChoice) c;
            root.children.put(c,createTrie(space, MSet, index + 1));
        }
        return root;
    }

}
