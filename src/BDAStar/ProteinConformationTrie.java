package BDAStar;

import java.util.Collection;
import java.util.Map;
import java.util.PriorityQueue;

public class ProteinConformationTrie {
    
    private BDAStarNode subroot;
    private ConformationMap children;
    private SolutionSpace space;
    
    public ProteinConformationTrie(SolutionSpace newSpace, ConformationMap map, Position p)
    {
    	space = newSpace;
    	children = map;
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
            return children.get(currentChoice).getAStarRoot(partial, positions, index+1);
        }
        return subroot;
    }
    
    public void insertConformation(Conformation partial, Conformation tree)
    {
       Position[] positions = partial.getPositions().toArray(new Position[]{});
       insertConformation(partial, tree, positions, 0);
    }

    private void insertConformation (Conformation partial, Conformation tree,
            Position[] positions, int index) {
        if(index < positions.length - 1)
        {
        	Position p = positions[index];
        	Choice currentChoice = partial.getChoiceAt(positions[index]);
        	if(children.get(currentChoice) == null)
        		children.put(currentChoice, new ProteinConformationTrie(space, 
        				space.createConformationMap(p), p));
            children.get(currentChoice).insertConformation(partial, tree, positions, index+1);            
        }
        if(subroot==null)
        	subroot = new BDAStarNode(null, null, space.getEmptyConformation());
        else
        	subroot.insertConformation(tree);
        
    }

    public void resort () {
    }
    
    public static ProteinConformationTrie createTrie(SolutionSpace space, Position[] MSet, int index)
    {
        Position currentPosition = MSet[index];
        ProteinConformationTrie root = new ProteinConformationTrie(space, space.createConformationMap(currentPosition), currentPosition);
        for(Choice c : space.getChoices(currentPosition))
        {
            ProteinChoice pc = (ProteinChoice) c;
            root.children.put(c,createTrie(space, MSet, index + 1));
        }
        return root;
    }

}