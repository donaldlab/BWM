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
    
    public BDAStarNode getAStarRoot(Conformation partial)
    {
    	return getAStarRoot(partial, 0);
    }
    
    private BDAStarNode getAStarRoot (Conformation partial, int index)
    {
        Collection<Position> positions = partial.getPositions();
        Position[] positionArray = positions.toArray(new Position[]{});
        return getAStarRoot(partial, positionArray, index);
    }

    private BDAStarNode getAStarRoot (Conformation partial, Position[] positions,
            int index) {
        if(index < positions.length)
        {
        	Choice currentChoice = partial.getChoiceAt(positions[index]);
            return children.get(currentChoice).getAStarRoot(partial, positions, index+1);
        }
        return subroot;
    }
    
    public void insertConformation(Conformation MConf, Conformation lambdaConf, Conformation suffix)
    {
       Position[] positions = MConf.getPositions().toArray(new Position[]{});
       insertConformation(MConf, lambdaConf, suffix, positions, 0);
    }

    private void insertConformation (Conformation MConf, Conformation lambdaConf, Conformation suffix,
            Position[] positions, int index) {
        if(index < positions.length)
        {
        	Position p = positions[index];
        	Choice currentChoice = MConf.getChoiceAt(positions[index]);
        	if(children.get(currentChoice) == null)
        		children.put(currentChoice, new ProteinConformationTrie(space, 
        				space.createConformationMap(p), p));
            children.get(currentChoice).insertConformation(MConf, lambdaConf, suffix, positions, index+1);            
        }
        if(subroot==null)
        	subroot = new BDAStarNode(null, null, space.getEmptyConformation());
        	subroot.insertConformation(lambdaConf, MConf, suffix);
        
    }

    public void resort () {
    }
    
    public static ProteinConformationTrie createTrie(SolutionSpace space, Position[] MSet, int index)
    {
        Position currentPosition = Position.NULL_POSITION;
        if(index < MSet.length)
        {
        	currentPosition = MSet[index];
        }
        ProteinConformationTrie root = new ProteinConformationTrie(space, space.createConformationMap(currentPosition), currentPosition);
        if(index < MSet.length)
        for(Choice c : space.getChoices(currentPosition))
        {
            root.children.put(c,createTrie(space, MSet, index + 1));
        }
        return root;
    }

}
