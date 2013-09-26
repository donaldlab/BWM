package BDAStar;

import java.util.Collection;
import java.util.PriorityQueue;

public class ProteinConformationTrie {
    
    private BDAStarNode subroot;
    private ProteinConformationTrie[][] children;
    
    public ProteinConformationTrie(BWMSolutionSpace space, Position p)
    {
        children = new ProteinConformationTrie[space.getAminoAcidsAtPosition(p).size()][];
        for(int aminoAcidIndex : space.getAminoAcidsAtPosition(p))
        {
            children[aminoAcidIndex] = new ProteinConformationTrie[space.getRotamersForAminoAcid(aminoAcidIndex)];
        }
    }
    
    public BDAStarNode getAStarRoot (Conformation partial, int index)
    {
        Collection<Position> positions = partial.getPositions();
        Position[] positionArray = positions.toArray(new Position[]{});
        return getAStarRoot(partial, positionArray, index);
    }

    private BDAStarNode getAStarRoot (Conformation partial, Position[] positionArray,
            int index) {
        if(index < positionArray.length - 1)
        {
            ProteinChoice choice = (ProteinChoice)partial.getChoiceAt(positionArray[index]);
            return children[choice.aminoAcid][choice.rotamer].getAStarRoot(partial, positionArray, index+1);
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
            ProteinChoice choice = (ProteinChoice)partial.getChoiceAt(positions[index]);
            children[choice.aminoAcid][choice.rotamer].insertConformation(partial, positions, index+1);            
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
            root.children[pc.aminoAcid][pc.rotamer] = createTrie(space, MSet, index + 1);
        }
        return root;
    }

}
