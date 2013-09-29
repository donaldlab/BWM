package BDAStar;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import kstar.PrunedRotamers;
import kstar.RotInfo;
import kstar.RotTypeMap;
import kstar.RotamerLibrary;
import kstar.RotamerSearch;

public class BWMSolutionSpace implements SolutionSpace {
    private EnergyFunction energyFunction;
    private Map<Position, Collection<ProteinChoice>> choices;
    private int[] designIndexToStrandIndex;
    private int[][] strandDesignIndices;
    private int[] designIndexToStrandResidueIndex;
    private int[] inverseResidueMap;
    private RotamerLibrary library;
    
    private int numRotForAA[];
    private int numAAAtPosition[];
    public Collection<ProteinChoice> getChoices (Position p){
        return choices.get(p);
    }

    public BWMSolutionSpace(PrunedRotamers<Boolean> lib, EnergyFunction e, int[] mutRes2Strand, 
    		int[][] strandMut, int[] mutRes2StrandMutIndex, int[] invResMap, RotamerLibrary rl)
    {
        super();
        energyFunction = e;
        designIndexToStrandIndex = mutRes2Strand;
        strandDesignIndices = strandMut;
        designIndexToStrandResidueIndex = mutRes2StrandMutIndex;
        choices = new HashMap<Position, Collection<ProteinChoice>>();
        inverseResidueMap = invResMap;
        library = rl;
        
        
        numAAAtPosition = new int[designIndexToStrandIndex.length];
        numRotForAA = new int[rl.getNumAAallowed()];
        /* We have to port over the rotamer library here, I think it's the RotamerSearch class. */
        /* TODO: 
         * 1. Convert library's contents into <Position, Collection<Choice>> Mapping.
         * 2. Optimize? Store?
         * 
         * NEW TODO:
         * 1. Get PrunedRotAtRes
         * 2. For everything in PruntedRotAtRes, produce a collection at each residue.
         */
        for(RotInfo<Boolean> r: lib)
        {
            if(!r.state)
                continue;
            Position p = positionFromPos(r.curPos);
            if(!choices.containsKey(p))
                choices.put(p, new ArrayList<ProteinChoice>());
            if(choices.get(p).size() < 4)
                choices.get(p).add(new ProteinChoice(r.curAA, r.curRot));
            numRotForAA[r.curAA]++;
            numAAAtPosition[p.pos]++;
        }
        
        /*
         * A position is defined by strand and residue number, or by molecule residue number.
         * In the code, there's a strandMutIndex which maps residues to the 0-based index used instead.
         * We want the strand number and the strand-relative residue number, I think. Though we'll store all of them for now.
         * 		int str = mutRes2Strand[i];
				int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
				int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
         */
    }

    @Override
    public ProteinConformation getEmptyConformation () {
        // TODO Auto-generated method stub
        return new ProteinConformation(energyFunction);
    }

    public Conformation createFromArray (int[] curState, RotTypeMap[][] rtm) {
        // TODO Auto-generated method stub
        ProteinConformation conf = getEmptyConformation();
        for(int i = 0; i < curState.length; i++)
        {
            RotTypeMap rotamerMap = rtm[i][curState[i]];
            int position = rotamerMap.pos;
            int aminoAcid = rotamerMap.aa;
            int rotamer = rotamerMap.rot;
            Position pos = positionFromPos(position);
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

    @Override
    public Position positionFromPos(Integer position) {
        int str = designIndexToStrandIndex[position];
        int strResNum = strandDesignIndices[str][designIndexToStrandResidueIndex[position]];
        return new ProteinPosition(str, strResNum, position);
    }
    
    public String printConformation(ProteinConformation p)
    {
    	return "";
    }

    public Set<? extends Position> MSetFromArray (LinkedHashSet<Integer> m) {
        /*
         * 1. For each integer, get the strand and sequence numbers
         * 2. create the corresponding ProteinPosition
         */
        HashSet<Position> MSet = new HashSet<Position>();
        for(int position: m)
        {
            MSet.add(positionFromPos(inverseResidueMap[position]));
        }
        return MSet;
    }

    public int getRotamersForAminoAcid (int AANum) {
        // TODO Auto-generated method stub
        return numRotForAA[AANum];
    }

    public int getAminoAcidsAtPosition (Position p) {
        // TODO Auto-generated method stub
        return numAAAtPosition[p.pos];
    }

	@Override
	public Collection<? extends Position> getPositions() {
		// TODO Auto-generated method stub
		return choices.keySet();
	}

	@Override
	public ConformationMap createConformationMap(Position p) {
		return new ProteinConformationMap(this, p);
	}



}
