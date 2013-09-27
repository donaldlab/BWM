package BDAStar;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

public class ProteinConformationMap implements Map<ProteinChoice, ProteinConformationTrie> {
	ProteinConformationTrie[][] values;
	

	public void initialize(BWMSolutionSpace space, Position p) {
		// TODO Auto-generated method stub
		values = new ProteinConformationTrie[space.getAminoAcidsAtPosition(p).size()][];
        for(int aminoAcidIndex : space.getAminoAcidsAtPosition(p))
        {
            values[aminoAcidIndex] = new ProteinConformationTrie[space.getRotamersForAminoAcid(aminoAcidIndex)];
        }
	}
	
	public void initialize(int numAA, int[] numRotamers)
	{
    	values = new ProteinConformationTrie[numAA][];
    	for(int i = 0; i < numAA; i++)
    	{
    		values[numAA] = new ProteinConformationTrie[numRotamers[i]];
    	}
	}

	@Override
	public void clear() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public boolean containsKey(Object key) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean containsValue(Object value) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Set<java.util.Map.Entry<ProteinChoice, ProteinConformationTrie>> entrySet() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ProteinConformationTrie get(Object key) {
		ProteinChoice choice = (ProteinChoice) key;
		return values[choice.aminoAcid][choice.rotamer];
	}

	@Override
	public boolean isEmpty() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Set<ProteinChoice> keySet() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ProteinConformationTrie put(ProteinChoice key,
			ProteinConformationTrie value) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void putAll(
			Map<? extends ProteinChoice, ? extends ProteinConformationTrie> m) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public ProteinConformationTrie remove(Object key) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int size() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public Collection<ProteinConformationTrie> values() {
		// TODO Auto-generated method stub
		return null;
	}

}
