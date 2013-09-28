package BDAStar;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

public class ProteinConformationMap extends ConformationMap {
	ProteinConformationTrie[][] values;
	
	public ProteinConformationMap(BWMSolutionSpace space, Position p)
	{
		initialize(space, p);
	}
	

	public void initialize(BWMSolutionSpace space, Position p) {
		values = new ProteinConformationTrie[space.getAminoAcidsAtPosition(p)][];
        for(int i = 0; i < space.getAminoAcidsAtPosition(p); i++)
        {
            values[i] = new ProteinConformationTrie[space.getRotamersForAminoAcid(i)];
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
	public ProteinConformationTrie put(Choice key,
			ProteinConformationTrie value) {
		ProteinChoice choice = (ProteinChoice)key;
		return values[choice.aminoAcid][choice.rotamer];
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
