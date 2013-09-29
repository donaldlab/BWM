package BDAStarTest;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

import BDAStar.BWMSolutionSpace;
import BDAStar.Choice;
import BDAStar.ConformationMap;
import BDAStar.Position;
import BDAStar.ProteinChoice;
import BDAStar.ProteinConformationTrie;

public class TestConformationMap extends ConformationMap {
	ProteinConformationTrie[] values;
	
	public TestConformationMap(int numChoices)
	{
		initialize(numChoices);
	}
	
	
	public void initialize(int numChoices)
	{
    	values = new ProteinConformationTrie[numChoices];
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
		Choice choice = (Choice) key;
		return values[choice.choice];
	}

	@Override
	public boolean isEmpty() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public ProteinConformationTrie put(Choice key,
			ProteinConformationTrie value) {
		values[key.choice] = value;
		return values[key.choice];
	}

	@Override
	public ProteinConformationTrie remove(Object key) {
		Choice choice = (Choice) key;
		ProteinConformationTrie out = values[choice.choice];
		return out;
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
