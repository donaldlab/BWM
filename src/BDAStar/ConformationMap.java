package BDAStar;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

public abstract class ConformationMap implements Map<Choice, ProteinConformationTrie> {

	@Override
	public void clear() {
		// TODO Auto-generated method stub

	}
	
	public void initialize(int numAA, int[] numRotamers)
	{
		
	}
	
	public void initialize(SolutionSpace space, Position p)
	{
		
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
	public Set<java.util.Map.Entry<Choice, ProteinConformationTrie>> entrySet() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ProteinConformationTrie get(Object key) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isEmpty() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Set<Choice> keySet() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public void putAll(
			Map<? extends Choice, ? extends ProteinConformationTrie> m) {
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
