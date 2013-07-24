package BDAStar;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

public abstract class AbstractConformation implements Conformation{

    private HashMap<Position, Choice> positions;
    public AbstractConformation()
    {
        positions = new HashMap<Position, Choice>();
    }
    
    public void append (Choice c) {
    }

    @Override
    public void delete (Position p) {
        // TODO Auto-generated method stub
        positions.remove(p);
    }

    @Override
    public double score () {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public Collection<Position> getPositions () {
        // TODO Auto-generated method stub
        return positions.keySet();
    }

    @Override
    public abstract Conformation join (Conformation nextConformation);

    @Override
    public void deleteLast () {
    	positions.remove(positions.size()-1);
    }

	@Override
	public void append(Position p, Choice c) {
		// TODO Auto-generated method stub
		positions.put(p, c);
	}

	@Override
	public Choice getChoiceAt(Position p) {
		return positions.get(p);
	}

}
