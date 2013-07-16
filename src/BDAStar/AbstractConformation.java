package BDAStar;

import java.util.ArrayList;
import java.util.Collection;

public class AbstractConformation implements Conformation{

    ArrayList<Position> positions;
    public AbstractConformation()
    {
        positions = new ArrayList<Position>();
    }
    
    @Override
    public void append (Choice c) {
        positions = new ArrayList<Position>();
    }

    @Override
    public void delete (Choice c) {
        // TODO Auto-generated method stub
        
    }

    @Override
    public double score () {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public Collection<Position> getPositions () {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public Conformation join (Conformation nextConformation) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public void deleteLast () {
        // TODO Auto-generated method stub
        
    }

}
