package BDAStar;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.Set;

public class TestConformation extends AbstractConformation {
 
    public TestConformation()
    {
        super();
    }
    
    public TestConformation(Conformation t)
    {
        for(Position p : t.getPositions())
        {
            positions.add(p);
        }
    }

    public void append (Choice c) {
        positions.add(new Position(c));
        // TODO Auto-generated method stub
        
    }

    public void delete (Choice c) {
        // TODO Auto-generated method stub
        int index = positions.lastIndexOf(c);
        if(index >= 0)
            positions.remove(index);
        
    }
    
    public void deleteLast()
    {
        positions.remove(positions.size()-1);
    }

    public double score () {
        // This will require a scoring function...
        return 0;
    }

    public Collection<Position> getPositions () {
        // TODO Auto-generated method stub
        /*
        if(positions.size() < 1)
        {
            System.out.println("Generating more positions...");
            for(int i = 0; i < 10; i ++)
                positions.add(new Position(i));
        }*/
        return positions;
    }

    //Generate a new copy Conformation joining the conformation with it.
    public TestConformation join (TestConformation nextConformation) {
        // TODO Auto-generated method stub
        return null;
    }

}
