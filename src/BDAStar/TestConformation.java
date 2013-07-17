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
    	super();
        for(Position p : t.getPositions())
        {
            append(p, t.getChoiceAt(p));
        }
    }

    public void append (Choice c) {
    	getPositions().add(new Position(getPositions().size(),c));
        
    }

    
    public double score () {
        // This will require a scoring function...
        return 0;
    }


    //Generate a new copy Conformation joining the conformation with it.
    public TestConformation join (TestConformation nextConformation) {
        // TODO Auto-generated method stub
        return null;
    }

}
