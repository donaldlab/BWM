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
    	int out = 0;
        for(Position p : getPositions())
        {
            out-=getChoiceAt(p).choice;
        }
        return out;
    }


    //Generate a new copy Conformation joining the conformation with it.
    public TestConformation join (TestConformation nextConformation) {
        // TODO Auto-generated method stub
        return null;
    }
    
    public String toString()
    {
    	String out = "[";
        for(Position p : getPositions())
        {
            out+=getChoiceAt(p).choice+",";
        }
        out+="]";
        return out;
    }

}
