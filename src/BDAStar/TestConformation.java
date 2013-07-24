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

    public TestConformation(Conformation testConformation,
			Conformation nextConformation) {
    	for(Position p : testConformation.getPositions())
    	{
    		if(nextConformation.getChoiceAt(p) != null){
    			System.out.println("WARNING: redundant position joining.");
    		}
    		System.out.println("Add "+testConformation.getChoiceAt(p).choice+" at "+p.pos);
    		append(p, testConformation.getChoiceAt(p));
    	}
    	for(Position p : nextConformation.getPositions())
    	{
    		if(testConformation.getChoiceAt(p) != null){
    			System.out.println("WARNING: redundant position joining.");
    		}
    		System.out.println("Add "+nextConformation.getChoiceAt(p).choice+" at "+p.pos);
    		
    		append(p, nextConformation.getChoiceAt(p));
    	}
	}

	public void append (Choice c) {
    	getPositions().add(new Position(getPositions().size(),c));
        
    }

    
    public double score () {
    	int out = 0;
        for(Position p : getPositions())
        {
            out+=getChoiceAt(p).choice;
        }
        return out;
    }


    //Generate a new copy Conformation joining the conformation with it.
    public Conformation join (Conformation nextConformation) {
        // TODO Auto-generated method stub
    	System.out.println("Joining "+this+" and "+nextConformation);
        return new TestConformation(this, nextConformation);
    }
    
    public String toString()
    {
    	String out = "[";
        for(Position p : getPositions())
        {
            out+=p.pos+": "+getChoiceAt(p).choice+", ";
        }
        out+="]";
        return out;
    }

}
