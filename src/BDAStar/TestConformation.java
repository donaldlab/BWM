package BDAStar;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
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
    		//	System.out.println("WARNING: redundant position joining.");
    		}
    		//System.out.println("Add "+testConformation.getChoiceAt(p).choice+" at "+p.pos);
    		append(p, testConformation.getChoiceAt(p));
    	}
    	for(Position p : nextConformation.getPositions())
    	{
    		if(testConformation.getChoiceAt(p) != null){
    		//	System.out.println("WARNING: redundant position joining.");
    		}
    		//System.out.println("Add "+nextConformation.getChoiceAt(p).choice+" at "+p.pos);
    		
    		append(p, nextConformation.getChoiceAt(p));
    	}
	}

	public void append (Choice c) {
    	getPositions().add(new Position(getPositions().size(),c));
        
    }

    
    public double score () {
    	double out = 0;
        for(Position p : getPositions())
        {
            if(p.pos == 5)
            {
            //System.out.println("Choice: "+getChoiceAt(p).choice+" Position: "+p.pos);
            //System.out.println("Adding "+1.0*Math.abs(getChoiceAt(p).choice-p.pos)/(Math.max(1, p.pos)));
            }
            out+=1.0*Math.abs(1.0*getChoiceAt(p).choice-p.pos)/(Math.max(1, p.pos));
        }
        return out;
    }
    
    public boolean equals(Object o)
    {
    	if(!o.getClass().equals(this.getClass()))
    		super.equals(o);
    	TestConformation t = (TestConformation) o;
    	if(t.getPositions().size() != getPositions().size())
    		return false;
    	boolean same = true;
    	
    	for(Position p: getPositions())
    	{
    		if(getChoiceAt(p) != t.getChoiceAt(p))
    			same = false;
    	}
    	return false;
    }


    //Generate a new copy Conformation joining the conformation with it.
    public Conformation join (Conformation nextConformation) {
        //System.out.println("Joining "+this+" and "+nextConformation);
        return new TestConformation(this, nextConformation);
    }
    
    public String toString()
    {
    	String out = "[";
    	Position[] c = getPositions().toArray(new Position[]{});
    	Arrays.sort(c);
        for(int i = 0; i < c.length; i++)
        {
            out+="("+c[i].pos+": "+getChoiceAt(c[i]).choice+"), ";
        }
        out+="]";
        return out;
    }

}
