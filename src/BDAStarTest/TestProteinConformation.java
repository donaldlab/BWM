package BDAStarTest;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.Set;

import BDAStar.AbstractConformation;
import BDAStar.Choice;
import BDAStar.Conformation;
import BDAStar.EnergyFunction;
import BDAStar.Position;
import BDAStar.ProteinConformation;
import kstar.InteractionGraph;
import kstar.Molecule;
import kstar.PairwiseEnergyMatrix;

public class TestProteinConformation extends ProteinConformation {

    private static EnergyFunction energyFunction;
    public static final double MAX_ENERGY = Double.MAX_VALUE;
    private double score = MAX_ENERGY;

    public TestProteinConformation(EnergyFunction e)
    {
        super(null);
    }

    public TestProteinConformation(Conformation t, EnergyFunction e)
    {
        super(null);
        for(Position p : t.getPositions())
        {
            if(getChoiceAt(p) != null && t.getChoiceAt(p) != getChoiceAt(p))
            {
                System.out.println("OVERWRITING CHOICE at "+p.pos +": "+getChoiceAt(p)+" with "+t.getChoiceAt(p));
            }
            append(p, t.getChoiceAt(p));
        }
    }

    public TestProteinConformation(Conformation testConformation,
            Conformation nextConformation) {
    	super(null);
        
        for(Position p : testConformation.getPositions())
        {
            if(testConformation.getChoiceAt(p) != null && nextConformation.getChoiceAt(p) != null && nextConformation.getChoiceAt(p) != testConformation.getChoiceAt(p))
            {

            	System.out.println("OVERWRITING CHOICE at "+p.pos +": "+testConformation.getChoiceAt(p)+" with "+nextConformation.getChoiceAt(p));
            }
        }
        
        for(Position p : testConformation.getPositions())
        {
            if(nextConformation.getChoiceAt(p) != null){
          //      System.out.println("WARNING: redundant position joining.");
            }
            //System.out.println("Add "+testConformation.getChoiceAt(p).choice+" at "+p.pos);
            append(p, testConformation.getChoiceAt(p));
        }
        for(Position p : nextConformation.getPositions())
        {
            if(testConformation.getChoiceAt(p) != null){
              //  System.out.println("WARNING: redundant position joining.");
            }
            //System.out.println("Add "+nextConformation.getChoiceAt(p).choice+" at "+p.pos);

            append(p, nextConformation.getChoiceAt(p));
        }
    }

    public void append (Choice c) 
    {
        getPositions().add(new Position(getPositions().size()));

    }


    public double score () {
        return score;
    }
    
    public void assignScore(double s)
    {
        score = s;
    }

    public boolean equals(Object o)
    {
        if(!o.getClass().equals(this.getClass()))
            super.equals(o);
        TestProteinConformation t = (TestProteinConformation) o;
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
        return new TestProteinConformation(this, nextConformation);
    }





    public String toString()
    {
        String out = "[";
        Position[] c = getPositions().toArray(new Position[]{});
        Arrays.sort(c);
        for(int i = 0; i < c.length; i++)
        {
            out+="("+c[i].pos+": "+getChoiceAt(c[i]).toString()+"), ";
        }
        out+="]";
        return out;
    }

    @Override
    public Conformation copy() {
        // TODO Auto-generated method stub
        return new TestProteinConformation(this, energyFunction);
    }

}
