package BDAStar;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public abstract class AbstractConformation implements Conformation{

    private HashMap<Position, Choice> positions;
    private double score;
    public AbstractConformation()
    {
        positions = new HashMap<Position, Choice>();
    }

    @Override
    public void delete (Position p) {
        positions.remove(p);
    }


    @Override
    public Collection<Position> getPositions () {
        return positions.keySet();
    }

    @Override
    public abstract Conformation join (Conformation nextConformation);

    public Conformation deletePositions(Collection<Position> toDelete)
    {
        Conformation out = copy();
        for(Position p: toDelete)
            out.delete(p);
        return out;
    }

    private Set<Position> copyPositionSet()
    {
        HashSet<Position> out = new HashSet<Position>();
        for(Position p : getPositions())
        {
            out.add(p);
        }
        return out;
    }

    @Override
    public Conformation extract(Set<? extends Position> target)
    {

        Collection<Position> toDelete = copyPositionSet();
        toDelete.removeAll(target);
        Conformation out = deletePositions(toDelete);
        return out;
    }

    public void assignScore(double s)
    {
        score = s;
    }

    @Override
    public void deleteLast () {
        positions.remove(positions.size()-1);
    }

    @Override
    public void append(Position p, Choice c) {
        // TODO Auto-generated method stub
        if(positions.containsKey(p) && !positions.get(p).equals(c))
        {
            System.out.println("FATAL OVERWRITE on "+p.pos+" "+positions.get(p).choice+" with "+c.choice);
        }
        positions.put(p, c);
    }

    @Override
    public Choice getChoiceAt(Position p) {
        return positions.get(p);
    }

}
