package BDAStar;

import java.util.Collection;
import java.util.Set;

public interface Conformation {
    public void append (Position p, Choice c);

    public void delete (Position p);
    
    public void deleteLast();

    public double score ();

    public Collection<Position> getPositions ();

    //Generate a new copy Conformation joining the conformation with it.
    public Conformation join (Conformation conformation);

	public Choice getChoiceAt(Position p);

	public Conformation copy();

    public void assignScore (double score);
    
    public Conformation extract(Set<? extends Position> target);

	public Conformation deletePositions(Collection<Position> toDelete);
}
