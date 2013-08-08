package BDAStar;

import java.util.Collection;

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

}
