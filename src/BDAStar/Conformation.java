package BDAStar;

import java.util.Collection;

public interface Conformation {
    public void append (Choice c);

    public void delete (Choice c);
    
    public void deleteLast();

    public double score ();

    public Collection<Position> getPositions ();

    //Generate a new copy Conformation joining the conformation with it.
    public Conformation join (Conformation conformation);

}
