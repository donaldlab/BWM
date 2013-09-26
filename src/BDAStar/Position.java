package BDAStar;

public class Position implements Comparable<Position>{
    public int pos = -1;
    
    public Position(int i)
    {
        pos = i;
    }
    
    public boolean equals(Position p)
    {
    	return p.pos == pos;
    }
    
    public int compareTo(Position p)
    {
        return pos - p.pos;
    }
    
    public int hashCode()
    {
    	return pos;
    }


}
