package BDAStar;

public class Position implements Comparable<Position>{
    public static final Position NULL_POSITION = new Position(-1);
	public int pos = -1;
    
    public Position(int i)
    {
        pos = i;
    }
    
    public boolean equals(Object o)
    {
    	if(o.getClass() != getClass())
    		return false;
    	Position p = (Position) o;
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
