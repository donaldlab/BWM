package BDAStar;

import java.util.Collection;

public class Position {
    int pos = 0;
    
    public Collection<Choice> getChoices()
    {
        return null;
    }
    
    public Position(int i, Choice c)
    {
        pos = i;
    }
    
    public Position(int i)
    {
        pos = i;
    }
    
    public boolean equals(Position p)
    {
    	return p.pos == pos;
    }

}
