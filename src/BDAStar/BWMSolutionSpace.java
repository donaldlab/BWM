package BDAStar;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import kstar.RotamerLibrary;

public class BWMSolutionSpace implements SolutionSpace {
    private Map<Position, Collection<Choice>> choices;
    public Collection<Choice> getChoices (Position p){
        return choices.get(p);
    }
    
    public BWMSolutionSpace(RotamerLibrary library)
    {
        /* We have to port over the rotamer library here, I think it's the RotamerSearch class. */
        /* TODO: 
         * 1. Convert library's contents into <Position, Collection<Choice>> Mapping.
         * 2. Optimize? Store?
         */
    }

    @Override
    public Conformation getEmptyConformation () {
        // TODO Auto-generated method stub
        return null;
    }

}
