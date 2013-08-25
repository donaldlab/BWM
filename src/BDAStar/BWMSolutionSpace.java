package BDAStar;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import kstar.RotamerLibrary;

public class BWMSolutionSpace implements SolutionSpace {
    private Map<Position, Collection<Choice>> choices;
    public Collection<Choice> getChoices (Position p){
        if(choices.size() < 1)
        {
            Collection<Choice> out = new ArrayList<Choice>();
            out.add(new Choice(0));
            out.add(new Choice(1));
        }
        return choices.get(p);
    }
    
    public BWMSolutionSpace(RotamerLibrary library)
    {
        /* We have to port over the rotamer library here, I think it's the RotamerSearch class. */
        /* TODO: 
         * 1. Convert library's contents into <Position, Collection<Choice>> Mapping.
         * 2. Optimize? Store?
         */
        choices = new HashMap<Position, Collection<Choice>>();
    }

    @Override
    public Conformation getEmptyConformation () {
        // TODO Auto-generated method stub
        return null;
    }

}
