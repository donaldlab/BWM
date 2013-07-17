package BDAStar;

import java.util.ArrayList;
import java.util.Collection;

public class SolutionSpace {

    public Collection<Choice> getChoices (Position p) {
        // TODO Auto-generated method stub
        ArrayList<Choice> choices = new ArrayList<Choice>();
        for(int i = 0; i< 2; i++)
            choices.add(new Choice(i));
        return choices;
    }

}
