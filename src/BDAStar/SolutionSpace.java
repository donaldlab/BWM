package BDAStar;

import java.util.ArrayList;
import java.util.Collection;

public interface SolutionSpace {

    public Collection<Choice> getChoices (Position p);
    public Conformation getEmptyConformation();

}
