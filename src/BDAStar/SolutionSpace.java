package BDAStar;

import java.util.ArrayList;
import java.util.Collection;

public interface SolutionSpace {

    public Collection<? extends Choice> getChoices (Position p);
    public Conformation getEmptyConformation();
    public Collection<? extends Position> getPositions();
	public ConformationMap createConformationMap(Position p);

}
