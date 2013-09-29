package BDAStar;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.Set;

public interface SolutionSpace {

    public Collection<? extends Choice> getChoices (Position p);
    public Conformation getEmptyConformation();
    public Collection<? extends Position> getPositions();
	public ConformationMap createConformationMap(Position p);
	public Set<? extends Position> MSetFromArray(LinkedHashSet<Integer> MSet);
	public Position positionFromPos(Integer integer);

}
