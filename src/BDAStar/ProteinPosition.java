package BDAStar;

public class ProteinPosition extends Position {

	private int strand;
	private int residue;
	public int designIndex;
    public ProteinPosition (int designInd, int strandNumber, int position) {
        super(position);
        designIndex = designInd;
        strand = strandNumber;
        residue = position;
    }
    
    public int hashCode()
    {
    	return strand*10000 + residue;
    }
    
    public boolean equals(Object o)
    {
    	if(!o.getClass().equals(getClass()))
    		return false;
    	ProteinPosition otherPosition = (ProteinPosition) o;
    	return otherPosition.strand == strand && otherPosition.residue == residue;
    }

}
