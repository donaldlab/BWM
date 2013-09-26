package BDAStar;

public class ProteinChoice extends Choice {

    public int aminoAcid;
    public int rotamer;
    public String[] aminoAcidNames;

    public ProteinChoice (int AA, int r) {
    	super(AA*100 + r);
        // TODO Auto-generated constructor stub
        aminoAcid = AA;
        rotamer = r;
    }
    
    public boolean equals(Object o)
    {
        if(o.getClass() != this.getClass())
            return false;
        ProteinChoice c = (ProteinChoice) o;
        return (c.aminoAcid == aminoAcid && c.rotamer == rotamer);
    }
    
    public int hashCode()
    {
        return choice;
    }

}
