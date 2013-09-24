package BDAStar;

public class ProteinChoice extends Choice {

    public int aminoAcid;
    public int rotamer;
    public String[] aminoAcidNames = {
    		"GLYCINE",
    		"ARGININE"
    };

    public ProteinChoice (int AA, int r) {
    	super(AA*100 + r);
        // TODO Auto-generated constructor stub
        aminoAcid = AA;
        rotamer = r;
    }
    
    public String toString()
    {
    	return aminoAcid+"-"+rotamer;
    }

}
