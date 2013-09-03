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

}
