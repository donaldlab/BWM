package BDAStar;

public class ProteinChoice extends Choice {

    public int aminoAcid;
    public int rotamer;
    public ProteinChoice (int AA, int r) {
    	super(AA);
        // TODO Auto-generated constructor stub
        aminoAcid = AA;
        rotamer = r;
    }

}
