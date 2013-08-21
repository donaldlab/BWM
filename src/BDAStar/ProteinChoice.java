package BDAStar;

public class ProteinChoice extends Choice {

    public int aminoAcid;
    public int rotamer;
    public ProteinChoice (int i, int AA, int r) {
        super(i);
        // TODO Auto-generated constructor stub
        aminoAcid = AA;
        rotamer = r;
    }

}
