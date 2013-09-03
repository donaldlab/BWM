package kstar;

//Maps state indices to residue/aa/rot tuples
public class RotTypeMap {
	public int pos = -1;
	public int aa = -1;
	public int rot = -1;
	public RotTypeMap(int p, int a, int r){
		pos = p;
		aa = a;
		rot = r;
	}
}