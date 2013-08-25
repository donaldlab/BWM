package kstar;

public class RotInfo<T>{
	public int curPos;
	public int curAA;
	public int curRot;
	public T state;
	
	public RotInfo(int curPos, int curAA, int curRot, T s) {
		this.curPos = curPos;
		this.curAA = curAA;
		this.curRot = curRot;
		this.state = s;
	}
	
	public String printCoord(){
		return "("+curPos+","+curAA+","+curRot+")";
	}
	
}