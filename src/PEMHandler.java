/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 1.0
	Copyright (C) 2001-2009 Bruce Donald Lab, Duke University
	
	OSPREY is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as 
	published by the Free Software Foundation, either version 3 of 
	the License, or (at your option) any later version.
	
	OSPREY is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Lesser General Public License for more details.
	
	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, see:
	      <http://www.gnu.org/licenses/>.
		
	There are additional restrictions imposed on the use and distribution
	of this open-source code, including: (A) this header must be included
	in any modification or extension of the code; (B) you are required to
	cite our papers in any publications that use this code. The citation
	for the various different modules of our software, together with a
	complete list of requirements and restrictions are found in the
	document license.pdf enclosed with this distribution.
	
	Contact Info:
			Bruce Donald
			Duke University
			Department of Computer Science
			Levine Science Research Center (LSRC)
			Durham
			NC 27708-0129 
			USA
			e-mail:   www.cs.duke.edu/brd/
	
	<signature of Bruce Donald>, 12 Apr, 2009
	Bruce Donald, Professor of Computer Science
*/

///////////////////////////////////////////////////////////////////////////////////////////////
//	PEMHandler.java
//
//	Version:           1.0
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//
///////////////////////////////////////////////////////////////////////////////////////////////

/**
* Written by Ivelin Georgiev (2004-2009)
* 
*/

/**
 * Manages operations on the pairwise energy matrices.
 * 
*/
public class PEMHandler {

	PEMHandler(){};
	
	/* 
	 * Initialize the pairwise energy matrices: each matrix has 6 dimensions;
	 * The first three dimensions correspond to the residue position, amino acid type, and rotamer identity
	 * 		for the first rotamer in the rotamer pair; the last three dimensions define the second rotamer in
	 * 		the rotamer pair;
	 * The first and fourth dimensions are of size (numInAS+1) if the input system does not have a ligand, or (numInAS+2) otherwise;
	 * 		The (numInAS+1)^th row in the first/fourth dimensions corresponds to the ligand entries (if present);
	 * 		The last row in the first dimension corresponds to the template energy (the last row in the fourth dimension is never used);
	 * The intra-rotamer and rot-shell energies for a given rotamer identity 'r' for amino acid type 'a' at residue position 'p'
	 * 		are stored in the two entries in [p][a][r][p][0][] (since there are no rotamer pairwise energies for rotamers at the 
	 * 		same residue position);
	 * NOTE: the difference between ligPresent and useLig is the following:
	 * 		ligPresent determines if a ligand is present in the input structure;
	 * 		useLig determines if the ligand will be used in the current computation (e.g., useLig will be false for SHL-AS)
	 */
	public float [][][][][][] initializePairEMatrix(int numInAS, boolean ligPresent, int resMut[], int residueMap[],
			RotamerSearch rs, boolean useLig, String ligType, boolean shellRun, boolean intraRun, RotamerLibrary rl, 
			RotamerLibrary grl, int numAAallowed, boolean initAll) {
		
		float eMatrix[][][][][][] = null;
		int numPos = numInAS + 1;
		if (ligPresent)
			numPos++;
		eMatrix = new float[numPos][][][][][];
		for (int p1=0; p1<numInAS; p1++){
			if (resMut[p1] == 1) {
				eMatrix[p1] = new float[numAAallowed][][][][];
				for (int a1=0; a1<rs.sysLR.getNumAllowable(residueMap[p1]); a1++){
					int curAAind1 = rs.sysLR.getIndexOfNthAllowable(residueMap[p1],a1);
					int numRot1 = rl.getNumRotForAAtype(curAAind1);
					if (numRot1==0) //ALA or GLY
						numRot1 = 1;
					eMatrix[p1][curAAind1] = new float[numRot1][][][];
					for (int r1=0; r1<numRot1; r1++){
						eMatrix[p1][curAAind1][r1] = new float[numPos][][];
						for (int p2=0; p2<numInAS; p2++){
							if ((initAll)||((!intraRun)&&(p2!=p1))){
								eMatrix[p1][curAAind1][r1][p2] = new float[numAAallowed][];
								for (int a2=0; a2<rs.sysLR.getNumAllowable(residueMap[p2]); a2++){
									int curAAind2 = rs.sysLR.getIndexOfNthAllowable(residueMap[p2],a2);
									int numRot2 = rl.getNumRotForAAtype(curAAind2);
									if (numRot2==0) //ALA or GLY
										numRot2 = 1;
									eMatrix[p1][curAAind1][r1][p2][curAAind2] = new float[numRot2];
									for (int r2=0; r2<numRot2; r2++){
										eMatrix[p1][curAAind1][r1][p2][curAAind2][r2] = 0.0f;
									}
								}
							}
							if (p2==p1) {
								eMatrix[p1][curAAind1][r1][p2] = new float[1][2];
							}
						}
						if (useLig){
							int p2 = numInAS;
							eMatrix[p1][curAAind1][r1][p2] = new float[numAAallowed][];
							int curAAind2 = grl.getAARotamerIndex(ligType);
							int numRot2 = grl.getNumRotamers(ligType);
							if (numRot2==0) //ALA or GLY
								numRot2 = 1;
							eMatrix[p1][curAAind1][r1][p2][curAAind2] = new float[numRot2];
							for (int r2=0; r2<numRot2; r2++){
								eMatrix[p1][curAAind1][r1][p2][curAAind2][r2] = 0.0f;
							}
						}
					}
				}
			}
		}
		if (useLig){ //ligand computation will be performed here
			int p1 = numInAS;
			eMatrix[p1] = new float[numAAallowed][][][][];
			int curAAind1 = grl.getAARotamerIndex(ligType);
			int numRot1 = grl.getNumRotamers(ligType);
			if (numRot1==0) //ALA or GLY
				numRot1 = 1;
			eMatrix[p1][curAAind1] = new float[numRot1][][][];
			for (int r1=0; r1<numRot1; r1++){
				eMatrix[p1][curAAind1][r1] = new float[numPos][][];
				for (int p2=0; p2<numInAS; p2++){
					eMatrix[p1][curAAind1][r1][p2] = new float[numAAallowed][];
					for (int a2=0; a2<rs.sysLR.getNumAllowable(residueMap[p2]); a2++){
						int curAAind2 = rs.sysLR.getIndexOfNthAllowable(residueMap[p2],a2);
						int numRot2 = rl.getNumRotForAAtype(curAAind2);
						if (numRot2==0) //ALA or GLY
							numRot2 = 1;
						eMatrix[p1][curAAind1][r1][p2][curAAind2] = new float[numRot2];
						for (int r2=0; r2<numRot2; r2++){
							eMatrix[p1][curAAind1][r1][p2][curAAind2][r2] = 0.0f;
						}
					}
				}
				eMatrix[p1][curAAind1][r1][p1] = new float[1][2];
			}
		}
		if (shellRun){
			eMatrix[numPos-1] = new float[1][1][1][1][1];
		}
		
		return eMatrix;
	}
	
	
	//Returns a new independent six-dimensional matrix that is a copy of fromMatrix[][][][][][]
	public float [][][][][][] copyMultiDimArray(float fromMatrix[][][][][][]){
		
		if (fromMatrix==null)
			return null;
		
		float toMatrix[][][][][][] = new float[fromMatrix.length][][][][][];
		for (int p1=0; p1<toMatrix.length; p1++){
			if (fromMatrix[p1]!=null){
				toMatrix[p1] = new float[fromMatrix[p1].length][][][][];
				for (int a1=0; a1<toMatrix[p1].length; a1++){
					if (fromMatrix[p1][a1]!=null){
						toMatrix[p1][a1] = new float[fromMatrix[p1][a1].length][][][];
						for (int r1=0; r1<toMatrix[p1][a1].length; r1++){
							if (fromMatrix[p1][a1][r1]!=null){
								toMatrix[p1][a1][r1] = new float[fromMatrix[p1][a1][r1].length][][];
								for (int p2=0; p2<toMatrix[p1][a1][r1].length; p2++){
									if (fromMatrix[p1][a1][r1][p2]!=null){
										toMatrix[p1][a1][r1][p2] = new float[fromMatrix[p1][a1][r1][p2].length][];
										for (int a2=0; a2<toMatrix[p1][a1][r1][p2].length; a2++){
											if (fromMatrix[p1][a1][r1][p2][a2]!=null){
												toMatrix[p1][a1][r1][p2][a2] = new float[fromMatrix[p1][a1][r1][p2][a2].length];
												System.arraycopy(fromMatrix[p1][a1][r1][p2][a2], 0, toMatrix[p1][a1][r1][p2][a2], 0, fromMatrix[p1][a1][r1][p2][a2].length);
											}
										}
									}
								}
							}
						}
					}
				}
			}				
		}
		
		return toMatrix;
	}
	
	//Called by slave nodes to generate cObj.compEE[] entries to return to the main node;
	//The two matrices should have the same structure (i.e., a computed entry in one matrix should also be computed in the other)
	public SamplingEEntries [] generateCompEE(float minEmatrix[][][][][][], float maxEmatrix[][][][][][]){
		
		if (minEmatrix==null || maxEmatrix==null) {
			System.out.println("ERROR: cannot generate compEE[] entries from a null PEM matrix.");
			System.exit(1);
		}
		
		else {
			SamplingEEntries compEE[] = new SamplingEEntries[100];
			int curEntry = 0;
			for (int p1=0; p1<minEmatrix.length; p1++){
				if (minEmatrix[p1]!=null){
					for (int a1=0; a1<minEmatrix[p1].length; a1++){
						if (minEmatrix[p1][a1]!=null){
							for (int r1=0; r1<minEmatrix[p1][a1].length; r1++){
								if (minEmatrix[p1][a1][r1]!=null){
									for (int p2=0; p2<minEmatrix[p1][a1][r1].length; p2++){
										if (minEmatrix[p1][a1][r1][p2]!=null){
											for (int a2=0; a2<minEmatrix[p1][a1][r1][p2].length; a2++){
												if (minEmatrix[p1][a1][r1][p2][a2]!=null){
													for (int r2=0; r2<minEmatrix[p1][a1][r1][p2][a2].length; r2++){
														if ( (minEmatrix[p1][a1][r1][p2][a2][r2]!=0.0f) || (maxEmatrix[p1][a1][r1][p2][a2][r2]!=0.0f) ) {
														
															compEE[curEntry] = new SamplingEEntries();
															compEE[curEntry].i1 = p1;
															compEE[curEntry].i2 = a1;
															compEE[curEntry].i3 = r1;
															compEE[curEntry].i4 = p2;
															compEE[curEntry].i5 = a2;
															compEE[curEntry].i6 = r2;
															compEE[curEntry].minE = minEmatrix[p1][a1][r1][p2][a2][r2];
															compEE[curEntry].maxE = maxEmatrix[p1][a1][r1][p2][a2][r2];
															
															curEntry++;
															
															if (curEntry>=compEE.length){
																SamplingEEntries tmp[] = new SamplingEEntries[compEE.length*2];
																System.arraycopy(compEE, 0, tmp, 0, compEE.length);
																compEE = tmp;
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}				
			}
			
			SamplingEEntries tmp[] = new SamplingEEntries[curEntry];
			System.arraycopy(compEE, 0, tmp, 0, curEntry);
			compEE = tmp;
			
			return compEE;
		}
		
		return null;
	}
}
