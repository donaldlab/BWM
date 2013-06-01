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
//SimpleMinimizer.java
//
//Version:           1.0
//
//
//authors:
//    initials    name                 organization                email
//   ---------   -----------------    ------------------------    ----------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//
///////////////////////////////////////////////////////////////////////////////////////////////

/* 
* Written by: Ryan Lilien (2001-2004);
* 		modified by Ivelin Georgiev (2004-2009)
*
* This class implements a simple energy minimization routine for the side-chains only
*  using the Amber96 forcefield, electrostatic, vdw, and dihedral energy
*  terms, and using the assumption that only certain torsions
*  can bend. Additionally there is a special residue, the
*  ligand, that can translate and globally rotate.
* 
* The system consists of one molecule containing: a protein strand, 
 *  a ligand strand (optional), and a cofactor strand (optional).
 * The protein strand does not have to contain sequential residues
 *  but it must be made of standard amino acids
 * The ligand strand can only be one 'thing'
 *  -if this 'thing' is an AA then the Penultimate Rotamer library is used
 *  -if this 'thing' is not an AA then a generic rotamer library is used
* 
*/



import java.io.Serializable;

/**
 * This class implements a simple energy minimization routine for the side-chains only. Handles the computation
 * of the Amber dihedral energy penalties (for minimizing away from the initial rotamer dihedrals).
 * Additionally there is a special residue, the ligand, that can translate and globally rotate.
 */
public class SimpleMinimizer implements Serializable {

	// If the debug flag is set to true then additional debug statements are
	//  printed to standard out.
	boolean debug = false;
	
	private int MAX_NUM_ATOMS_DISTAL = 30; // the maximum number of atoms distal to the dihedral; this is increased as needed for big ligands


	Molecule m = null;
	Amber96ext a96ff = null;
		// the molecule and forcefield (typical stuff)
	int sysStrNum = -1;
	int ligStrNum = -1;
		// the system and ligand strand numbers
		// if ligStrNum == -1 then there is no ligand
	StrandRotamers sysRH = null;
	StrandRotamers ligRH = null;
		// the system and ligand rotamer handlers, these
		//  are required so that we can get a handle on
		//  the dihedrals
	int numSysDihedrals = 0;
	int numLigDihedrals = 0;
		// the number of system and ligand dihedrals
	int sysDihedralAtNums[][] = null;
		// array of dihedrals, this array is n by 4
		//  and has the moleculeAtomNumbers for each
		//  of the 4 atoms of the dihedral
	int sysDihedralDistal[][] = null;
		// array of dihedrals, the first index is over all
		//  dihedrals, each row consists of an atomList
		//  that can be passed to the setTorsion
	int sysNumAtomsDistal[] = null;
		// the number of atoms distal for each dihedral,
		//  ie. the first x entries in a row of 
		//  sysDihedralDistal are valid atoms
	int ligDihedralAtNums[][] = null;
		// array of dihedrals, this array is n by 4
		//  and has the moleculeAtomNumbers for each
		//  of the 4 atoms of the dihedral
	int ligDihedralDistal[][] = null;
		// similar to sysDihedralArray but for the
		//  dihedrals of the ligand residue
	int ligNumAtomsDistal[] = null;
		// the number of atoms distal for each dihedral,
		//  ie. the first x entries in a row of 
		//  ligDihedralDistal are valid atoms

	float initialAngleStepSize = 1.0f;
		// initial angular step size (in degrees)
	double sysDihedDiff[] = null;
		// current finite difference step for the given dihedral
	double ligDihedDiff[] = null;
		// current finite difference step for the given dihedral
	float tempCoords[] = null;
		// temporary holding location for coordinates as we
		//  muck with a dihedral, note that the size of this
		//  array is defined by MAX_NUM_ATOMS_DISTAL * 3
	double sysCumulativeDihedStep[] = null;
		// the cumulative change in angle taken by the given
		//  dihedral. this is kept track of because we limit
		//  the total movement
	double ligCumulativeDihedStep[] = null;
		// the cumulative change in angle taken by the given
		//  dihedral. this is kept track of because we limit
		//  the total movement
	double ligStartCOM[] = new double[3];
		// the position of the ligand's initial COM
	double ligCurCOM[] = new double[3];
		// the ligand's current COM

	double maxMovement = 9.0;
		// maximum degrees by which a torsion can
		//  cumulatively change
	float ligRotStep = 0.25f;
		// step size for the rigid rotation
		//  of the ligand
	float ligTransStep = 0.10f;
		// step size in � for the rigid ligand
		//  translation
	double ligMaxTrans = 1.2;
		// the maximum ligand translation allowed
	
	int numFlexRes = -1;
		// the number of flexible residues
	int flexResAtomList[][] = new int[0][0];
		// each row is a residue and contains
		//  the moleucleatomnumbers located more
		//  distal than the 3rd atom of any dihedral
	int flexResListSize[] = new int[0];
		// the number of valid elements in each
		//  row of flexResAtomList
	int sysDihedToResNum[] = new int[0];
		// mapping from dihedral number to residue
		//  number
	int ligResNum = -1;
		// the residue index of the ligand in the
		//  flexResAtomList, flexResListSize
	
	// Section for dealing with dihedral energies
	// I know this shouldn't be in the minimizer but
	//  there's no other good way to do this.
	int sDihedPN[] = new int[0];
	double sDihedWeight[] = new double[0];
	int sDihedLocalNum[] = new int[0];
	int sDihedNumTerms = 0;
		// each of these arrays contains the PN and Weight for
		//  a dihedral term, each physical dihedral may contain
		//  multiple 'paths' (ie. set of four atoms that goes
		//  through the same central two). And each path may
		//  contain multiple terms (due to AMBER).
		// the local number is the number of the dihedral in
		//  the sysCumulativeDihedStep array
		// the dihednumterms is the total number of terms
	int lDihedPN[] = new int[0];
	double lDihedWeight[] = new double[0];
	int lDihedLocalNum[] = new int[0];
	int lDihedNumTerms = 0;
		// the same as the previous 4 variables, but for the
		//  ligand.
	boolean doDihedEnergy = false;
	// If true then dihedral energies are computed and
	//  used in minimization

	
	SimpleMinimizer(){
		// not much to do here
	
	}
	
		
	// Initializes the Minimizer for use with a system without a ligand
	public void initialize(Molecule theM, int theSysStrNum, Amber96ext theA96ff,
		StrandRotamers theSysRH, int curAANum[], boolean doDihedral, RotamerLibrary rl){
	
		// snag the local variables
		m = theM;
		a96ff = theA96ff;
		doDihedEnergy = doDihedral;
		sysRH = theSysRH;
		sysStrNum = theSysStrNum;
	
		ligStrNum = -1;
		numLigDihedrals = 0;
		ligResNum = -1;
		for(int i=0;i<3;i++){
			ligCurCOM[i] = 0.0;
			ligStartCOM[i] = 0.0;
		}
			
		numSysDihedrals = 0;		
		// setup system dihedrals
		// First determine how many there are
		for(int i=0;i<m.strand[sysStrNum].numberOfResidues;i++){
			if(m.strand[sysStrNum].residue[i].flexible)
				numSysDihedrals += rl.getNumDihedrals(curAANum[i]);
		}		
		
		sysDihedralAtNums = new int[numSysDihedrals][4];
		sysDihedralDistal = new int[numSysDihedrals][MAX_NUM_ATOMS_DISTAL];
		sysNumAtomsDistal = new int[numSysDihedrals];
		tempCoords = new float[MAX_NUM_ATOMS_DISTAL * 3];

		// Count number of flexible residues
		numFlexRes = 0;
		for(int i=0;i<m.strand[sysStrNum].numberOfResidues;i++){
			if(m.strand[sysStrNum].residue[i].flexible)
				numFlexRes++;
		}
		
		// 2 is added to numFlexRes so that there is room for the ligand at the
		//  end if there is a ligand present, if there is no ligand then it
		//  doesn't really matter as we won't look at that row
		// The first ligand term includes nonbonded terms for computing dihedrals
		// The second ligand term includes nonbonded terms for computing the gradient
		//  and thus includes terms for all atoms
		flexResAtomList = new int[numFlexRes+2][MAX_NUM_ATOMS_DISTAL];
		flexResListSize = new int[numFlexRes+2];
		sysDihedToResNum = new int[numSysDihedrals];
		
		// Now determine which atoms are involved with each system dihedral
		int atoms[] = new int[4];
		int tmpAtomList[] = new int[MAX_NUM_ATOMS_DISTAL];
		int at3num = -1, at2num = -1;
		int numDihed = 0;
		int curNumFlex = 0;
		Residue localRes = null;
		for(int i=0;i<m.strand[sysStrNum].numberOfResidues;i++){
			localRes = m.strand[sysStrNum].residue[i];
			if(localRes.flexible){
				for(int j=0;j<rl.getNumDihedrals(curAANum[i]);j++){
					sysDihedToResNum[numDihed] = curNumFlex;
					atoms = rl.getDihedralInfo(m,sysStrNum,i,j);
					// note: atoms are residue relative numbering
					sysDihedralAtNums[numDihed][0] = localRes.atom[atoms[0]].moleculeAtomNumber;
					sysDihedralAtNums[numDihed][1] = localRes.atom[atoms[1]].moleculeAtomNumber;
					sysDihedralAtNums[numDihed][2] = localRes.atom[atoms[2]].moleculeAtomNumber;
					sysDihedralAtNums[numDihed][3] = localRes.atom[atoms[3]].moleculeAtomNumber;
					
					at2num = localRes.atom[atoms[2]].moleculeAtomNumber;
					at3num = localRes.atom[atoms[3]].moleculeAtomNumber;					
					for(int k=0;k<MAX_NUM_ATOMS_DISTAL;k++){
						tmpAtomList[k]=1;
					}
					tmpAtomList[atoms[1]]=0;
					tmpAtomList[atoms[2]]=0;
					sysRH.getAtomsMoreDistal(m,localRes.moleculeResidueNumber,m.atom[at2num],tmpAtomList);
					sysNumAtomsDistal[numDihed]=0;
					for(int k=0;k<MAX_NUM_ATOMS_DISTAL;k++){
						if ((tmpAtomList[k]==2) && (k != atoms[3])){
							sysDihedralDistal[numDihed][sysNumAtomsDistal[numDihed]]=localRes.atom[k].moleculeAtomNumber;
							sysNumAtomsDistal[numDihed] += 1;
						}
					}
					if (j==0){
						// If this is the first dihedral for this residue snag information
						flexResListSize[curNumFlex] = sysNumAtomsDistal[numDihed]+1;
						flexResAtomList[curNumFlex][0] = at3num;
						for(int k=0;k<sysNumAtomsDistal[numDihed];k++){
							flexResAtomList[curNumFlex][k+1] = sysDihedralDistal[numDihed][k];
						}
					}
					numDihed++;
				}
				curNumFlex++;
			}
		}
	}	
	
	// Initializes the Minimizer for use with a system and a ligand
	public void initialize(Molecule theM, int theSysStrNum, int theLigStrNum,
		Amber96ext theA96ff, StrandRotamers theSysRH, StrandRotamers theLigRH, int
		curAANum[], int curLigAANum, boolean doDihedral, RotamerLibrary rl, RotamerLibrary grl){
		
		MAX_NUM_ATOMS_DISTAL = Math.max(MAX_NUM_ATOMS_DISTAL, theM.strand[theLigStrNum].residue[0].numberOfAtoms);
	
		// First call the other initialize interface to grab the
		//  system dihedrals
		initialize(theM,theSysStrNum,theA96ff,theSysRH,curAANum,doDihedral,rl);
		doDihedEnergy = doDihedral;
		
		ligRH = theLigRH;
		ligStrNum = theLigStrNum;

		// Now compute the ligand dihedrals
		// First determine how many there are
		numLigDihedrals = grl.getNumDihedrals(curLigAANum);
		
		ligDihedralAtNums = new int[numLigDihedrals][4];
		ligDihedralDistal = new int[numLigDihedrals][MAX_NUM_ATOMS_DISTAL];
		ligNumAtomsDistal = new int[numLigDihedrals];

		// numFlexRes, flexResAtomList, flexResListSize are setup in the system
		//  initialize code. The ligand information at numFlexRes are the dihedral
		//  relevant nonbonded terms, the ligand information at numFlexRes+1 are
		//  used in computing the gradient and include all atoms of the ligand not
		//  just those distal to the first dihedral
	
		// Now determine which atoms are involved with each dihedral
		int atoms[] = new int[4];
		int tmpAtomList[] = new int[MAX_NUM_ATOMS_DISTAL];
		int at3num = -1, at2num = -1;
		int numDihed = 0;
		Residue localRes = m.strand[ligStrNum].residue[0];
		// Snag all ligand atom information for flexResAtomList, flexResListSize
		ligResNum = numFlexRes;  // numFlexRes was set in the other initialize call
		int ligResNumPP = ligResNum + 1;
		flexResListSize[ligResNumPP] = localRes.numberOfAtoms;
		for(int k=0;k<flexResListSize[ligResNumPP];k++){
			flexResAtomList[ligResNumPP][k] = localRes.atom[k].moleculeAtomNumber;
		}

		if(m.strand[ligStrNum].residue[0].flexible){
			for(int j=0;j<grl.getNumDihedrals(curLigAANum);j++){
				// ligDihedToResNum[numDihed] = 0; // all ligand dihedrals are from residue 0
				atoms = grl.getDihedralInfo(m,ligStrNum,0,j);
				// note: atoms are residue relative numbering
				ligDihedralAtNums[numDihed][0] = localRes.atom[atoms[0]].moleculeAtomNumber;
				ligDihedralAtNums[numDihed][1] = localRes.atom[atoms[1]].moleculeAtomNumber;
				ligDihedralAtNums[numDihed][2] = localRes.atom[atoms[2]].moleculeAtomNumber;
				ligDihedralAtNums[numDihed][3] = localRes.atom[atoms[3]].moleculeAtomNumber;
				
				at2num = localRes.atom[atoms[2]].moleculeAtomNumber;
				at3num = localRes.atom[atoms[3]].moleculeAtomNumber;
					
				for(int k=0;k<MAX_NUM_ATOMS_DISTAL;k++){
					tmpAtomList[k]=1;
				}
				tmpAtomList[atoms[1]]=0;
				tmpAtomList[atoms[2]]=0;
				ligRH.getAtomsMoreDistal(m,localRes.moleculeResidueNumber,m.atom[at2num],tmpAtomList);
				ligNumAtomsDistal[numDihed]=0;
				for(int k=0;k<MAX_NUM_ATOMS_DISTAL;k++){
					if ((tmpAtomList[k]==2) && (k != atoms[3])){
						ligDihedralDistal[numDihed][ligNumAtomsDistal[numDihed]]=localRes.atom[k].moleculeAtomNumber;
						ligNumAtomsDistal[numDihed] += 1;
					}
				}
				if (j==0){
					// If this is the first dihedral for this residue snag information
					flexResListSize[ligResNum] = ligNumAtomsDistal[numDihed]+1;
					flexResAtomList[ligResNum][0] = at3num;
					for(int k=0;k<ligNumAtomsDistal[numDihed];k++){
						flexResAtomList[ligResNum][k+1] = ligDihedralDistal[numDihed][k];
					}
				}
				numDihed++;
			}
		}
	}
	
	public int getNumTotDihed(){
		return numSysDihedrals + numLigDihedrals;
	}
	
	// Sets the initial angle step size
	public void setInitialAngleStepSize(float num){	
		initialAngleStepSize = num;
	}
	
	// Setup the dihedral angle terms for both the system and ligand
	public boolean setupDihedralTerms(){
		
		// At this point we assume sysDihedralAtNums, ligDihedralAtNums,
		//  numSysDihedrals, and numLigDihedrals are all correct.
		
		int tmpSSize = 1000;
		int tmpLSize = 1000;
		sDihedPN = new int[tmpSSize];
		lDihedPN = new int[tmpLSize];
		sDihedWeight = new double[tmpSSize];
		lDihedWeight = new double[tmpLSize];
		sDihedLocalNum = new int[tmpSSize];
		lDihedLocalNum = new int[tmpLSize];
		sDihedNumTerms = 0;
		lDihedNumTerms = 0;
		int atom2 = 0, atom3 = 0;
		double fConst[] = new double[5];
		double eqAngle[] = new double[5];
		int terms[] = new int[5];
		int multiplic[] = new int[1];
		
		int newDihedPN[] = null;
		double newDihedWeight[] = null;
		int newDihedLocalNum[] = null;
		
		// SYSTEM DIHEDRALS
		for(int i=0;i<numSysDihedrals;i++) {
			// Get center two atoms
			atom2 = sysDihedralAtNums[i][1];
			atom3 = sysDihedralAtNums[i][2];
			for(int q=0;q<m.connected[atom2][0];q++) {
				if (m.connected[atom2][q+1] != atom3) {
					for(int w=0;w<m.connected[atom3][0];w++) {
						if (m.connected[atom3][w+1] != atom2) {
							// At this point 'q',atom2,atom3,'w' is a dihedral
							if (!(a96ff.getTorsionParameters(m.atom[m.connected[atom2][q+1]].type,m.atom[atom2].type,
															 m.atom[atom3].type,m.atom[m.connected[atom3][w+1]].type,
															 fConst,eqAngle,terms,multiplic))) {
								System.out.println("WARNING: Could not find torsion parameters for " + m.connected[atom2][q+1] + "," + atom2 + "," + atom3 +"," + m.connected[atom3][w+1]);
								return false;
							}
							for(int r=0;r<=multiplic[0];r++) {
								sDihedNumTerms++;
								if (sDihedNumTerms > tmpSSize) {
									// increase the size of the arrays
									newDihedPN = new int[tmpSSize + 500];
									newDihedWeight = new double[tmpSSize + 500];
									newDihedLocalNum = new int[tmpSSize + 500];
									System.arraycopy(sDihedPN,0,newDihedPN,0,tmpSSize);
									System.arraycopy(sDihedWeight,0,newDihedWeight,0,tmpSSize);
									System.arraycopy(sDihedLocalNum,0,newDihedLocalNum,0,tmpSSize);
									sDihedPN = newDihedPN;
									sDihedWeight = newDihedWeight;
									sDihedLocalNum = newDihedLocalNum;
									tmpSSize += 500;
								}
								sDihedPN[sDihedNumTerms] = terms[r];
								sDihedWeight[sDihedNumTerms] = fConst[r];
								sDihedLocalNum[sDihedNumTerms] = i;
							}
						}
					}
				}
			}
		}
		
		// Shrink the size of the array
		newDihedPN = new int[sDihedNumTerms];
		newDihedWeight = new double[sDihedNumTerms];
		newDihedLocalNum = new int[sDihedNumTerms];
		System.arraycopy(sDihedPN,0,newDihedPN,0,sDihedNumTerms);
		System.arraycopy(sDihedWeight,0,newDihedWeight,0,sDihedNumTerms);
		System.arraycopy(sDihedLocalNum,0,newDihedLocalNum,0,sDihedNumTerms);
		sDihedPN = newDihedPN;
		sDihedWeight = newDihedWeight;
		sDihedLocalNum = newDihedLocalNum;
		
		// LIGAND DIHEDRALS
		for(int i=0;i<numLigDihedrals;i++) {
			// Get center two atoms
			atom2 = ligDihedralAtNums[i][1];
			atom3 = ligDihedralAtNums[i][2];
			for(int q=0;q<m.connected[atom2][0];q++) {
				if (m.connected[atom2][q+1] != atom3) {
					for(int w=0;w<m.connected[atom3][0];w++) {
						if (m.connected[atom3][w+1] != atom2) {
							// At this point 'q',atom2,atom3,'w' is a dihedral
							if (!(a96ff.getTorsionParameters(m.atom[m.connected[atom2][q+1]].type,m.atom[atom2].type,
															 m.atom[atom3].type,m.atom[m.connected[atom3][w+1]].type,
															 fConst,eqAngle,terms,multiplic))) {
								System.out.println("WARNING: Could not find torsion parameters for " + m.connected[atom2][q+1] + "," + atom2 + "," + atom3 +"," + m.connected[atom3][w+1]);
								return false;
							}
							for(int r=0;r<=multiplic[0];r++) {
								lDihedNumTerms++;
								if (lDihedNumTerms > tmpLSize) {
									// increase the size of the arrays
									newDihedPN = new int[tmpLSize + 500];
									newDihedWeight = new double[tmpLSize + 500];
									newDihedLocalNum = new int[tmpLSize + 500];
									System.arraycopy(lDihedPN,0,newDihedPN,0,tmpLSize);
									System.arraycopy(lDihedWeight,0,newDihedWeight,0,tmpLSize);
									System.arraycopy(lDihedLocalNum,0,newDihedLocalNum,0,tmpLSize);
									lDihedPN = newDihedPN;
									lDihedWeight = newDihedWeight;
									lDihedLocalNum = newDihedLocalNum;
									tmpLSize += 500;
								}
								lDihedPN[lDihedNumTerms] = terms[r];
								lDihedWeight[lDihedNumTerms] = fConst[r];
								lDihedLocalNum[lDihedNumTerms] = i;
							}
						}
					}
				}
			}
		}
		
		// Shrink the size of the array
		newDihedPN = new int[lDihedNumTerms];
		newDihedWeight = new double[lDihedNumTerms];
		newDihedLocalNum = new int[lDihedNumTerms];
		System.arraycopy(lDihedPN,0,newDihedPN,0,lDihedNumTerms);
		System.arraycopy(lDihedWeight,0,newDihedWeight,0,lDihedNumTerms);
		System.arraycopy(lDihedLocalNum,0,newDihedLocalNum,0,lDihedNumTerms);
		lDihedPN = newDihedPN;
		lDihedWeight = newDihedWeight;
		lDihedLocalNum = newDihedLocalNum;

		return true;
	}
	
	
	// Computes and returns the dihedral energy for the dihedrals
	//  that are allowed to change during minimization (both for
	//  the system and the ligand). Each dihedral is assumed to
	//  start in a energy well of the appropriate dihedral energy
	//  terms from the AMBER forcefield.
	// Should be called only after calling minimize()
	public double computeDihedEnergy() {
	
		double dihedE = 0.0;
		double d2r = 3.14159265 / 180.0;

		// System dihedral terms first
		for(int i=0;i<sDihedNumTerms;i++){
			dihedE += sDihedWeight[i] * (1 - Math.cos(sDihedPN[i]*sysCumulativeDihedStep[sDihedLocalNum[i]]*d2r));
		}
		// Ligand dihedral terms
		for(int i=0;i<lDihedNumTerms;i++){
			dihedE += lDihedWeight[i] * (1 - Math.cos(lDihedPN[i]*ligCumulativeDihedStep[lDihedLocalNum[i]]*d2r));
		}
		
		//System.out.println("computeDihedEnergy "+dihedE);
		
		return(dihedE);
	}
	
	// This function sets up the dihedral energy computation for the
	//	given dihedral; the computation is done in computeOneDihedEnergyDiffHelper()
	public double computeOneDihedEnergyDiff(boolean isLigDihed, int dihedNumForCur, double dihedChange){
		
		int sDNum, lDNum;
		double sDChange, lDChange;
		if (!isLigDihed) { //system dihedral
			sDNum = dihedNumForCur;
			sDChange = dihedChange;
			lDNum = -1;
			lDChange = 0.0;
		}
		else { //ligand dihedral
			sDNum = -1;
			sDChange = 0.0;
			lDNum = dihedNumForCur;
			lDChange = dihedChange;
		}
		
		return computeOneDihedEnergyDiffHelper(sDNum, sDChange, lDNum, lDChange);
	}	

	// Computes the difference in energy for specified dihedral
	//  initial - change. 
	// System dihedral sDNum is changed by sDChange and ligand
	//  dihedral lDNum is changed by lDChange. These changes are
	//  only temporary and used to compute an energy, that is the
	//  sys(lig)CumulativeDihedStep is not changed. If sDNum or
	//  lDNum are -1 then they are not changed.
	public double computeOneDihedEnergyDiffHelper(int sDNum, double sDChange,
										int lDNum, double lDChange) {
		
		double dihedE = 0.0;
		double d2r = 3.14159265 / 180.0;
		
		// System dihedral terms first
		if (sDNum >= 0) {
			for(int i=0;i<sDihedNumTerms;i++) {
				if (sDihedLocalNum[i] == sDNum)
					dihedE -= sDihedWeight[i] * (1 - Math.cos(sDihedPN[i]*(sDChange)*d2r));
			}
		}
		// Ligand dihedral terms
		if (lDNum >= 0) {
			for(int i=0;i<lDihedNumTerms;i++) {
				if (lDihedLocalNum[i] == lDNum)
					dihedE -= lDihedWeight[i] * (1 - Math.cos(lDihedPN[i]*(lDChange)*d2r));
			}
		}

		return(dihedE);
	}

	// Returns the cross product: a x b
	public double[] crossProduct(double a[], double b[]){

		double f[] = new double[3];
		f[0] = a[1]*b[2] - a[2]*b[1];
		f[1] = a[2]*b[0] - a[0]*b[2];
		f[2] = a[0]*b[1] - a[1]*b[0];
		return(f);
	}
	
	// Returns the dot product of a and b
	public double dotProduct(double a[], double b[]){
		return(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
	}

	public double getMagnitude(double a[]){
		double sum = 0.0;
		for(int i=0;i<a.length;i++)
			sum += (a[i]*a[i]);
		return(Math.sqrt(sum));
	}
	
	// Computes the torque around the specified dihedral
	public double[] computeTorque(int diArray[]){
	
		// remember that the indices of the dihedralArray are
		// [0] residueNum
		// [1] numberOfAtoms (n)
		// [2] atom1
		// [3] atom2
		// [4] atom3
		// [5] atom4 (the distal atom of the dihedral)
		// ...
		// [n+1] atomn
		// the result is a float array where the first
		//  three values are the axis and the fourth
		//  is the torque
		// interestingly since our force is kcal/mol/A
		//  the torque is kcal*A/mol/A or kcal/mol
		
		double result[] = new double[4];
		double torque[] = new double[3];
		double temp[] = new double[3];
		double R[] = new double[3];
		double F[] = new double[3];
		int atBase = -1;
		Residue localRes = m.residue[diArray[0]];
		int at3Base = localRes.atom[diArray[4]].moleculeAtomNumber * 3;
		int at2Base = localRes.atom[diArray[3]].moleculeAtomNumber * 3;
		
		// compute the rotation (dihedral) axis
		result[0] = m.actualCoordinates[at3Base] - m.actualCoordinates[at2Base];
		result[1] = m.actualCoordinates[at3Base+1] - m.actualCoordinates[at2Base+1];
		result[2] = m.actualCoordinates[at3Base+2] - m.actualCoordinates[at2Base+2];
		// normalize the axis to unit length
		double axisMag = getMagnitude(result);
		result[0] /= axisMag;
		result[1] /= axisMag;
		result[2] /= axisMag;
		
		torque[0] = 0.0;
		torque[1] = 0.0;
		torque[2] = 0.0;

		// loop over all involved atoms
		for(int i=5;i<diArray[1]+2;i++){
			atBase = localRes.atom[diArray[i]].moleculeAtomNumber * 3;
			// R is the vector from atom3 to the atom of interest
			R[0] = m.actualCoordinates[atBase] - m.actualCoordinates[at3Base];
			R[1] = m.actualCoordinates[atBase+1] - m.actualCoordinates[at3Base+1];
			R[2] = m.actualCoordinates[atBase+2] - m.actualCoordinates[at3Base+2];
			// F is the current gradient on the atom of interest
			F[0] = m.gradient[atBase];
			F[1] = m.gradient[atBase+1];
			F[2] = m.gradient[atBase+2];
			temp = crossProduct(R,F);
			torque[0] += temp[0];
			torque[1] += temp[1];
			torque[2] += temp[2];
		}
		// Project the sum of the torques onto the axis of rotation
		//  which is stored in the first 3 elements of result
		result[3] = dotProduct(torque,result);

		return(result);
	}
	
	//Sets up and performs ligand translation/rotation
	private void doLigTransRot(double ligTorque[], double ligTrans[], double lligRotStep, double lligTransStep, double lligMaxTrans){
		
		float bckpLigCoords[] = backupLigCoord(); //backup actual ligand coordinates
		double initE[] = a96ff.calculateTotalEnergy(m.actualCoordinates, -1); //compute energy before translation/rotation
	
		a96ff.calculateGradient(-1); //calculate the gradient (perhaps eventually just calculate part of the gradient)
		//a96ff.calculateEVGradientPartWithArrays(ligResNumPP);
		computeLigTorqueTrans(ligTorque, ligTrans);
		applyLigTorqueTrans(ligTorque, lligRotStep, ligTrans, lligTransStep, lligMaxTrans);
		
		double minE[] = a96ff.calculateTotalEnergy(m.actualCoordinates, -1); //compute energy after translation/rotation
		if (initE[0]<minE[0]) { //restore initial ligand actual coordinates if trans/rot increses energy
			restoreLigCoord(bckpLigCoords);
		}
	}

	// Computes the torque around the ligand
	public void computeLigTorqueTrans(double ligTorque[], double ligTrans[]){
	
		Strand thisStrand = m.strand[ligStrNum];
		double centOfMass[] = m.getStrandCOM(ligStrNum);
		int atBase = -1;
		double R[] = new double[3];
		double F[] = new double[3];
		double temp[];
		ligTorque[0] = 0.0;
		ligTorque[1] = 0.0;
		ligTorque[2] = 0.0;
		ligTrans[0] = 0.0;
		ligTrans[1] = 0.0;
		ligTrans[2] = 0.0;

		// loop over all involved atoms
		for(int i=0;i<thisStrand.numberOfResidues;i++){
			for(int j=0;j<thisStrand.residue[i].numberOfAtoms;j++){
				atBase = thisStrand.residue[i].atom[j].moleculeAtomNumber * 3;
				// R is the vector from the COM to the current atom
				R[0] = m.actualCoordinates[atBase] - centOfMass[0];
				R[1] = m.actualCoordinates[atBase+1] - centOfMass[1];
				R[2] = m.actualCoordinates[atBase+2] - centOfMass[2];
				// F is the current gradient for the atom of interest
				F[0] = m.gradient[atBase];
				F[1] = m.gradient[atBase+1];
				F[2] = m.gradient[atBase+2];
				temp = crossProduct(R,F);
				ligTorque[0] += temp[0];
				ligTorque[1] += temp[1];
				ligTorque[2] += temp[2];
				ligTrans[0] += F[0];
				ligTrans[1] += F[1];
				ligTrans[2] += F[2];
			}
		}

	}


	// Applies the specified torque and rigid body translation to the
	//  atoms in the ligand
	// The total translation can not be more than lligMaxTrans
	// The torque is applied with a magnitude of lligRotStep
	private void applyLigTorqueTrans(double ligTorque[], double lligRotStep,
		double ligTrans[], double lligTransStep, double lligMaxTrans){
	
		// normalize the ligand translation
		// ligTrans will now have magnitude lligTransStep
		// we change sign because the gradient points uphill
		double scale = getMagnitude(ligTrans)/lligTransStep;
		if (scale != 0.0) {
			ligTrans[0] /= -scale;
			ligTrans[1] /= -scale;
			ligTrans[2] /= -scale;
		}
		
		// determine the new COM if we took this step
		double tempCOM[] = new double[3];
		tempCOM[0] = ligCurCOM[0] + ligTrans[0];
		tempCOM[1] = ligCurCOM[1] + ligTrans[1];
		tempCOM[2] = ligCurCOM[2] + ligTrans[2];
		
		// compute how large of a step this would be from the
		//  start point
		double totalMovement [] = new double[3];
		totalMovement[0] = tempCOM[0] - ligStartCOM[0];
		totalMovement[1] = tempCOM[1] - ligStartCOM[1];
		totalMovement[2] = tempCOM[2] - ligStartCOM[2];
		scale = getMagnitude(totalMovement);
		
		// if the step would take us too far away then
		//  scale it back
		if (scale > lligMaxTrans){
			scale = scale / lligMaxTrans;
			totalMovement[0] /= scale;
			totalMovement[1] /= scale;
			totalMovement[2] /= scale;
		}
		
		// compute the translation to get us to the new COM
		float theTranslation[] = new float[3];
		theTranslation[0] = (float)(ligStartCOM[0] + totalMovement[0] - ligCurCOM[0]);
		theTranslation[1] = (float)(ligStartCOM[1] + totalMovement[1] - ligCurCOM[1]);
		theTranslation[2] = (float)(ligStartCOM[2] + totalMovement[2] - ligCurCOM[2]);

		// update the current COM
		ligCurCOM[0] = ligStartCOM[0]+totalMovement[0];
		ligCurCOM[1] = ligStartCOM[1]+totalMovement[1];
		ligCurCOM[2] = ligStartCOM[2]+totalMovement[2];
		
		// Apply torque
		// Since our step size is fixed we just need to know the
		//  vector to rotate around
		// we change sign because the gradient points uphill
		double mag = getMagnitude(ligTorque);
		if (mag != 0.0) {
			ligTorque[0] /= -mag;
			ligTorque[1] /= -mag;
			ligTorque[2] /= -mag;
			// The last parameter says not to update the atom coordintes,
			//  only update the actualCoords
			if (Math.abs(lligRotStep)>0.0001)   // this should be in
				m.rotateStrandAroundCOM(ligStrNum,(float)ligTorque[0],(float)ligTorque[1],(float)ligTorque[2],(float)lligRotStep,false);
		}
		// Apply translation, don't update atom coordinates
		// The last parameter says not to update the atom coordintes,
		//  only update the actualCoords
		m.translateStrand(ligStrNum,theTranslation[0],theTranslation[1],theTranslation[2],false);
		
	}

	private void storeCoord(int index, int mAtNum){
	
		tempCoords[index*3] = m.actualCoordinates[mAtNum*3];
		tempCoords[index*3 + 1] = m.actualCoordinates[mAtNum*3 + 1];
		tempCoords[index*3 + 2] = m.actualCoordinates[mAtNum*3 + 2];
	}

	private void restoreCoord(int index, int mAtNum){
		
		m.actualCoordinates[mAtNum*3] = tempCoords[index*3];
		m.actualCoordinates[mAtNum*3 + 1] = tempCoords[index*3 + 1];
		m.actualCoordinates[mAtNum*3 + 2] = tempCoords[index*3 + 2];
	}
	
	private void applyDihedStep(int diAtArray[], int atomList[], int alSize,
		double dihedDiff){
		
		// Perform the rotation
		// at1, at2, and at3 don't move, at4 and the atoms
		//  in the atomList rotate
		// The last parameter says not to update the atom coordintes,
		//  only update the actualCoords
		if (Math.abs(dihedDiff)>0.0001)
			m.changeTorsion(diAtArray[0],diAtArray[1],diAtArray[2],
				diAtArray[3],dihedDiff,atomList,alSize,false);
	}	
	
////////////////////////////////////////////////////////////////////////////////
//  Minimization Section
////////////////////////////////////////////////////////////////////////////////

	// For the given dihedral, compute initial energy, save coordinates,
	//  modify dihedral, compute new energy, compute difference in
	//  energy, restore coords, return step size and direction of lower energy
	private float computeDihedDiff(int diAtArray[], int atomList[],
		int alSize, int AANum, float stepSize, boolean isLigDihed, int dihedNumForCur){
		
		double initialEnergy[], secondEnergy[];
		
		double cumulStep = 0.0;
		if (!isLigDihed)
			cumulStep = sysCumulativeDihedStep[dihedNumForCur];
		else
			cumulStep = ligCumulativeDihedStep[dihedNumForCur];
		
		// Store coordinates
		storeCoord(0,diAtArray[3]);
		for(int i=0;i<alSize;i++)
			storeCoord(i+1,atomList[i]);
		
		// Compute first partial energy
		initialEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,AANum);
		if (doDihedEnergy)
			initialEnergy[0] += computeOneDihedEnergyDiff(isLigDihed,dihedNumForCur,0.0f);
		
		// Apply a rotation
		// at1, at2, and at3 don't move, at4 and the atoms
		//  in the atomList rotate
		m.changeTorsion(diAtArray[0],diAtArray[1],diAtArray[2],
			diAtArray[3],stepSize,atomList,alSize,false);
			
		// Compute second energy
		secondEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,AANum);
		if (doDihedEnergy){
			secondEnergy[0] += computeOneDihedEnergyDiff(isLigDihed,dihedNumForCur,cumulStep+stepSize);
		}
		
		// Restore coordinates
		restoreCoord(0,diAtArray[3]);
		for(int i=0;i<alSize;i++)
			restoreCoord(i+1,atomList[i]);
		
		
		double thirdEnergy[];
		//Store coordinates
		storeCoord(0,diAtArray[3]);
		for(int i=0;i<alSize;i++)
			storeCoord(i+1,atomList[i]);
		
		// Apply a rotation
		// at1, at2, and at3 don't move, at4 and the atoms
		//  in the atomList rotate
		m.changeTorsion(diAtArray[0],diAtArray[1],diAtArray[2],
			diAtArray[3],-stepSize,atomList,alSize,false);
			
		// Compute second energy
		thirdEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,AANum);
		if (doDihedEnergy)
			thirdEnergy[0] += computeOneDihedEnergyDiff(isLigDihed,dihedNumForCur,cumulStep-stepSize);
		
		//Restore coordinates
		restoreCoord(0,diAtArray[3]);
		for(int i=0;i<alSize;i++)
			restoreCoord(i+1,atomList[i]);
		
		if ((initialEnergy[0] > secondEnergy[0])&&(initialEnergy[0] > thirdEnergy[0])){
			if ((initialEnergy[0] - secondEnergy[0])>(initialEnergy[0] - thirdEnergy[0]))
				return stepSize;
			else
				return -stepSize;
		}
		else if(initialEnergy[0] > secondEnergy[0])
			return stepSize;
		else if(initialEnergy[0] > thirdEnergy[0])
			return -stepSize;
		else
			return 0.0f;
	}

	// This function updates the total amount (degrees) in which the given dihedral
	// has moved. It limits this movement to be within +/- lmaxMovement
	private void updateCumulative(double cumulativeDihedStep[], double newDihedDiff[], int index, 
		double lmaxMovement){
		
		if ((cumulativeDihedStep[index] + newDihedDiff[index]) > lmaxMovement)
			newDihedDiff[index] = lmaxMovement - cumulativeDihedStep[index];
		if ((cumulativeDihedStep[index] + newDihedDiff[index]) < -lmaxMovement)
			newDihedDiff[index] = -lmaxMovement - cumulativeDihedStep[index];
		cumulativeDihedStep[index] += newDihedDiff[index];		
	}
	
	//Clears the molecule gradient
	private void clearMolGradient() {
		int natomsx3 = m.numberOfAtoms * 3;
		m.gradient = new double[natomsx3];
		for(int i=0; i<natomsx3; i++){
			m.gradient[i] = 0;
		}
	}
	
	// Performs a simple steepest descent minimization
	// Assumptions:
	//  -a96ff is current
	//     all atoms of appropriate residues have been assigned
	//     initializeEVCalculation has been called
	//  -m.actualCoords contains current atomic coordinates
	// The user should have set numMinimizationSteps and
	//  initialAngleStepSize (else defaults will be used)
	// Uses specific precomputed nonbonded arrays for each residue (this makes things
	//  run faster)
	// If ligandOnly is true, then only the ligand is allowed to minimize, while
	//		the system residues are fixed
	public void minimize(int numSteps, boolean ligandOnly){
	
		if ((ligandOnly)&&(ligStrNum == -1)){
			return;
		}
	
		float step = initialAngleStepSize;
		double lmaxMovement = maxMovement;
			// maximum degrees by which a torsion can
			//  cumulatively change
		float lligRotStep = ligRotStep;
			// step size for the rigid rotation
			//  of the ligand
		float lligTransStep = ligTransStep;
			// step size in � for the rigid ligand
			//  translation
		double lligMaxTrans = ligMaxTrans;
			// the maximum ligand translation allowed
		
		int i=0;
		boolean done = false;
		double ligTorque[] = new double[3];
		double ligTrans[] = new double[3];

		sysDihedDiff = new double[numSysDihedrals];
		ligDihedDiff = new double[numLigDihedrals];
		sysCumulativeDihedStep = new double[numSysDihedrals];
		ligCumulativeDihedStep = new double[numLigDihedrals];

		if(ligStrNum != -1){
			// get the staring COM
			ligStartCOM = m.getStrandCOM(ligStrNum);
			ligCurCOM[0] = ligStartCOM[0];
			ligCurCOM[1] = ligStartCOM[1];
			ligCurCOM[2] = ligStartCOM[2];
		}
				
		// Initialize the dihedral movement arrays
		for(int j=0;j<numSysDihedrals;j++){
			sysDihedDiff[j] = 0.0;
			sysCumulativeDihedStep[j] = 0.0;
		}
		for(int j=0;j<numLigDihedrals;j++){
			ligDihedDiff[j] = 0.0;
			ligCumulativeDihedStep[j] = 0.0;
		}
		
		// If computing dihedral energies initialize them
		if (doDihedEnergy){
			if (!setupDihedralTerms()) //could not initialize dihed energies
				System.exit(1);
		}
		
		float deltaStep = step / numSteps;
		float deltaLigRotStep = lligRotStep / numSteps;
		float deltaLigTransStep = lligTransStep / numSteps;		
		
		int ligResNumber;
		if (!ligandOnly){
			// numFlexRes, flexResAtomList, and flexResListSize include the ligand if one exists
			if(ligStrNum != -1)
				a96ff.setupPartialArrays(numFlexRes+2,MAX_NUM_ATOMS_DISTAL,flexResAtomList,
					flexResListSize);
			else
				a96ff.setupPartialArrays(numFlexRes,MAX_NUM_ATOMS_DISTAL,flexResAtomList,
					flexResListSize);
			
			ligResNumber = ligResNum;
		}
		else //flag Amber not to use partial arrays
			ligResNumber = -1;

		while(!done){
			if (!ligandOnly){
				for(int j=0;j<numSysDihedrals;j++) {
					sysDihedDiff[j] = computeDihedDiff(sysDihedralAtNums[j],sysDihedralDistal[j],
						sysNumAtomsDistal[j],sysDihedToResNum[j], step, false, j);
					updateCumulative(sysCumulativeDihedStep,sysDihedDiff,j,lmaxMovement);
					applyDihedStep(sysDihedralAtNums[j],sysDihedralDistal[j],sysNumAtomsDistal[j],sysDihedDiff[j]);
				}
			}
				
			for(int j=0;j<numLigDihedrals;j++) {
				ligDihedDiff[j] = computeDihedDiff(ligDihedralAtNums[j],ligDihedralDistal[j],
					ligNumAtomsDistal[j],ligResNumber,step,true,j);
				updateCumulative(ligCumulativeDihedStep,ligDihedDiff,j,lmaxMovement);
				applyDihedStep(ligDihedralAtNums[j],ligDihedralDistal[j],ligNumAtomsDistal[j],ligDihedDiff[j]);
			}

			//Translate and rotate the ligand
			if(ligStrNum != -1)				
				doLigTransRot(ligTorque, ligTrans, lligRotStep, lligTransStep, lligMaxTrans);

			i++;
			if(i>=numSteps){
				done = true;
			}

			step -= deltaStep;
			lligRotStep -= deltaLigRotStep;
			lligTransStep -= deltaLigTransStep;
		}
		
		clearMolGradient(); //after minimization is done, clear the molecule gradient
		
		if(debug){
			// Display movement
			if (!ligandOnly){
				for(int j=0;j<numSysDihedrals;j++){
					System.out.print(sysCumulativeDihedStep[j] + " ");
				}
			}
			for(int j=0;j<numLigDihedrals;j++){
				System.out.print(ligCumulativeDihedStep[j] + " ");
			}
			System.out.println();
		}
	}
////////////////////////////////////////////////////////////////////////////////
//	 End of Minmization Section
////////////////////////////////////////////////////////////////////////////////
	
	//Backup the actual ligand coordinates (assumes there is only 1 residue in the ligand strand)
	private float [] backupLigCoord(){		
		float bckpLigCoords[] = new float[m.strand[ligStrNum].numberOfAtoms*3];
		for (int i=0; i<m.strand[ligStrNum].numberOfAtoms; i++){
			int curLigAtom = m.strand[ligStrNum].residue[0].atom[i].moleculeAtomNumber;
			bckpLigCoords[i*3] = m.actualCoordinates[curLigAtom*3];
			bckpLigCoords[i*3+1] = m.actualCoordinates[curLigAtom*3+1];
			bckpLigCoords[i*3+2] = m.actualCoordinates[curLigAtom*3+2];
		}
		return bckpLigCoords;
	}
	
	//Restore the backup actual ligand coordinates (assumes there is only 1 residue in the ligand strand)
	private void restoreLigCoord(float bckpLigCoords[]){			
		for (int i=0; i<m.strand[ligStrNum].numberOfAtoms; i++){
			int curLigAtom = m.strand[ligStrNum].residue[0].atom[i].moleculeAtomNumber;
			m.actualCoordinates[curLigAtom*3] = bckpLigCoords[i*3];
			m.actualCoordinates[curLigAtom*3+1] = bckpLigCoords[i*3+1];
			m.actualCoordinates[curLigAtom*3+2] = bckpLigCoords[i*3+2];
		}
	}
}
