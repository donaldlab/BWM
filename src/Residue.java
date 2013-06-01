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

////////////////////////////////////////////////////////////////////////////////////////////////
// Residue.java
//
//  Version:           1.0
//
//
// authors:
//    initials    name                 organization                email
//   ---------   -----------------    ------------------------    ----------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//
///////////////////////////////////////////////////////////////////////////////////////////////

/*
 *
 * Major changes were made by Ryan Lilien (2001-2004);
 * 		other changes made by Ivelin Georgiev (2004-2009)
 *
 * Rewritten by Ryan Lilien based on code by Neill White
 * Many functions have been added, others removed, most have had 
 *  at least some parts rewritten. Code rewrites have often been
 *  major to fix bugs or add functionality.
 * 
 * Based on software copyrighted, 1999, by Neill White. 
 *  The author hereby grants permission to use, copy, modify, and re-distribute
 *  this software and its documentation for any purpose, provided
 *  that existing copyright notices are retained in all copies and that this
 *  notice is included verbatim in any distributions. No written agreement,
 *  license, or royalty fee is required for any of the authorized uses.
 *  Modifications to this software may be distributed provided that
 *  the nature of the modifications are clearly indicated.
 *
 */

import java.util.StringTokenizer;
import java.io.Serializable;

/**
 * Handles functions and data associated with residues.
 */
public class Residue implements Serializable {

	String name = "";			// The short residue name "GLU"
	String fullName = "";	// The long version of the residue name ie "GLU A  47"
	int	numberOfAtoms=0;	// Number of atoms
	int	moleculeResidueNumber=-1;	// The molecule relative residue number
	int	strandResidueNumber=-1;		// The strand relative residue number
	int	strandNumber=-1;	// The number of the strand containing this residue
	Atom	atom[];					// Array of atoms in this residue
	private boolean energyEvalSC = true;	// Should the residue side-chain be included in computing energy
	private boolean energyEvalBB = true;	// Should the residue backbone be included in computing energy
	boolean ffAssigned = false;	// Are forcefield types assigned for atoms in this residue
	boolean flexible = false;	// Is this residue flexible
	
		
	Residue(){
		atom = new Atom[1];
	}

	Residue(String resname){
		atom = new Atom[1];
		name = resname;
	}
	
	// Displays some residue info to System.out for debugging
	public void printResidueInfo(){
		System.out.println("name = *"+name+"*");
		System.out.println("fullName = *"+fullName+"*");
		System.out.println("numberOfAtoms = *"+numberOfAtoms+"*");
		System.out.println("moleculeResidueNumber = *"+moleculeResidueNumber+"*");
		System.out.println("strandResidueNumber = *"+strandResidueNumber+"*");
		System.out.println("strandNumber = *"+strandNumber+"*");
	}
	
	// Adds specified atom to this residue
	// Appropriately updates all Residue and Atom fields except for
	//  Atom.moleculeAtomNumber which should be updated in the molecule class
	// *This function should rarely be called. Atoms should be added via the
	//  molecule containing this residue. This function is important however
	//  as it is called by the molecule and strand addAtom functions
	public int addAtom(Atom newAtom){
		int newAtomNumber = numberOfAtoms + 1;
		int newAtomNumberx3 = newAtomNumber * 3;

    Atom largerAtomArray[] = new Atom[newAtomNumber];
    System.arraycopy(atom, 0, largerAtomArray, 0, atom.length);
    atom = largerAtomArray;
    atom[numberOfAtoms] = newAtom;
    
		newAtom.residueAtomNumber = numberOfAtoms;
		newAtom.moleculeResidueNumber = moleculeResidueNumber;
		newAtom.strandResidueNumber = strandResidueNumber;
		newAtom.strandNumber = strandNumber;

		return(numberOfAtoms++);
	}

	// Deletes the specified atom
	// *This function should rarely be called. Atoms should be deleted via the
	//  molecule containing this residue. This function is important however
	//  as it is called by the molecule and strand deleteAtom functions
	public int deleteAtom(int atomNumber){
		int newAtomNumber = numberOfAtoms - 1;
		int newAtomNumberx3 = newAtomNumber * 3;
		int atomNumberx3 = atomNumber * 3;
		int moleculeAtomNumber = atom[atomNumber].moleculeAtomNumber;
		
		// If there is no molecule atom number then this is a lone
		//  residue (ie. it's not in a molecule, so use the simple
		//  numbering as the moleculeAtomNumber
		if (moleculeAtomNumber == -1)
			moleculeAtomNumber = atomNumber;
		// update residueAtomNumbers of higher numbered atoms
		for(int i=atomNumber; i<numberOfAtoms; i++)
			atom[i].residueAtomNumber -= 1;

		Atom smallerAtomArray[] = new Atom[numberOfAtoms-1];
    System.arraycopy(atom, 0, smallerAtomArray, 0, atomNumber);
		if (atomNumber<newAtomNumber)
   	  System.arraycopy(atom, atomNumber+1, smallerAtomArray,
                    	atomNumber, atom.length-atomNumber-1);
		atom = smallerAtomArray;

		numberOfAtoms--;

		// update bond indices
		for(int i=0;i<numberOfAtoms;i++) {
			for(int m=0; m<atom[i].numberOfBonds; m++){
				if (atom[i].bond[m] == moleculeAtomNumber){
					atom[i].deleteBond(m);
					m--;
				}
				else if (atom[i].bond[m] > moleculeAtomNumber){
					atom[i].bond[m] -= 1;
				}
			}
		}

		return(numberOfAtoms);
	}

	// This function takes a residue and renumbers its
	//  internal indices such that the first atom is atom 0;
	// It adjusts atom members moleculeAtomNumber and residueAtomNumber;
	// It resets atom bond[] and numberOfBonds, so the bond information must be updated 
	//  using the determineBonds() method from Molecule.java;
	// This function is used when we pull a residue out
	//  from one molecule and we want to add it back
	//  into another one.
	public void renumberResidue(){
		for(int i=0;i<numberOfAtoms;i++){
			atom[i].moleculeAtomNumber = i;
			atom[i].residueAtomNumber = i;
			atom[i].bond = null;	
			atom[i].numberOfBonds = 0;
		}
	}


	// This function returns the xth token in string s
	private String getToken(String s, int x) {
	
		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");
		
		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
			  st.nextToken();
			else {
				System.out.println("ERROR: Unable to access parameter " + x);
				return(new String(""));
			}
		}
		
		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));
	} // end getToken

	// This function returns the residue number (not the sequential one,
	//  but the 'true' one from the pdb file).  If the pdb doesn't have
	//  one (unlikely) then use the sequential numbering
	public int getResNumber() {
		if (fullName.length() > 5)
			return( (new Integer(getToken(fullName.substring(5),1)).intValue()) );
		return (moleculeResidueNumber+1);
	}

	// This function rotates the specified atoms in the residue
	//  by thetaDeg degrees around axis dx, dy, dz
	// at1 is the pivot atom (the 3rd atom if we were doing dihedrals)
	//  where the 4th atom+ would be the moving atoms
	// atomList is a list of atom numbers to rotate
	// numAtoms is the number of atoms to rotate (size of atomList)
	public void rotateResidue(Atom at1, double dx, double dy,
		double dz, double thetaDeg, int atomList[], int numAtoms) {
	
		float fx,fy,fz, tx,ty,tz;
		fx = (new Double(dx)).floatValue();
		fy = (new Double(dy)).floatValue();
		fz = (new Double(dz)).floatValue();
		Atom pseudoAtom = new Atom("a", fx, fy, fz);
	
		int atomNumber;
		if (numAtoms == 0)
			return;
	
		int numberOfCoordinatesx3 = numAtoms*3;
		float temporaryCoordinates[] = new float[numberOfCoordinatesx3];
		int qx3;
		for(int q=0;q<numAtoms;q++) {
			qx3 = q*3;
			temporaryCoordinates[qx3+0]=atom[atomList[q]].coord[0] - at1.coord[0];
			temporaryCoordinates[qx3+1]=atom[atomList[q]].coord[1] - at1.coord[1];
			temporaryCoordinates[qx3+2]=atom[atomList[q]].coord[2] - at1.coord[2];
		}

		float[][] rot_mtx = new float[3][3];
		RotMatrix rM = new RotMatrix();
		rM.getRotMatrix(fx,fy,fz,(float) thetaDeg,rot_mtx);

		for(int q=0;q<numAtoms;q++) {
			qx3 = q*3;
			tx = temporaryCoordinates[qx3+0];
			ty = temporaryCoordinates[qx3+1];
			tz = temporaryCoordinates[qx3+2];

			atom[atomList[q]].coord[0] = tx * rot_mtx[0][0] + ty * rot_mtx[0][1] + tz * rot_mtx[0][2] + at1.coord[0];
			atom[atomList[q]].coord[1] = tx * rot_mtx[1][0] + ty * rot_mtx[1][1] + tz * rot_mtx[1][2] + at1.coord[1];
			atom[atomList[q]].coord[2] = tx * rot_mtx[2][0] + ty * rot_mtx[2][1] + tz * rot_mtx[2][2] + at1.coord[2];
		}

	}
	
	//Returns the distance (the minimum distance between a pair of non-hydrogen side-chain atoms) between the side-chains of this residue and res2;
	//If bbAt is false, then only side-chain atoms are considered
	public float getDist(Residue res2, boolean bbAt){
		
		float minDist = (float)Math.pow(10, 10);
		
		for (int a1=0; a1<numberOfAtoms; a1++){
			
			String a1type = atom[a1].elementType;
			
			if ( (!a1type.equalsIgnoreCase("H")) && ( bbAt || (!atom[a1].getIsBBatom()) ) ){
				
				for (int a2=0; a2<res2.numberOfAtoms; a2++){
					
					String a2type = res2.atom[a2].elementType;
					
					if ( (!a2type.equalsIgnoreCase("H")) && (bbAt || (!res2.atom[a2].getIsBBatom()) ) ){
						
						float curDist = getDist(atom[a1],res2.atom[a2]);
						minDist = Math.min(minDist, curDist);
					}
				}
			}
		}
		
		return minDist;
	}
	
	//Returns the distance between the two atoms
	private float getDist(Atom a1, Atom a2){
			
		float rijx, rijy, rijz, rij, rij2;
		
		rijx = a1.coord[0] - a2.coord[0];
		rijy = a1.coord[1] - a2.coord[1];
		rijz = a1.coord[2] - a2.coord[2];
		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		rij = (float)Math.sqrt(rij2);
		
		return rij;
	}
	
	public void setEnergyEval(boolean scEval, boolean bbEval){
		energyEvalSC = scEval;
		energyEvalBB = bbEval;
	}
	
	public boolean getEnergyEvalSC(){
		return energyEvalSC;
	}
	
	public boolean getEnergyEvalBB(){
		return energyEvalBB;
	}
	
	//Returns the moleculeAtomNumber for the atom of this residue that has the same name as the parameter n
	public int getAtomNameToMolnum(String n){
		
		for (int i=0; i<numberOfAtoms; i++){
			if (atom[i].name.equalsIgnoreCase(n))
				return atom[i].moleculeAtomNumber;
		}
		return -1;
	}
}
