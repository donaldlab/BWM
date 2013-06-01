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

////////////////////////////////////////////////////////////////////////////////////////////
// Molecule.java
//
//  Version:           1.0
//
//
// authors:
//    initials    name            organization                email
//   ---------   --------------  ------------------------    ------------------------------
//     RHL        Ryan Lilien     Dartmouth College           ryan.lilien@dartmouth.edu
//
////////////////////////////////////////////////////////////////////////////////////////////

/*
 *
 * Major changes were made by Ryan Lilien (2001-2004) and Ivelin Georgiev (2004-2009)
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

import java.util.*;
import java.io.Serializable;

/**
 * Handles functions and data associated with molecules. Handles rotations/translations of (parts of) molecules.
 * Manages the data associated with a molecule; handles changes to the molecule (e.g., coordinate changes,
 * deletion or addition of residues, etc.). Determines the bond information for the molecule.
 *
 */
public class Molecule implements Serializable{

	// All adding and subtracting of strands, residues, and atoms are "top down".
	// So, adding a strand, residue, or atom, should be done to a MOLECULE,
	//  which, in turn will propagate down and update the other classes. E.g.,
	//  if I want to add an atom, I call m.addAtom( stuff ), which will then 
	//  update the local arrays, then call strand.addAtom() which will update
	//  the strand info, then call residue.addAtom() which will update the
	//  residue info.
	 
	// Proper functioning of all ILMM classes relies upon a proper increasing
	//  numerical ordering of the atoms, residues, and strands. This means for
	//  all i < j all atoms in residue[i] are numbered less than all atoms in
	//  residue[j]. Additionally all residues in strand[i] are numbered less
	//  than all residues in strand[j]

	// If the debug flag is set to true then additional debug statements are
	//  printed to standard out. I'm hoping that by making it a public static
	//  final variable that the compiler will be smart and not compile/include
	//  code that is unreachable.
	public static final boolean debug = false;
	// Connect_Increment_Size is the chunk by which the connected arrays grow
	//  when they run out of space
	public static final int CONNECT_INCREMENT_SIZE = 6000;
	public static final int NB_INCREMENT_SIZE = 82000;
	public static final int MAX_ATOM_CONNECTIONS = 7;
	
	String	name;
	int	numberOfStrands = 0;
	int	numberOfResidues = 0;
	int	numberOfAtoms = 0;
	int	numberOfAtomsx3 = 0;
	Strand strand[] = null;			// Array of strands in this molecule
	Residue	residue[] = null;		// Array of residues in this molecule
	Atom atom[] = null;					// Array of atoms in this molecule
	Atom origAtom[] = null;				// Array of atoms to store a copy of the molecule atoms
	float backupCoordinates[] = null;	// Backup of atom coordinates
	float	actualCoordinates[] = null;		// Local copy of atom coordinates
	boolean connectivityValid = false;	// Are the connectivity arrays up to date
	boolean connectivity12Valid = false;	// Are connected and connected12 up to date
	int	connected[][] = null;		// 2D array of connected atoms
	int	connected12[] = null;		// 1D array of atoms directly connected
	int	connected13[] = null;		// 1D array of atoms 1-3 connected (1 atom between)
	int	connected14[] = null;		// 1D array of atoms 1-4 connected (2 atoms between)
	int	nonBonded[] = null;			// List of non bonded pairs
	boolean bondedMatrix[][] = null;
		// entries are only in half of bondedMatrix
		// so bondedMatrix{i][j] is only valid where i<j
	int bondedMatrixSize = 0;
		// bondedMatrix size is bondedMatrixSize x bondedMatrixSize
	int	numberOf12Connections = 0;
	int	numberOf13Connections = 0;
	int	numberOf14Connections = 0;
	int	numberNonBonded = 0;
	double	gradient[] = null;	// The computed gradient for this molecule (x,y,z) for each atom

	// Generic constructor
	Molecule(){
		strand = new Strand[1];
		residue = new Residue[1];
		residue[0] = new Residue();
		strand[0] = new Strand(residue[0]);
		atom = new Atom[1];
		actualCoordinates = new float[1];
		name = "Untitled";
	}

	// Prints molecule information
	public String toString(){
		System.out.println("name = "+name );
		System.out.println("numberOfStrands = "+numberOfStrands );
		for ( int i = 0; i < numberOfStrands; i++ )
			System.out.println( "\t" + i + ": " + strand[ i ].name );	
		System.out.println("numberOfResidues = "+numberOfResidues );
		for ( int i = 0; i < numberOfResidues; i++ )
			System.out.println( "\t" + i + ": " + residue[ i ].name );	
		System.out.println("numberOfAtoms = "+numberOfAtoms );
		for( int i = 0; i < numberOfAtoms; i++ ){
			System.out.print( "("+i+")"+atom[i].name + " : " +
				actualCoordinates[ i*3 ] + ", " + 
				actualCoordinates[ i*3+1 ] + ", " + 
				actualCoordinates[ i*3+2 ] );
			System.out.print( "; bonded to: " );
			for ( int j = 0; j < atom[ i ].numberOfBonds; j++ )
				System.out.print( atom[ i ].bond[ j ] + ", " );
			System.out.println();
		}
		return(new String(""));
	}

	// This function checks a molecule to make sure that it's legit
	// There are various types of checking
	// 1: Check Numbering
	//     Make sure Residues and Strands are numerically increasing
	//      with no gaps
	//     Make sure Atoms are numerically increasing with no gaps
	//     Make sure the following fields in each atom are correct:
	//      moleculeAtomNumber, residueAtomNumber,
	//      moleculeResidueNumber, strandResidueNumber,
	//      strandNumber
	//     Make sure the following fields in each strand and residue
	//      are correct: Strand: numberOfResidues, numberOfAtoms
	//      Residue: numberOfAtoms, moleculeResidueNumber,
	//      strandResidueNumber, strandNumber
	// 2: Check Bonds
	//     Make sure all bonds incident on an atom in strand i stay
	//      in strand i
	//     Make sure there are no bonds to atoms that don't exist
	//     Make sure each atom's numberOfConnections is accurate
	//     Make sure the molecule's bond member is accurate
	// vlevel is the verbosity level
	//  0: print nothing
	//  1: print only errors
	//  2: print everything
	//  
	public int checkMolecule(int level, int vlevel){
		boolean correct = true;
		boolean correct2 = true;
		boolean correct3 = true;
		
		if(level==1){
			int lastSnum, lastRnum, lastAnum;
			if(vlevel>=2){
				System.out.println();
				System.out.println("Checking molecule based numbering of Strands, Residues, and Atoms...");
			}
			lastSnum = lastRnum = lastAnum = -1;
			for(int i=0; i<numberOfStrands; i++){
				if (strand[i].number != ++lastSnum){
					if(vlevel>=1)
						System.out.println("  Strand " + i + " is misnumbered");
					correct=false;
				}
				for(int j=0; j<strand[i].numberOfResidues; j++){
					if(strand[i].residue[j].moleculeResidueNumber != ++lastRnum){
						if(vlevel>=1)
							System.out.println("  Residue " + j + " in strand " + i + " is misnumbered");
						correct=false;
					}
					for(int k=0; k<strand[i].residue[j].numberOfAtoms; k++){
						if(strand[i].residue[j].atom[k].moleculeAtomNumber != ++lastAnum){
							if(vlevel>=1)
								System.out.println("  Atom " + k + " in residue " + j + " in strand " + i + " is misnumbered");
							correct=false;
						}
					}
				}
			}
			if (correct){
				if(vlevel>=2)
					System.out.println("  Strand, Residue, and Atom -- molecule relative numbers are correct");
			}
			else {
				if(vlevel>=1)
					System.out.println("  *** Strand, Residue, and Atom -- molecule relative numbers have errors");
			}
			if(vlevel>=2)
				System.out.println("Checking strand based numbering of Residues...");
			
			correct = true;
			for(int i=0; i<numberOfStrands; i++){
				lastRnum = -1;
				for(int j=0; j<strand[i].numberOfResidues; j++){
					if(strand[i].residue[j].strandResidueNumber != ++lastRnum){
						if(vlevel>=1)
							System.out.println("  Residue " + j + " in strand " + i + " is misnumbered");
						correct=false;
					}
				}
			}
			if (correct){
				if(vlevel>=2)
					System.out.println("  Residue -- strand relative numbers are correct");
			}
			else {
				if(vlevel>=1)
					System.out.println("  *** Residue -- strand relative numbers have errors");
			}
			
			if(vlevel>=2)
				System.out.println("Checking residue based numbering of Atoms...");
				
			correct = true;
			for(int i=0; i<numberOfStrands; i++){
				for(int j=0; j<strand[i].numberOfResidues; j++){
					lastAnum=-1;
					for(int k=0; k<strand[i].residue[j].numberOfAtoms; k++){
						if(strand[i].residue[j].atom[k].residueAtomNumber != ++lastAnum){
							if(vlevel>=1)
								System.out.println("  Atom " + k + " in residue " + j + " in strand " + i + " is misnumbered");
							correct=false;
						}
					}
				}
			}
			if (correct){
				if(vlevel>=2)
					System.out.println("  Atom -- residue relative numbers are correct");
			}
			else {
				if(vlevel>=1)
					System.out.println("  *** Atom -- residue relative numbers have errors");
			}
			
			if(vlevel>=2)
				System.out.println("Checking atom.strandNumber...");
				
			correct = true;
			for(int i=0; i<numberOfStrands; i++){
				for(int j=0; j<strand[i].numberOfResidues; j++){
					for(int k=0; k<strand[i].residue[j].numberOfAtoms; k++){
						if(strand[i].residue[j].atom[k].strandNumber != i){
							if(vlevel>=1)
								System.out.println("  Atom " + k + " in residue " + j + " of strand " + i + " has misnumbered .strandNumber");
							correct=false;
						}
					}
				}
			}
			if (correct){
				if(vlevel>=2)
					System.out.println("  atom.strandNumbers are correct");
			}
			else {
				if(vlevel>=1)
					System.out.println("  *** atom.strandNumbers have errors");
			}

			if(vlevel>=2)
				System.out.println("Checking atom.strandResidueNumber...");
				
			correct = true;
			for(int i=0; i<numberOfStrands; i++){
				for(int j=0; j<strand[i].numberOfResidues; j++){
					for(int k=0; k<strand[i].residue[j].numberOfAtoms; k++){
						if(strand[i].residue[j].atom[k].strandResidueNumber != j){
							if(vlevel>=1)
								System.out.println("  Atom " + k + " in residue " + j + " of strand " + i + " has misnumbered .strandResidueNumber");
							correct=false;
						}
					}
				}
			}
			if (correct){
				if(vlevel>=2)
					System.out.println("  atom.strandResidueNumbers are correct");
			}
			else {
				if(vlevel>=1)
					System.out.println("  *** atom.strandResidueNumbers have errors");
			}

			if(vlevel>=2)
				System.out.println("Checking atom.moleculeResidueNumber...");
				
			correct = true;
			for(int i=0; i<numberOfResidues; i++){
				for(int j=0; j<residue[i].numberOfAtoms; j++){
					if(residue[i].atom[j].moleculeResidueNumber != i){
						if(vlevel>=1)
							System.out.println("  Atom " + j + " in residue " + i + " of molecule has misnumbered .moleculeResidueNumber");
						correct=false;
					}
				}
			}
			if (correct){
				if(vlevel>=2)
					System.out.println("  atom.moleculeResidueNumbers are correct");
			}
			else {
				if(vlevel>=1)
					System.out.println("  *** atom.moleculeResidueNumbers have errors");
			}

			
			if(vlevel>=2)
				System.out.println("Checking residue.strandNumber...");
				
			correct = true;
			for(int i=0; i<numberOfStrands; i++){
				for(int j=0; j<strand[i].numberOfResidues; j++){
					if(strand[i].residue[j].strandNumber != i){
						if(vlevel>=1)
							System.out.println("  Residue " + j + " of strand " + i + " has misnumbered .strandNumber");
						correct=false;
					}
				}
			}
			if (correct){
				if(vlevel>=2)
					System.out.println("  residue.strandNumbers are correct");
			}
			else {
				if(vlevel>=1)
					System.out.println("  *** residue.strandNumbers have errors");
			}

			if(vlevel>=2)
				System.out.println("Checking residue.numberOfAtoms...");
				
			correct = true;
			for(int i=0; i<numberOfStrands; i++){
				for(int j=0; j<strand[i].numberOfResidues; j++){
					if(strand[i].residue[j].numberOfAtoms > 1){
						if(strand[i].residue[j].atom.length != strand[i].residue[j].numberOfAtoms){
							if(vlevel>=1)
								System.out.println("  Residue " + j + " of strand " + i + " has mismatched residue.atom.length residue.numberOfAtoms");
							correct=false;
						}
					}
				}
			}
			if (correct){
				if(vlevel>=2)
					System.out.println("  residue.numberOfAtoms are correct");
			}
			else {
				if(vlevel>=1)
					System.out.println("  *** residue.numberOfAtoms have errors");
			}

			if(vlevel>=2)
				System.out.println("Checking strand.numberOfResidues...");
			correct = true;
			for(int i=0; i<numberOfStrands; i++){
				if(strand[i].numberOfResidues > 1){
					if(strand[i].residue.length != strand[i].numberOfResidues){
						if(vlevel>=1)
							System.out.println("  Strand " + i + " has mismatched strand.residue.length strand.numberOfResidues");
						correct=false;
					}
				}
			}
			if (correct){
				if(vlevel>=2)
					System.out.println("  strand.numberOfResiues are correct");
			}
			else {
				if(vlevel>=1)
					System.out.println("  *** strand.numberOfResiues have errors");
			}

			if(vlevel>=2)
				System.out.println("Checking strand.numberOfAtoms...");
			correct = true;
			for(int i=0; i<numberOfStrands; i++){
				int atCnt=0;
				for(int j=0; j<strand[i].numberOfResidues; j++)
					atCnt += strand[i].residue[j].numberOfAtoms;
				if (atCnt != strand[i].numberOfAtoms){
					if(vlevel>=1)
						System.out.println("  Strand " + i + " has a strand.numberOfAtoms that doesn't match with its subresidues");
					correct=false;
				}
			}
			if (correct){
				if(vlevel>=2)
					System.out.println("  strand.numberOfAtoms are correct");
			}
			else {
				if(vlevel>=1)
					System.out.println("  *** strand.numberOfAtoms have errors");
			}

			if(vlevel>=2)
				System.out.println();		
			return(0);
		}
		else if(level==2){
			correct = true;
			correct2 = true;
			correct3 = true;
			if(vlevel>=2){
				System.out.println();
				System.out.println("Checking bonds...");
			}
			for(int i=0;i<numberOfAtoms;i++){
				if(atom[i].bond == null){
					if(atom[i].numberOfBonds != 0){
						if(vlevel>=1)
							System.out.println("  atomMoleculeNumber " + i + " has a bond array with size different than numberOfBonds");
						correct2=false;
					}
				}
				else {
					if (atom[i].bond.length != atom[i].numberOfBonds){
						if(vlevel>=1)
							System.out.println("  atomMoleculeNumber " + i + " has a bond array with size different than numberOfBonds");
						correct2=false;
					}
				}
				
				for(int j=0;j<atom[i].numberOfBonds;j++){
					if (atom[atom[i].bond[j]].strandNumber != atom[i].strandNumber){
						if(vlevel>=1)
							System.out.println("  atomMoleculeNumber " + i + " has a bond outside its strand");
						correct3=false;
					}
				}
			}
			
			if (correct){
				if(vlevel>=2)
					System.out.println("  Atom numberOfConnections are correct");
			}
			else {
				if(vlevel>=1)
					System.out.println("  *** Atom numberOfConnections have errors");
			}
			if (correct2){
				if(vlevel>=2)
					System.out.println("  Atom bond array sizes agree with number of bonds");
			}
			else {
				if(vlevel>=1)
					System.out.println("  *** Atom bond array size doesn't agree with number of bonds (single,double,triple bonds)");
			}
			if (correct3){
				if(vlevel>=2)
					System.out.println("  Atom bonds are all correctly intrastrand");
			}
			else {
				if(vlevel>=1)
					System.out.println("  *** Atom bonds are not intrastrand");
			}
			
			if(vlevel>=2)
				System.out.println();
			return(0);		
		}
		
		
		return(0);
	}

	// This function goes through a molecule and searches for a few situations
	//  that could happen in normal functioning but that shouldn't exist in the
	//  real world. If the input boolean fix is true this function attempts to
	//  fix the situation.
	// checkLevel
	//  1: check for bonds between residue sidechains in protein strands
	// The return value is 0 if there were no problems, 1 if there was a problem
	//  we couldn't fix, and 2 if there was a problem and we fixed it
	public int sanityCheck(int checkLevel, boolean fix){
	
		Atom at1;
		int result = 0;
		int at2res = -1;
		boolean cantfix = false;
		
		if(checkLevel==1){
			for(int i=0; i<numberOfStrands; i++){
				if(strand[i].isProtein){
					for(int j=0; j<strand[i].numberOfResidues; j++){
						for(int k=0; k<strand[i].residue[j].numberOfAtoms; k++){
							at1 = strand[i].residue[j].atom[k];
							for(int l=0; l<at1.numberOfBonds; l++){
								at2res = atom[at1.bond[l]].moleculeResidueNumber;
								if(at2res != at1.moleculeResidueNumber){
									if(!(at1.name.equals("C") || at1.name.equals("N"))){
										result = 1;
										int at2num = atom[at1.bond[l]].moleculeAtomNumber;
										if(debug)
											System.out.print("Error found, bond between " + at1.moleculeAtomNumber + " and " + at2num);
										if (fix) {
											deleteBondBetween(at1.moleculeAtomNumber, at2num);
											result = 2;
											if(debug)
												System.out.print("...fixed");
										}
										if(debug)
											System.out.println();
									}
								}
							}
						}
					}
				}
			}	
		
		}
	
		if(cantfix)
			result=1;
		return(result);
	}

	// If numResidues == 0, numAtoms > 0 then
	//   space is made for atoms in the atom[] and actualCoordinates[] after the
	//   last atom of the residue at strand strandPos residue residuePos
	// If numResidues == 1, numAtoms > 0 then
	//   space is made for atoms and 1 new residue in atom[], actualCoordinates[],
	//   and residue[] after the end of strand strandPos (residuePos is not used)
	// The molecule numbers of all atoms and residues above strandPos, residuePos
	//   are updated.
	// The molecule number of atoms and residues have NOT been updated.
	// The structures residue.atom[] and strand.residue[] are not changed.
	// strandPos is molecule based (ie. moleculeResidueNumber)
	private void expandMolecule(int numResidues, int numAtoms, int strandPos, int residuePos){
	
		// check for valid numResidues
		if (numResidues > 1)
			return;
			
		// first do things common to both cases
		Strand theStrand = strand[strandPos];
		Residue theResidue;
		int moleculeAtomInsertPoint = 0;
		int moleculeResInsertPoint = 0;

		if (numResidues == 1){
			theResidue = new Residue();
			if (numberOfResidues > 0){
				for(int i=0;i<strandPos+1;i++)
					moleculeResInsertPoint += strand[i].numberOfResidues;
				if (moleculeResInsertPoint == 0)
					moleculeAtomInsertPoint = 0;
				else
					moleculeAtomInsertPoint = residue[moleculeResInsertPoint-1].atom[residue[moleculeResInsertPoint-1].numberOfAtoms-1].moleculeAtomNumber+1;
			}
		}
		else{
			theResidue = residue[residuePos];
			moleculeAtomInsertPoint = theResidue.atom[0].moleculeAtomNumber + theResidue.numberOfAtoms;
		}

		for(int i=0;i<numberOfAtoms;i++){
			// update atom.moleculeAtomNumber
			if (atom[i].moleculeAtomNumber >= moleculeAtomInsertPoint)
				atom[i].moleculeAtomNumber += numAtoms;
			// if (numResidues==1) update atom.moleculeResidueNumber
			if (numResidues==1)
				if (atom[i].moleculeResidueNumber >= moleculeResInsertPoint)
					atom[i].moleculeResidueNumber += 1;
			// update bond references
			for(int k=0; k<atom[i].numberOfBonds; k++)
				if (atom[i].bond[k] >= moleculeAtomInsertPoint)
					atom[i].bond[k] += numAtoms;
		}

		if (numResidues==1){
			for(int i=0;i<numberOfResidues;i++){
				if (residue[i].moleculeResidueNumber >= moleculeResInsertPoint)
					residue[i].moleculeResidueNumber += 1;
			}		
		}
		
		// add the atom and residue space in atom[] and residue[]		
		int newAtomNumber = numberOfAtoms + numAtoms;
		Atom largerAtomArray[] = new Atom[newAtomNumber];
		System.arraycopy(atom, 0, largerAtomArray, 0, moleculeAtomInsertPoint);
		if (moleculeAtomInsertPoint+numAtoms < newAtomNumber)
			System.arraycopy(atom, moleculeAtomInsertPoint, largerAtomArray,
				numAtoms+moleculeAtomInsertPoint, newAtomNumber-(moleculeAtomInsertPoint+numAtoms));
		atom = largerAtomArray;
		
		if (numResidues==1) {
			Residue largerResidueArray[] = new Residue[numberOfResidues + 1];
			System.arraycopy(residue, 0, largerResidueArray, 0, moleculeResInsertPoint);
			if (moleculeResInsertPoint < numberOfResidues)
				System.arraycopy(residue, moleculeResInsertPoint, largerResidueArray,
					moleculeResInsertPoint+1, numberOfResidues-moleculeResInsertPoint);
			residue = largerResidueArray;
		}			
		// handle coordinates
		int newAtomNumberx3 = newAtomNumber * 3;
		int moleculeAtomInsertPointx3 = moleculeAtomInsertPoint * 3;
		float largerCoordinateArray[] = new float[newAtomNumberx3];	
		System.arraycopy(actualCoordinates, 0, largerCoordinateArray, 0,
			moleculeAtomInsertPointx3);
		if (moleculeAtomInsertPoint+numAtoms < newAtomNumber)
			System.arraycopy(actualCoordinates, moleculeAtomInsertPointx3,
				largerCoordinateArray, numAtoms*3+moleculeAtomInsertPointx3,
				newAtomNumberx3-(moleculeAtomInsertPointx3+numAtoms*3));
		actualCoordinates = largerCoordinateArray;
	}

	// Adds strand (using below fn) - if no name, call it the strandnum
	public int addStrand(){
		return(addStrand(new String().valueOf(numberOfStrands)));
	}

	// Adds a new empty strand with name "strandName" to molecule and updates array. 
	// The strand number is just the next number in the strand array
	// All molecule and strand fields are appropriately updated
	public int addStrand(String strandName){
		connectivityValid = false;
		connectivity12Valid = false;
		Strand largerStrandArray[] = new Strand[numberOfStrands + 1];
		System.arraycopy(strand,0,largerStrandArray,0,strand.length);
		strand = largerStrandArray;
		strand[numberOfStrands] = new Strand(strandName);
		strand[numberOfStrands].number = numberOfStrands;
		return(numberOfStrands++);
	}

	// Adds an existing strand to this molecule
	// All molecule, strand, residue, and atom fields are appropriately updated
	public int addStrand(Strand stnd) {
	
		connectivityValid = false;
		connectivity12Valid = false;
		
		int newAtomNumber = numberOfAtoms + stnd.numberOfAtoms;
		int newAtomNumberx3 = newAtomNumber * 3;
		int newResidueNumber = numberOfResidues + stnd.numberOfResidues;

		// Add the strand itself
		Strand largerStrandArray[] = new Strand[numberOfStrands + 1];
		System.arraycopy(strand, 0, largerStrandArray, 0, strand.length);
		strand = largerStrandArray;
		strand[numberOfStrands] = stnd;
		strand[numberOfStrands].number = numberOfStrands;

		// Add the residues to the molecule array
		Residue largerResidueArray[] = new Residue[newResidueNumber];
		System.arraycopy(residue, 0, largerResidueArray, 0, residue.length);
		for(int i=numberOfResidues, j=0;i<newResidueNumber;i++,j++) {
			largerResidueArray[i] = stnd.residue[j];
		}
		// Shift new residue numbers up by numberOfResidues, update strandNumber
		for(int i=numberOfResidues, j=0; i<newResidueNumber; i++,j++) {
			largerResidueArray[i].moleculeResidueNumber = numberOfResidues + j;
			largerResidueArray[i].strandNumber = numberOfStrands;
		}
		residue = largerResidueArray;
		
		// Copy over and update Atoms
		Atom largerAtomArray[] = new Atom[newAtomNumber];
		System.arraycopy(atom, 0, largerAtomArray, 0, atom.length);
		int q = numberOfAtoms;
		for(int i=0; i<stnd.numberOfResidues; i++) {
			for(int j=0; j<stnd.residue[i].numberOfAtoms; j++) {
				largerAtomArray[q] = stnd.residue[i].atom[j];
				largerAtomArray[q].moleculeAtomNumber = q;
				largerAtomArray[q].moleculeResidueNumber = numberOfResidues + i;
				largerAtomArray[q].strandNumber = numberOfStrands;
				q++;
			}
		}
		atom = largerAtomArray;

		// Handle coordinates (and update bond numbers)
		float largerCoordinateArray[] = new float[ newAtomNumberx3 ];	
		System.arraycopy(actualCoordinates, 0, largerCoordinateArray, 0,
			actualCoordinates.length);
		actualCoordinates = largerCoordinateArray;
		int ix3;
		for(int i=numberOfAtoms; i<newAtomNumber; i++) {
			ix3 = i * 3;
			actualCoordinates[ix3] = atom[i].coord[0];
			actualCoordinates[ix3+1] = atom[i].coord[1];
			actualCoordinates[ix3+2] = atom[i].coord[2];
			// These 3 bond updates assume each strand starts with new numbering
			for(int k=0; k<atom[i].numberOfBonds; k++)
				atom[i].bond[k] += numberOfAtoms;
		}

		numberOfAtoms = newAtomNumber;
		numberOfAtomsx3 = newAtomNumberx3;
		numberOfResidues = newResidueNumber;
		numberOfStrands++;
		// Update connectivity tables
		establishConnectivity(false);
		return(numberOfStrands);
	}
	
	//Calls the addResidue() version below with (updateBonds==false);
	//Only call this version if:
	//	1) determineBonds() is called immediately after this call, or
	//	2) determineBonds() is called later and the atoms bonds[] array is not used immediately, or
	//	3) the atoms bonds[] array will not be used
	public void addResidue(int strandNumber, Residue newResidue){
		addResidue(strandNumber, newResidue, false);
	}

	// Public method which adds a residue and all atoms to the molecule
	// All molecule, strand, residue, and atom fields are appropriately updated;
	// If updatedBonds is true, then the bonds[] array for each atom in newResidue is computed
	public void addResidue(int strandNumber, Residue newResidue, boolean updateBonds){
		
		if (strandNumber>=numberOfStrands){
			addStrand();
			strandNumber = numberOfStrands-1;
		}
		
		connectivityValid = false;
		connectivity12Valid = false;
		
		int newAtomNumber = numberOfAtoms + newResidue.numberOfAtoms;

		int newResPosition = 0;
		int newAtomPosition = 0;

		if (numberOfResidues > 0){
			for(int i=0;i<strandNumber+1;i++)
				newResPosition += strand[i].numberOfResidues;
			if (newResPosition == 0)
				newAtomPosition = 0;
			else
				newAtomPosition = residue[newResPosition-1].atom[residue[newResPosition-1].numberOfAtoms-1].moleculeAtomNumber+1;
		}

		// make room for the new residue and atoms in atom[],
		//  actualCoordinates[], and residue[]
		// The 4th parameter is the newResiduePosition but it isn't used
		//  by expandMolecule if the first paramter is 1
		expandMolecule(1,newResidue.numberOfAtoms,strandNumber,newResPosition-1);
		// appropriately add residue
		newResidue.moleculeResidueNumber = newResPosition;
		residue[newResPosition]=newResidue;

		// appropriately add atoms		
		int q = newAtomPosition;
		int qx3 = q*3;
		for(int i=0; i<newResidue.numberOfAtoms; i++, q++) {
			atom[q] = newResidue.atom[i];
			atom[q].moleculeAtomNumber = q;
			atom[q].moleculeResidueNumber = newResPosition;
			qx3 = q * 3;
			actualCoordinates[qx3] = atom[q].coord[0] = 
				newResidue.atom[i].coord[0];
			actualCoordinates[qx3+1] = atom[q].coord[1] = 
				newResidue.atom[i].coord[1];
			actualCoordinates[qx3+2] = atom[q].coord[2] = 
				newResidue.atom[i].coord[2];
		}

		strand[strandNumber].addResidue(newResidue);
		numberOfAtoms = newAtomNumber;
		numberOfAtomsx3 = numberOfAtoms * 3;
		numberOfResidues++;
		if (updateBonds) //update the bonds
			determineBonds();
	}


	// Public method to add an atom to a specified strand and residue
	//  The residue number passed in is strand relative
	// All molecule, strand, residue, and atom fields are appropriately updated
	public int addAtom(int strandNumber, int strandResidueNumber, Atom newAtom){
		Residue theResidue = strand[strandNumber].residue[strandResidueNumber];
		return(addAtom(theResidue.moleculeResidueNumber,newAtom));
	}

	// Public method to add an atom to a specified residue
	//  The residue number passed in is molecule relative
	// All molecule, strand, residue, and atom fields are appropriately updated
	public int addAtom(int moleculeResidueNumber, Atom newAtom){
		Residue theResidue = residue[moleculeResidueNumber];
		int newAtomPosition = theResidue.atom[theResidue.numberOfAtoms-1].moleculeAtomNumber+1;
		
		newAtom.moleculeAtomNumber = newAtomPosition;
		
		int newAtomNumber = numberOfAtoms + 1;
		int newAtomNumberx3 = newAtomNumber * 3;
		
		connectivityValid = false;
		connectivity12Valid = false;

		// make room for the new atoms in atom[] and actualCoordinates[]
		expandMolecule(0,1,residue[moleculeResidueNumber].strandNumber,moleculeResidueNumber);
		atom[newAtomPosition]=newAtom;
		int newAtPosx3 = (newAtomPosition+1)*3;
		actualCoordinates[newAtPosx3 - 3] = newAtom.coord[0];
		actualCoordinates[newAtPosx3 - 2] = newAtom.coord[1];
		actualCoordinates[newAtPosx3 - 1] = newAtom.coord[2];

		// strand.addAtom() does additional bookkeeping and will call residue.addAtom
		strand[residue[moleculeResidueNumber].strandNumber].addAtom( 
			residue[moleculeResidueNumber].strandResidueNumber, newAtom);

		numberOfAtomsx3 = newAtomNumberx3;
		return( numberOfAtoms++ );
	}

	// Public method to delete a specified strand and all residue and
	//  atom information contained in the strand
	// All molecule, strand, residue, and atom fields are appropriately updated
	// This function relies upon a proper increasing numerical ordering of the
	//  atoms, residues, and strands. This means for all i < j all atoms in
	//  residue[i] are numbered less than all atoms in residue[j]. Additionally
	//  all residues in strand[i] are numbered less than all residues in strand[j]
	public int deleteStrand(int strandNumber){
		
		connectivityValid = false;
		connectivity12Valid = false;
		
		int newAtomNumber = numberOfAtoms - strand[strandNumber].numberOfAtoms;
		int newAtomNumberx3 = newAtomNumber*3;
		int newResidueNumber = numberOfResidues - strand[strandNumber].numberOfResidues;
		int oldStrandNumAtoms = strand[strandNumber].numberOfAtoms;
		int oldStrandNumResidues = strand[strandNumber].numberOfResidues;
		for(int i=strandNumber+1; i<numberOfStrands; i++)
			strand[i].number -= 1;

		Atom smallerAtomArray[] = new Atom[newAtomNumber];
		float smallerActualCoordinateArray[] = new float[newAtomNumberx3];
		Residue smallerResidueArray[] = new Residue[newResidueNumber];
		Strand smallerStrandArray[] = new Strand[numberOfStrands-1];

		// delete all strand atoms and atom info from Molecule arrays
		for(int i=0, j=0, jx3=0, ix3=0; i<numberOfAtoms; i++){
			if(!(atom[i].strandNumber == strandNumber)){
				ix3 = i*3;
				jx3 = j*3;
				smallerActualCoordinateArray[jx3] = actualCoordinates[ix3];
				smallerActualCoordinateArray[jx3+1] = actualCoordinates[ix3+1];
				smallerActualCoordinateArray[jx3+2] = actualCoordinates[ix3+2];
				smallerAtomArray[j++] = atom[i];
			}
		}
		atom = smallerAtomArray;
		actualCoordinates = smallerActualCoordinateArray;
	
		/* update residue array */
		for(int i=0, j=0; i<numberOfResidues; i++){
			if(!(residue[i].strandNumber == strandNumber))
				smallerResidueArray[j++] = residue[i];
		}
		residue = smallerResidueArray;

		numberOfAtoms = newAtomNumber;
		numberOfAtomsx3 = newAtomNumberx3;
		numberOfResidues = newResidueNumber;
		System.arraycopy(strand, 0, smallerStrandArray, 0, strandNumber);
		if (strandNumber<numberOfStrands-1)
			System.arraycopy(strand, strandNumber + 1, smallerStrandArray,
				strandNumber, strand.length - strandNumber - 1);
		strand = smallerStrandArray;
		
		int curAtNum=0;
		for(int i=0;i<(numberOfStrands-1);i++) {
			for(int j=0;j<strand[i].numberOfResidues;j++) {
				strand[i].residue[j].strandNumber = i;
				if (i>=strandNumber)
					strand[i].residue[j].moleculeResidueNumber -= oldStrandNumResidues;
				for(int k=0;k<strand[i].residue[j].numberOfAtoms;k++) {
					strand[i].residue[j].atom[k].moleculeAtomNumber = curAtNum++;	
					strand[i].residue[j].atom[k].strandNumber = i;
					// This is ok because there should not be bonds between strands
					if (i>=strandNumber) {
						strand[i].residue[j].atom[k].moleculeResidueNumber -= oldStrandNumResidues;
						for(int m=0; m<strand[i].residue[j].atom[k].numberOfBonds; m++)
							strand[i].residue[j].atom[k].bond[m] -= oldStrandNumAtoms;
					}
				}		
			}
		}
		
		establishConnectivity(false);
		return(--numberOfStrands);
	}

	// Public interface for deleting a residue
	public int deleteResidue(int moleculeResidueNumber){

		connectivityValid = false;
		connectivity12Valid = false;
		int strandNumber = residue[moleculeResidueNumber].strandNumber;
		// if this is the only residue in the strand delete the strand
		if(strand[strandNumber].numberOfResidues == 1) {
			return(deleteStrand(strandNumber));
		}

		int oldResNumAtoms = residue[moleculeResidueNumber].numberOfAtoms;
		int newAtomNumber = numberOfAtoms - oldResNumAtoms;
		int newAtomNumberx3 = newAtomNumber * 3;
		int newNumberOfResidues = numberOfResidues - 1;
		int lowAtomIndex = residue[moleculeResidueNumber].atom[0].moleculeAtomNumber;
		int highAtomIndex = lowAtomIndex + oldResNumAtoms-1;
		Atom smallerAtomArray[] = new Atom[newAtomNumber];
		float smallerActualCoordinateArray[] = new float[newAtomNumberx3];
		Residue smallerResidueArray[] = new Residue[newNumberOfResidues];

	    // copy over good atoms and atom info
	    for(int i=0, j=0, jx3=0, ix3=0; i<numberOfAtoms; i++){
			if (!(atom[i].moleculeResidueNumber==moleculeResidueNumber)){
		        ix3 = i*3;
		        jx3 = j*3;
		        smallerActualCoordinateArray[jx3] = actualCoordinates[ix3];
		        smallerActualCoordinateArray[jx3+1] = actualCoordinates[ix3+1];
		        smallerActualCoordinateArray[jx3+2] = actualCoordinates[ix3+2];
		        smallerAtomArray[j++] = atom[i];
			}
		}
	    atom = smallerAtomArray;
		actualCoordinates = smallerActualCoordinateArray;

		// update atom member variables
		for(int i=0; i<newAtomNumber; i++){
			if (atom[i].moleculeResidueNumber > moleculeResidueNumber){
				atom[i].moleculeAtomNumber -= oldResNumAtoms;
				atom[i].moleculeResidueNumber -= 1;
			}
			// only check atoms in another strand, the atoms
			//  in the strand with the residue being deleted
			//  will handle its own atoms
			if (atom[i].strandNumber != strandNumber){
				for(int m=0; m<atom[i].numberOfBonds; m++){
					if (atom[i].bond[m] >= lowAtomIndex){
						if (atom[i].bond[m] <= highAtomIndex){
							atom[i].deleteBond(m);
							m--;
						}
						else
							atom[i].bond[m] -= oldResNumAtoms;
					}
				}
			}
		}

		// delete residue from strand, does additional bookkeeping
		strand[strandNumber].deleteResidue(residue[moleculeResidueNumber].strandResidueNumber);

		System.arraycopy(residue,0,smallerResidueArray,0,moleculeResidueNumber);
		if (moleculeResidueNumber < numberOfResidues-1)
			System.arraycopy(residue, moleculeResidueNumber+1, smallerResidueArray, 
				moleculeResidueNumber, residue.length-moleculeResidueNumber-1);
		residue = smallerResidueArray;

		for(int i=0; i<(numberOfResidues-1); i++){
			if(i>=moleculeResidueNumber){
				residue[i].moleculeResidueNumber -= 1;
			}
		}

		numberOfAtoms = newAtomNumber;
		numberOfAtomsx3 = newAtomNumberx3;

		return(--numberOfResidues);
	}

	// One public interface for deleting atoms
	// Deletes an atom given a residue number and an atom name as a string
	// If multiple atoms match this name in this residue, only the first
	//  one is deleted
	public int deleteAtom(String atomName, int moleculeResidueNumber){
		for(int i=0; i<residue[moleculeResidueNumber].numberOfAtoms; i++)
			if(residue[moleculeResidueNumber].atom[i].name.equals(atomName))
				return(deleteAtom(residue[moleculeResidueNumber].atom[ 
					i].moleculeAtomNumber));
		return(0);
	}

	// One public interface for deleting atoms
	// Deletes an atom given a molecule relative strand number, the strand
	//  relative residue number, and a residue relative atom number
	public int deleteAtom(int strandNumber, int strandResidueNumber, int residueAtomNumber){
		return(deleteAtom(strand[strandNumber].residue[strandResidueNumber].atom[
			residueAtomNumber].moleculeAtomNumber));
	}

	// One public interface for deleting atoms
	// Deletes an atom given a molecule relative residue number and a residue
	//  relative atom number
	public int deleteAtom(int moleculeResidueNumber, int residueAtomNumber){
		return(deleteAtom(residue[moleculeResidueNumber].atom[residueAtomNumber].moleculeAtomNumber));
	}

	// One public interface for deleting atoms
	// Deletes an atom given a molecule residue atom number
	public int deleteAtom(int moleculeAtomNumber){

		connectivityValid = false;
		connectivity12Valid = false;

		int strandNumber = atom[moleculeAtomNumber].strandNumber;
		int moleculeResidueNumber = atom[moleculeAtomNumber].moleculeResidueNumber;
		
		if(residue[moleculeResidueNumber].numberOfAtoms == 1){
			return(deleteResidue(moleculeResidueNumber));
		}

		int newAtomNumber = numberOfAtoms-1;
		int newAtomNumberx3 = newAtomNumber * 3;
		int moleculeAtomNumberx3 = moleculeAtomNumber * 3;
		Atom atomToDelete = atom[moleculeAtomNumber];

    // copy over good atoms and atom info
		Atom smallerAtomArray[] = new Atom[numberOfAtoms-1];
		System.arraycopy(atom, 0, smallerAtomArray, 0, moleculeAtomNumber);
		if (moleculeAtomNumber<newAtomNumber)
			System.arraycopy(atom, moleculeAtomNumber+1, smallerAtomArray, 
				moleculeAtomNumber, atom.length-moleculeAtomNumber-1);
		atom = smallerAtomArray;

		float smallerCoordinateArray[] = new float[newAtomNumberx3];
		System.arraycopy(actualCoordinates, 0, smallerCoordinateArray, 0,
			moleculeAtomNumberx3);
		if (moleculeAtomNumber<newAtomNumber)
			System.arraycopy(actualCoordinates, moleculeAtomNumberx3 + 3, 
				smallerCoordinateArray, moleculeAtomNumberx3, newAtomNumberx3 - 
				moleculeAtomNumberx3); 
		actualCoordinates = smallerCoordinateArray;

		// update atom member variables
		for(int i=0; i<newAtomNumber; i++){
      if (atom[i].moleculeAtomNumber > moleculeAtomNumber){
				atom[i].moleculeAtomNumber -= 1;
			}
			if (atom[i].strandNumber != strandNumber){
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
		}

		// delete residue from strand, does additional bookkeeping
		strand[strandNumber].deleteAtom(atomToDelete.strandResidueNumber,
			atomToDelete.residueAtomNumber);

		numberOfAtoms = newAtomNumber;	
		numberOfAtomsx3 = numberOfAtoms * 3;
		return(numberOfAtoms);
	}


	// This function deletes both directions of the bond between
	//  the two specified atoms. Both atomi and atomj are
	//  molecule based atom numbers
	public void deleteBondBetween(int atomi, int atomj){
	
		connectivity12Valid = false;
		connectivityValid = false;
		
		atom[atomi].deleteBondTo(atomj);
		atom[atomj].deleteBondTo(atomi);
	}


	// This functions adds both directions of a bond between
	//  two specificed atoms. Both atomi and atomj are
	//  molecule based atom numbers
	public void addBondBetween(int atomi, int atomj){
	
		connectivity12Valid = false;
		connectivityValid = false;
		
		atom[atomi].addBond(atomj);
		atom[atomj].addBond(atomi);
	}
	

	// This function copies the molecule local actualCoordinates back to the
	//  corresponding atoms. The actualCoordinates array is a local temporary
	//  copy of the atom coordinates which allows us to muck with the coordinates
	//  without directly involving the atoms
	public void resolveCoordinates(){
		int ix3 = -3;
		for(int i=0; i<numberOfAtoms; i++){
			ix3 += 3;
			atom[i].coord[0] = actualCoordinates[ix3];
			atom[i].coord[1] = actualCoordinates[ix3+1];
			atom[i].coord[2] = actualCoordinates[ix3+2];
		}
	}

	// As before but copies coordinates back only for the specified atom
	public void resolveCoordinates(int atomNumber){
		int atomNumberx3 = atomNumber * 3;
		atom[atomNumber].coord[0] = actualCoordinates[atomNumberx3];
		atom[atomNumber].coord[1] = actualCoordinates[atomNumberx3 + 1];
		atom[atomNumber].coord[2] = actualCoordinates[atomNumberx3 + 2];
	}

	// This function copies the coordinates of the specified atom
	//  into the local actualCoorinates array. The atom itself
	//  is passed in
	public void updateCoordinates(Atom updateAtom){
		int atomNumberx3 = updateAtom.moleculeAtomNumber * 3;
			
		actualCoordinates[atomNumberx3] = updateAtom.coord[0];	
		actualCoordinates[atomNumberx3 + 1] = updateAtom.coord[1];	
		actualCoordinates[atomNumberx3 + 2] = updateAtom.coord[2];	
	}

	// This function copies the coordinates of the specified atom
	//  into the local actualCoorinates array. The atom number
	//  is passed in
	public void updateCoordinates(int atomNumber){
		Atom updateAtom = atom[atomNumber];
		int atomNumberx3 = atomNumber * 3;
		
		actualCoordinates[atomNumberx3] = updateAtom.coord[0];	
		actualCoordinates[atomNumberx3 + 1] = updateAtom.coord[1];	
		actualCoordinates[atomNumberx3 + 2] = updateAtom.coord[2];	
	}

	// This function copies the coordinates of the specified atom
	//  into the local actualCoorinates array. All atoms are
	//  updated
	public void updateCoordinates(){
		for(int i=0;i<numberOfAtoms;i++)
			updateCoordinates(atom[i]);
	}

	// Goes through each atom and makes sure its moleculeResidueNumber points to the right residue
	public void updateResidueNumbers() {
		for(int q=0;q<this.numberOfResidues;q++) {
			for(int w=0;w<residue[q].numberOfAtoms;w++)
				residue[q].atom[w].moleculeResidueNumber = q;
		}
	}


	// Goes through and updates the number of atoms in each strand as well as
	//  the number of atoms in the molecule as a whole
	public void updateNumAtoms() {
		numberOfAtoms = 0;
		for(int i=0;i<numberOfStrands;i++){
			strand[i].numberOfAtoms = 0;
			for(int j=0;j<strand[i].numberOfResidues;j++){
				strand[i].numberOfAtoms += strand[i].residue[j].numberOfAtoms;
			}
			numberOfAtoms += strand[i].numberOfAtoms;
		}
	}
	
	// Goes through each atom and updates the moleculeatomnumber
	//  used when reading pdb files or such and we need to number all the
	//  atoms we've just read in
	public void updateMoleculeAtomNumbers() {
		for(int q=0;q<this.numberOfAtoms;q++) {
			atom[q].moleculeAtomNumber = q;
		}
	}

	// Complete bonds created by amino acid template. Bonds
	//  that are only unidirectional are made bidirectional
	public void completeBonds(int resnum) {
		for (int i=0;i<residue[resnum].numberOfAtoms;i++) {
			int curAtNum = residue[resnum].atom[i].moleculeAtomNumber;
			for (int j=0;j<residue[resnum].atom[i].numberOfBonds;j++) {
				if (!(atom[residue[resnum].atom[i].bond[j]].bondedTo(curAtNum)))
					atom[residue[resnum].atom[i].bond[j]].addBond(curAtNum);
			}
		}
	}

	// Adds a bond from atomA of residueA to atomB of residueB. The residue
	//  numbering is molecule based.
	public void connectResidues(int residueA, String atomAName, int residueB, String atomBName){
		Atom atomA = null, atomB = null;
		for(int i=0; i<residue[residueA].numberOfAtoms; i++){
			if (residue[residueA].atom[i].name.equals(atomAName)){
				atomA = residue[residueA].atom[i];
			}
		}
		for(int j=0; j<residue[residueB].numberOfAtoms; j++){
			if (residue[residueB].atom[j].name.equals(atomBName)){
				atomB = residue[residueB].atom[j];
			}
		}
		if ((atomA != null) && (atomB != null)) {
			atomA.addBond(atomB.moleculeAtomNumber);
			atomB.addBond(atomA.moleculeAtomNumber);
		}
	}

	// This function sets a specified torsion by applying the appropriate
	//  rotation to specified atoms
	// aXnum are the four atoms of the dihedral, atom4 is distal and is rotated
	// torsion is the new angle
	// atomList is a list of atoms to change
	// numberOfElements is the number of elements in atomList
	public void setTorsion(int a1num, int a2num, int a3num, int a4num, double torsion, 
		int atomList[], int numberOfElements){
		setTorsion(a1num,a2num,a3num,a4num,torsion,atomList,numberOfElements,true);
	}

	public void changeTorsion(int a1num, int a2num, int a3num, int a4num, double deltaTorsion, 
		int atomList[], int numberOfElements){
		changeTorsion(a1num,a2num,a3num,a4num,deltaTorsion,atomList,numberOfElements,true);
	}

	// This function sets a specified torsion by applying the appropriate
	//  rotation to specified atoms
	// aXnum are the four atoms of the dihedral, atom4 is distal and is rotated
	// torsion is the new angle
	// atomList is a list of atoms to change
	// numberOfElements is the number of elements in atomList
	// If updateAtoms is true the actualCoordinates are copied
	//  back into the atom coordinates
	// This function does some crazy stuff with atom coordinates based on the
	//  functions it needs to call and where they expect coordinates to be.
	// First the four atoms atom.coord are stored, then the atom.coords are set
	//  to the actualCoordinate values, then the torsion is done on all the
	//  atom's acutalCoordinates (this requires atom.coords for the 4 atoms
	//  which is why resolveCoordinates is called), finally if we don't want
	//  to updateAtom coords then the atom.coords are restored
	public void setTorsion(int a1num, int a2num, int a3num, int a4num, double torsion, 
		int atomList[], int numberOfElements, boolean updateAtoms){

		// Aligns dihedral on the +x-axis, modifies the angle then rotates it back

		float backupCoords[][] = new float[4][3];

		Atom a1 = atom[a1num];  
		Atom a2 = atom[a2num];
		Atom a3 = atom[a3num];
		Atom a4 = atom[a4num];

		backupCoords[0][0] = a1.coord[0];
		backupCoords[0][1] = a1.coord[1];
		backupCoords[0][2] = a1.coord[2];
		backupCoords[1][0] = a2.coord[0];
		backupCoords[1][1] = a2.coord[1];
		backupCoords[1][2] = a2.coord[2];
		backupCoords[2][0] = a3.coord[0];
		backupCoords[2][1] = a3.coord[1];
		backupCoords[2][2] = a3.coord[2];
		backupCoords[3][0] = a4.coord[0];
		backupCoords[3][1] = a4.coord[1];
		backupCoords[3][2] = a4.coord[2];

		resolveCoordinates(a1num);  // update a1 coords
		resolveCoordinates(a2num);  // update a2 coords
		resolveCoordinates(a3num);  // update a3 coords
		resolveCoordinates(a4num);  // update a4 coords
		
		Atom pseudoAtom2 = new Atom("a2", a2.coord[0] - a3.coord[0],
			a2.coord[1] - a3.coord[1], a2.coord[2] - a3.coord[2]);
		double originalTorsion = a4.torsion(a1, a2, a3);

		changeTorsion(a1num, a2num, a3num, a4num, torsion-originalTorsion,atomList,
			numberOfElements, updateAtoms);
	
		if (!updateAtoms){
			// restore the 4 atoms coordinates
			a1.coord[0] = backupCoords[0][0];	
			a1.coord[1] = backupCoords[0][1];
			a1.coord[2] = backupCoords[0][2];
			a2.coord[0] = backupCoords[1][0];
			a2.coord[1] = backupCoords[1][1];
			a2.coord[2] = backupCoords[1][2];
			a3.coord[0] = backupCoords[2][0];
			a3.coord[1] = backupCoords[2][1];
			a3.coord[2] = backupCoords[2][2];
			a4.coord[0] = backupCoords[3][0];
			a4.coord[1] = backupCoords[3][1];
			a4.coord[2] = backupCoords[3][2];
		}
	}


	// This function changes a specified torsion by applying the appropriate
	//  rotation to specified atoms
	// See the comment description for setTorsion
	public void changeTorsion(int a1num, int a2num, int a3num, int a4num, double deltaTorsion, 
		int atomList[], int numberOfElements, boolean updateAtoms){

		// Aligns dihedral on the +x-axis, modifies the angle then rotates it back
		int at2x3 = a2num * 3;
		int at3x3 = a3num * 3;
		int at4x3 = a4num * 3;
		Atom pseudoAtom2 = new Atom("a2", actualCoordinates[at2x3] - actualCoordinates[at3x3],
																	actualCoordinates[at2x3+1] - actualCoordinates[at3x3+1],
																	actualCoordinates[at2x3+2] - actualCoordinates[at3x3+2]);
		int numberOfCoordinates = 0;
		int atomListLength = 0;
		int atomNumber;
		if (numberOfElements == 0){
			atomListLength = 0;
			numberOfCoordinates = 1;
		}
		else{
			atomListLength = numberOfElements;
			numberOfCoordinates = atomListLength + 1;
		}
		int numberOfCoordinatesx3 = numberOfCoordinates * 3;
		// array of atom numbers to be rotated
		float temporaryCoordinates[] = new float[numberOfCoordinatesx3];
		// store atom 4 as first entry in array
		temporaryCoordinates[0] = actualCoordinates[at4x3] - actualCoordinates[at3x3];
		temporaryCoordinates[1] = actualCoordinates[at4x3+1] - actualCoordinates[at3x3+1];
		temporaryCoordinates[2] = actualCoordinates[at4x3+2] - actualCoordinates[at3x3+2];
		// translate all other atoms by -a3
		for (int i = 0, ix3, atx3; i<atomListLength; i++){
			atomNumber = atomList[i];	
			if ((atomNumber == a1num) || (atomNumber == a2num) || 
				(atomNumber == a3num) || (atomNumber == a4num))
				continue;
			atx3 = atomNumber * 3;   
			ix3 = i * 3;
			temporaryCoordinates[ix3 + 3] = actualCoordinates[atx3] - 
				actualCoordinates[at3x3];
			temporaryCoordinates[ix3 + 4] = actualCoordinates[atx3 + 1] - 
				actualCoordinates[at3x3+1];
			temporaryCoordinates[ix3 + 5] = actualCoordinates[atx3 + 2] - 
				actualCoordinates[at3x3+2];
		}
		// get angle info from atom2
		RotMatrix theRotation = new RotMatrix();
		double a2YAngle = pseudoAtom2.angleAboutYAxis();
		theRotation.yAxisRotate( -a2YAngle, pseudoAtom2.coord, 1);
		double a2ZAngle = pseudoAtom2.angleAboutZAxis();
		
		theRotation.yAxisRotate( -a2YAngle, temporaryCoordinates, numberOfCoordinates);
		theRotation.zAxisRotate( 180 - a2ZAngle, temporaryCoordinates, numberOfCoordinates);
		theRotation.xAxisRotate( deltaTorsion, temporaryCoordinates, numberOfCoordinates);
		theRotation.zAxisRotate( a2ZAngle -180, temporaryCoordinates, numberOfCoordinates);
		theRotation.yAxisRotate( a2YAngle, temporaryCoordinates, numberOfCoordinates);
		theRotation.translate(actualCoordinates[at3x3], actualCoordinates[at3x3+1],
			actualCoordinates[at3x3+2], temporaryCoordinates, numberOfCoordinates);
			
		for(int i=0, ix3, atx3; i<atomListLength; i++){
			atomNumber = atomList[i];	
			if ((atomNumber == a1num) || (atomNumber == a2num) || 
				(atomNumber == a3num ) || (atomNumber == a4num))
				continue;
			atx3 = atomNumber * 3;
			ix3 = i * 3;
			
			actualCoordinates[atx3] = temporaryCoordinates[ix3 + 3];
			actualCoordinates[atx3 + 1] = temporaryCoordinates[ix3 + 4];
			actualCoordinates[atx3 + 2] = temporaryCoordinates[ix3 + 5];
			if (updateAtoms)
				resolveCoordinates(atomNumber);
		}
		actualCoordinates[at4x3] = temporaryCoordinates[0];
		actualCoordinates[at4x3 + 1] = temporaryCoordinates[1]; 
		actualCoordinates[at4x3 + 2] = temporaryCoordinates[2]; 
		if (updateAtoms)
			resolveCoordinates(a4num);
	}


	// Returns the geometric center of a molecule *not* center of
	//  mass although for large molecules they are essentially
	//  the same
	public float[] getCenter(){
		float xCenter = 0.0f, yCenter = 0.0f, zCenter = 0.0f;
		numberOfAtomsx3 = numberOfAtoms * 3;
		for(int i=0; i<numberOfAtoms; i+=3){
			xCenter += atom[i].coord[0];
			yCenter += atom[i].coord[1];
			zCenter += atom[i].coord[2];
		}
		xCenter /= numberOfAtoms;	
		yCenter /= numberOfAtoms;	
		zCenter /= numberOfAtoms;	
		float[] cent = new float[3];
		cent[0] = xCenter;
		cent[1] = yCenter;
		cent[2] = zCenter;
		return cent;
	}

	// Returns the center of mass of the molecule
	public float[] getCenterOfMass(){
		float centOfMass[] = new float[3];
		float xCenter = 0.0f, yCenter = 0.0f, zCenter = 0.0f;
		double totalMass = 0.0;
		double tmpMass = 0.0;
		
		for(int j=0;j<numberOfResidues;j++) {
			for(int i=0;i<residue[j].numberOfAtoms;i++){
				tmpMass = residue[j].atom[i].mass;
				totalMass += tmpMass;
				xCenter += residue[j].atom[i].coord[0] * tmpMass;
				yCenter += residue[j].atom[i].coord[1] * tmpMass;
				zCenter += residue[j].atom[i].coord[2] * tmpMass;
			}
		}
		xCenter /= totalMass;
		yCenter /= totalMass;
		zCenter /= totalMass;
		centOfMass[0] = xCenter;
		centOfMass[1] = yCenter;
		centOfMass[2] = zCenter;
		return(centOfMass);
	}


	// Returns the center of mass of the given strand
	public double[] getStrandCOM(int strNum){
		
		double centOfMass[] = new double[3];
		double xCenter = 0.0, yCenter = 0.0, zCenter = 0.0;
		double totalMass = 0.0;
		double tmpMass = 0.0;
		
		int baseNum = strand[strNum].residue[0].atom[0].moleculeAtomNumber * 3;
		for(int j=0;j<strand[strNum].numberOfResidues;j++) {
			for(int i=0;i<strand[strNum].residue[j].numberOfAtoms;i++){
				tmpMass = strand[strNum].residue[j].atom[i].mass;
				totalMass += tmpMass;
				xCenter += actualCoordinates[baseNum] * tmpMass;
				yCenter += actualCoordinates[baseNum+1] * tmpMass;
				zCenter += actualCoordinates[baseNum+2] * tmpMass;
				baseNum += 3;
			}
		}
		xCenter /= totalMass;
		yCenter /= totalMass;
		zCenter /= totalMass;
		centOfMass[0] = xCenter;
		centOfMass[1] = yCenter;
		centOfMass[2] = zCenter;
		return(centOfMass);
	}
	
	// Returns the center of mass of the given residue (molecule-relative numbering)
	public double[] getResidueCOM(int resNum){
		
		double centOfMass[] = new double[3];
		double xCenter = 0.0, yCenter = 0.0, zCenter = 0.0;
		double totalMass = 0.0;
		double tmpMass = 0.0;
		
		int baseNum = residue[resNum].atom[0].moleculeAtomNumber * 3;
		for(int i=0;i<residue[resNum].numberOfAtoms;i++){
			tmpMass = residue[resNum].atom[i].mass;
			totalMass += tmpMass;
			xCenter += actualCoordinates[baseNum] * tmpMass;
			yCenter += actualCoordinates[baseNum+1] * tmpMass;
			zCenter += actualCoordinates[baseNum+2] * tmpMass;
			baseNum += 3;
		}
		
		xCenter /= totalMass;
		yCenter /= totalMass;
		zCenter /= totalMass;
		centOfMass[0] = xCenter;
		centOfMass[1] = yCenter;
		centOfMass[2] = zCenter;
		return(centOfMass);
	}
	
	//Computes the geometric center (from the actualCoordinates[]) for residue resNum (molecule-relative numbering)
	public double [] getResidueGC(int resNum){
		
		Residue r = residue[resNum];
		
		float x=0, y=0, z=0;
		for(int i=0; i<r.numberOfAtoms; i++){
			int curAtom = r.atom[i].moleculeAtomNumber;
			x += actualCoordinates[curAtom*3];
			y += actualCoordinates[curAtom*3+1];
			z += actualCoordinates[curAtom*3+2];
		}
		
		double com[] = new double[3];
		com[0] = (x / r.numberOfAtoms);
		com[1] = (y / r.numberOfAtoms);
		com[2] = (z / r.numberOfAtoms);
		return com;
	}

	// Translate the specified residue by dx, dy, dz
	// resNum is molecule relative
	public void translateResidue(int resNum, float dx, float dy,
		float dz){
		translateResidue(resNum,dx,dy,dz,true);
	}
	
	// Translate the specified residue by dx, dy, dz
	// resNum is molecule relative
	// If updateAtoms is true the actualCoordinates are copied
	//  back into the atom coordinates
	public void translateResidue(int resNum, float dx, float dy,
		float dz, boolean updateAtoms){
		translateResidue(residue[resNum].strandNumber,
			residue[resNum].strandResidueNumber,dx,dy,dz,updateAtoms);
	}
	// Translate the specified residue by dx, dy, dz
	// resNum is strand relative
	public void translateResidue(int strNum, int resNum, float dx,
		float dy, float dz){
		translateResidue(strNum,resNum,dx,dy,dz,true);		
	}
		
	// Translate the specified residue by dx, dy, dz
	// resNum is strand relative
	// If updateAtoms is true the actualCoordinates are copied
	//  back into the atom coordinates
	public void translateResidue(int strNum, int resNum, float dx,
		float dy, float dz, boolean updateAtoms){
		
		int index = strand[strNum].residue[resNum].atom[0].moleculeAtomNumber * 3;
		for(int i=0;i<strand[strNum].residue[resNum].numberOfAtoms;i++){
			actualCoordinates[index++] += dx;
			actualCoordinates[index++] += dy;
			actualCoordinates[index++] += dz;
			if (updateAtoms)
				resolveCoordinates(strand[strNum].residue[resNum].atom[i].moleculeAtomNumber);
		}
	}

	// Translate the specified strand by dx, dy, dz
	public void translateStrand(int strNum, float dx, float dy,
		float dz){
		translateStrand(strNum,dx,dy,dz,true);
	}

	// Translate the specified strand by dx, dy, dz
	// If updateAtoms is true the actualCoordinates are copied
	//  back into the atom coordinates
	public void translateStrand(int strNum, float dx, float dy,
		float dz, boolean updateAtoms){
		for(int i=0;i<strand[strNum].numberOfResidues;i++){
			translateResidue(strNum,i,dx,dy,dz,updateAtoms);
		}		
	}
	
	// Translates the molecule by the amount dx, dy, dz
	public void translateMolecule(float dx, float dy, float dz){
		translateMolecule(dx,dy,dz,true);
	}	
	
	// Translates the molecule by the amount dx, dy, dz
	// If updateAtoms is true the actualCoordinates are copied
	//  back into the atom coordinates
	public void translateMolecule(float dx, float dy, float dz, boolean
		updateAtoms){
		for(int i=0; i<numberOfAtomsx3; i=i+3){
			actualCoordinates[i] += dx;
			actualCoordinates[i+1] += dy;
			actualCoordinates[i+2] += dz;
		}
		if (updateAtoms)
			resolveCoordinates();
	}


	public void rotateStrandAroundCOM(int ligStrNum, float dx, float dy,
		float dz, float theta){
		rotateStrandAroundCOM(ligStrNum,dx,dy,dz,theta,true);
	}

	public void rotateStrandAroundCOM(int ligStrNum, float dx, float dy,
		float dz, float theta, boolean updateAtoms){
		
		double theCOM[] = getStrandCOM(ligStrNum);
		
		rotateStrand(ligStrNum,dx,dy,dz,theCOM[0],theCOM[1],theCOM[2],
			theta, updateAtoms);
	}

	public void rotateStrand(int ligStrNum, float dx, float dy, float dz,
		double cx, double cy, double cz, float theta){
		
		rotateStrand(ligStrNum,dx,dy,dz,cx,cy,cz,theta,true);
	}
	
	public void rotateStrand(int ligStrNum, float dx, float dy, float dz,
		double cx, double cy, double cz, float thetaDeg, boolean updateAtoms){

		double tx,ty,tz;

		float[][] rot_mtx = new float[3][3];
		RotMatrix rM = new RotMatrix();
		rM.getRotMatrix(dx,dy,dz,thetaDeg,rot_mtx);

		Strand theStrand = strand[ligStrNum];
		int i=theStrand.residue[0].atom[0].moleculeAtomNumber;
		int baseNum = i * 3;
		
		for(int w=0;w<theStrand.numberOfResidues;w++) {
			for(int q=0;q<theStrand.residue[w].numberOfAtoms;q++) {
				tx=actualCoordinates[baseNum] - cx;
				ty=actualCoordinates[baseNum+1] - cy;
				tz=actualCoordinates[baseNum+2] - cz;

				actualCoordinates[baseNum] = (float)(tx * rot_mtx[0][0] + ty * rot_mtx[0][1] + tz * rot_mtx[0][2] + cx);
				actualCoordinates[baseNum+1] = (float)(tx * rot_mtx[1][0] + ty * rot_mtx[1][1] + tz * rot_mtx[1][2] + cy);
				actualCoordinates[baseNum+2] = (float)(tx * rot_mtx[2][0] + ty * rot_mtx[2][1] + tz * rot_mtx[2][2] + cz);

				if (updateAtoms)
					resolveCoordinates(i);
				baseNum += 3;
				i++;
			}
		}
	}
	
	//Rotates the given residue (molecule-relative numbering) by thetaDeg degrees
	//		around axis dx, dy, dz (around the point cs,cy,cz)
	public void rotateResidue(int resNum, float dx, float dy, float dz, 
			double cx, double cy, double cz, float thetaDeg, boolean updateAtoms){

		double tx,ty,tz;

		float[][] rot_mtx = new float[3][3];
		RotMatrix rM = new RotMatrix();
		rM.getRotMatrix(dx,dy,dz,thetaDeg,rot_mtx);

		int i = residue[resNum].atom[0].moleculeAtomNumber;
		int baseNum = i * 3;
		
		for(int q=0;q<residue[resNum].numberOfAtoms;q++) {
			
			tx=actualCoordinates[baseNum] - cx;
			ty=actualCoordinates[baseNum+1] - cy;
			tz=actualCoordinates[baseNum+2] - cz;

			actualCoordinates[baseNum] = (float)(tx * rot_mtx[0][0] + ty * rot_mtx[0][1] + tz * rot_mtx[0][2] + cx);
			actualCoordinates[baseNum+1] = (float)(tx * rot_mtx[1][0] + ty * rot_mtx[1][1] + tz * rot_mtx[1][2] + cy);
			actualCoordinates[baseNum+2] = (float)(tx * rot_mtx[2][0] + ty * rot_mtx[2][1] + tz * rot_mtx[2][2] + cz);

			if (updateAtoms)
				resolveCoordinates(i);
			baseNum += 3;
			i++;
		}
	}
	
	//Rotates the given atomList (molecule-relative numbering) by thetaDeg degrees
	//		around axis dx, dy, dz (around the point cs,cy,cz)
	public void rotateAtomList(int atomList[], float dx, float dy, float dz, 
			double cx, double cy, double cz, float thetaDeg, boolean updateAtoms){

		double tx,ty,tz;

		float[][] rot_mtx = new float[3][3];
		RotMatrix rM = new RotMatrix();
		rM.getRotMatrix(dx,dy,dz,thetaDeg,rot_mtx);
		
		for(int q=0;q<atomList.length;q++) {
			
			int baseNum = atomList[q]*3;
			
			tx=actualCoordinates[baseNum] - cx;
			ty=actualCoordinates[baseNum+1] - cy;
			tz=actualCoordinates[baseNum+2] - cz;

			actualCoordinates[baseNum] = (float)(tx * rot_mtx[0][0] + ty * rot_mtx[0][1] + tz * rot_mtx[0][2] + cx);
			actualCoordinates[baseNum+1] = (float)(tx * rot_mtx[1][0] + ty * rot_mtx[1][1] + tz * rot_mtx[1][2] + cy);
			actualCoordinates[baseNum+2] = (float)(tx * rot_mtx[2][0] + ty * rot_mtx[2][1] + tz * rot_mtx[2][2] + cz);

			if (updateAtoms)
				resolveCoordinates(atomList[q]);
		}
	}
	
	//Rotates the given atom i (molecule-relative numbering) by thetaDeg degrees
	//		around axis dx, dy, dz (around the point cs,cy,cz)
	public void rotateAtom(int i, float dx, float dy, float dz, 
			double cx, double cy, double cz, float thetaDeg, boolean updateAtoms){

		double tx,ty,tz;
		int baseNum = i*3;

		float[][] rot_mtx = new float[3][3];
		RotMatrix rM = new RotMatrix();
		rM.getRotMatrix(dx,dy,dz,thetaDeg,rot_mtx);
			
		tx=actualCoordinates[baseNum] - cx;
		ty=actualCoordinates[baseNum+1] - cy;
		tz=actualCoordinates[baseNum+2] - cz;

		actualCoordinates[baseNum] = (float)(tx * rot_mtx[0][0] + ty * rot_mtx[0][1] + tz * rot_mtx[0][2] + cx);
		actualCoordinates[baseNum+1] = (float)(tx * rot_mtx[1][0] + ty * rot_mtx[1][1] + tz * rot_mtx[1][2] + cy);
		actualCoordinates[baseNum+2] = (float)(tx * rot_mtx[2][0] + ty * rot_mtx[2][1] + tz * rot_mtx[2][2] + cz);

		if (updateAtoms)
			resolveCoordinates(i);
	}

	
	// This function rotates the entire molecule by thetaDeg
	//  degrees around axis dx, dy, dz (around the center of mass)
	public void rotateAroundCOM(double dx, double dy, double dz,
		double thetaDeg) {

		// Get center of mass
		float[] centOfMass = getCenterOfMass();
	
		rotateMolecule(dx,dy,dz,thetaDeg,centOfMass[0],centOfMass[1],centOfMass[2]);
	}	

	
	// This function rotates the entire molecule by thetaDeg
	//  degrees around axis dx, dy, dz (around the point cx,cy,cz)
	public void rotateMolecule(double dx, double dy, double dz,
		double thetaDeg, float cx, float cy, float cz) {
	
		float fx,fy,fz, tx,ty,tz;
		fx = (new Double(dx)).floatValue();
		fy = (new Double(dy)).floatValue();
		fz = (new Double(dz)).floatValue();

		float[][] rot_mtx = new float[3][3];
		RotMatrix rM = new RotMatrix();
		rM.getRotMatrix(fx,fy,fz,(float) thetaDeg,rot_mtx);

		for(int w=0;w<numberOfResidues;w++) {
			for(int q=0;q<residue[w].numberOfAtoms;q++) {
				tx=residue[w].atom[q].coord[0] - cx;
				ty=residue[w].atom[q].coord[1] - cy;
				tz=residue[w].atom[q].coord[2] - cz;
				
				residue[w].atom[q].coord[0] = tx * rot_mtx[0][0] + ty * rot_mtx[0][1] + tz * rot_mtx[0][2] + cx;
				residue[w].atom[q].coord[1] = tx * rot_mtx[1][0] + ty * rot_mtx[1][1] + tz * rot_mtx[1][2] + cy;
				residue[w].atom[q].coord[2] = tx * rot_mtx[2][0] + ty * rot_mtx[2][1] + tz * rot_mtx[2][2] + cz;

				updateCoordinates(residue[w].atom[q]);
			}
		}
	}


	// Places the geometric center of the molecule at the origin
	public void center(){
		float theCent[];
		theCent = getCenter();
		translateMolecule(-theCent[0],-theCent[1],-theCent[2]);
	}

	// Places the center of mass of the molecule at the origin
	public void centerByMass(){
		float theCOM[];
		theCOM = getCenterOfMass();
		translateMolecule(-theCOM[0],-theCOM[1],-theCOM[2]);
	}
	
	// Marks the bondedMatrix for atoms i and j
	// If i < j an entry is put at [i][j] else it's added at [j][i]
	private void markBonded(int i, int j){
		if(i<j)
			bondedMatrix[i][j]=true;
		else
			bondedMatrix[j][i]=true;
	}

	// This function updates the connected, connected12, connected13,
	//  and connected14 arrays based on atom connectivity.
	// If just12 == true then only the 12 connections are computed
	// *Note that connectivity calculations check single, double, and triple bonds.
	public void establishConnectivity(boolean just12){
		// make lists of 1-2, 1-3, 1-4 connectivities for refinement
		connected = new int[numberOfAtoms][MAX_ATOM_CONNECTIONS];
		connected12 = new int[1];
		numberOf12Connections = 0;
		connected13 = new int[1];
		numberOf13Connections = 0;
		connected14 = new int[1];
		numberOf14Connections = 0;
		numberNonBonded = 0;
		nonBonded = new int[NB_INCREMENT_SIZE];
		int atomi = -1, atomj = -1, atomk = -1, atomx = -1;

		if(!just12){
			// Allocate the bondedMatrix
			if (bondedMatrixSize < numberOfAtoms){
				bondedMatrix = new boolean[numberOfAtoms][numberOfAtoms];
				bondedMatrixSize = numberOfAtoms;
			}
			// Clear the bonded matrix
			for(int i=0;i<numberOfAtoms;i++){
				for(int j=i+1;j<numberOfAtoms;j++){
					bondedMatrix[i][j] = false;
				}
			}
		}

		if(debug)
			System.out.println("Establishing connectivity...");

		// update 1-2 connections
		for(int i=0; i<numberOfAtoms; i++){
			for(int j=i+1; j<numberOfAtoms; j++){
				if(isConnected(i,j)){
					connect12(i,j);
					connect(i,j);
					connect(j,i);
					if (!just12)
						markBonded(i,j);
				}
			}
		}

		// trim down the connected array to the actual size used
		int smallerConnectArray[] = new int[numberOf12Connections * 2];
		System.arraycopy(connected12, 0, smallerConnectArray, 0, numberOf12Connections * 2);
		connected12 = smallerConnectArray;

		connectivity12Valid = true;
		
		if (just12)
			return;

		// update 1-3 connections, chain off existing 1-2 connections
		int numberOf12Connectionsx2 = numberOf12Connections * 2 - 1;
		for(int i=0; i<numberOf12Connectionsx2; i=i+2){
			atomi = connected12[i];
			atomj = connected12[i+1];
			for(int k=1;k<connected[atomi][0]+1;k++){
				atomx = connected[atomi][k];
				if(atomx != atomj){
					connect13(atomx, atomi, atomj);
					markBonded(atomj,atomx);
				}
			}
			for(int k=1;k<connected[atomj][0]+1;k++){
				atomx = connected[atomj][k];
				if(atomx != atomi){
					connect13(atomx, atomj, atomi);
					markBonded(atomi,atomx);
				}
			}
		}

		// trim down the connected array to the actual size used
		smallerConnectArray = new int[numberOf13Connections * 3];
		System.arraycopy(connected13, 0, smallerConnectArray, 0, numberOf13Connections * 3);
		connected13 = smallerConnectArray;

		// update 1-4 connections
		int numberOf13Connectionsx3 = numberOf13Connections * 3 - 2;
		for(int i=0; i<numberOf13Connectionsx3; i=i+3){
			atomi = connected13[i];
			atomj = connected13[i+1];
			atomk = connected13[i+2];
			for(int l=1;l<connected[atomi][0]+1;l++){
				atomx = connected[atomi][l];
				if ((atomx != atomj) && (atomx != atomk)){
					if (isNOT13Connected(atomx,atomk) && isNOTAlready14Connected(atomx,atomk)){
						connect14(atomx, atomi, atomj, atomk);
						markBonded(atomx,atomk);
					}
				}
			}
			for(int l=1;l<connected[atomk][0]+1;l++){
				atomx = connected[atomk][l];
				if ((atomx != atomi) && (atomx != atomj)){
					if (isNOT13Connected(atomx,atomi) && isNOTAlready14Connected(atomx,atomi)){
						connect14(atomi, atomj, atomk, atomx);
						markBonded(atomi,atomx);
					}
				}
			}
		}

		// trim down the connected array to the actual size used
		smallerConnectArray = new int[numberOf14Connections * 4];
		System.arraycopy(connected14, 0, smallerConnectArray, 0, numberOf14Connections * 4);
		connected14 = smallerConnectArray;

		// Search through all pairs for nonbonded pairs
		int numberOfNBConnectionsx2 = 0;
		for(int i=0; i<numberOfAtoms; i++) {
			for(int j=i+1; j<numberOfAtoms; j++) {
				if (bondedMatrix[i][j] == false) {
					numberNonBonded++;
					numberOfNBConnectionsx2 = numberNonBonded*2;
					if (numberOfNBConnectionsx2 >= nonBonded.length) {
						int largerConnectArray[] = new int[numberOfNBConnectionsx2 + NB_INCREMENT_SIZE*2];
						System.arraycopy(nonBonded, 0, largerConnectArray, 0, nonBonded.length);
						nonBonded = largerConnectArray;
					}
					nonBonded[numberOfNBConnectionsx2 - 2] = i;
					nonBonded[numberOfNBConnectionsx2 - 1] = j;
				}
			}
		}

		// trim down the nonbonded array to the actual size used
		smallerConnectArray = new int[numberOfNBConnectionsx2];
		System.arraycopy(nonBonded, 0, smallerConnectArray, 0, numberOfNBConnectionsx2);
		nonBonded = smallerConnectArray;
		
		connectivityValid = true;
	}

	// Adds an entry to the atom1 row in the connected array
	//  indiating its connection to atom2
	private void connect(int atom1, int atom2){
		if (debug) {
			if (connected[atom1][0] >= 6){
				System.out.println("Debug: Error in molecule.connect. Atom " + atom1 + " has too many bonds");
				System.out.println("molecAtNum: " +atom[atom1].moleculeAtomNumber + " resAtNum: " + atom[atom1].residueAtomNumber);
			}
		}
		int numberOfConnections = connected[atom1][0] + 1;
		connected[atom1][numberOfConnections] = atom2;
		connected[atom1][0] = numberOfConnections;
	}

	// Returns true if atom1 is bonded to atom2, if they
	//  aren't bonded then false is returned
	public boolean isConnected(int atom1, int atom2){
		
		if ( atom[atom1].bondedTo(atom2) )
			return(true);
		else if ( atom[atom2].bondedTo(atom1) )
			return(true);
		
		return(false);
	}

	// Returns true if atom1 and atom2 are not 1-3 connected
	// Looking through two levels of direct connectivity is faster
	//  than looking though the whole connected13 array
	private boolean isNOT13Connected(int atom1, int atom2){
		int tmpAt = 0;
		for(int i=0; i<atom[atom1].numberOfBonds; i++){
			tmpAt = atom[atom1].bond[i];
			for(int j=0; j<atom[tmpAt].numberOfBonds; j++){
				if(atom2 == atom[tmpAt].bond[j])
					return false;
			}
		}
		return true;
	}

	// Returns true if atom1 and atom2 are not 1-4 connected
	// Here we can't look through three levels of direct connectivity
	//  because we essentially just want to know if these atoms
	//  are already in the 1-4 connected array, we know that they are
	//  1-4 connected, but have they been found yet
	private boolean isNOTAlready14Connected(int atom1, int atom2){
		int ix4=0;
		for(int i=0;i<numberOf14Connections;i++){
			if(atom1 == connected14[ix4]){
				if(atom2 == connected14[ix4+3])
					return false;
			}
			else if(atom2 == connected14[ix4]){
				if(atom1 == connected14[ix4+3])
					return false;
			}
			ix4+=4;
		}
		return true;
	}

	// Adds an entry in the connected12 table for the bond
	//  between atom1 and atom2
	private void connect12(int atom1, int atom2){
		numberOf12Connections++;
		int numberOf12Connectionsx2 = numberOf12Connections * 2;

		if(numberOf12Connectionsx2 > connected12.length){
			// if we've run out of space, expand the connected12 array
			int largerConnectArray[] = new int[numberOf12Connectionsx2 + CONNECT_INCREMENT_SIZE*2];
			System.arraycopy(connected12, 0, largerConnectArray, 0, connected12.length);
			connected12 = largerConnectArray;
		}
		if (atom1>atom2){
			int temp = atom2;
			atom2 = atom1;
			atom1 = temp;
		}
		connected12[numberOf12Connectionsx2 - 2] = atom1;
		connected12[numberOf12Connectionsx2 - 1] = atom2;
	}

	// Adds an entry in the connected13 table for the connection
	//  between atom1, atom2, and atom3
	private void connect13(int atom1, int atom2, int atom3){
		if (atom1 > atom3){
			int temp = atom3;
			atom3 = atom1;
			atom1 = temp;
		}  
		// check for duplicate entries
		for(int i=0; i<numberOf13Connections; i++){
			int ix3 = i * 3;
			if (connected13[ix3] == atom1)
				if ((connected13[ix3 + 1] == atom2) &&
					(connected13[ix3 + 2] == atom3))
					return;
		}
		numberOf13Connections++;
		int numberOf13Connectionsx3 = numberOf13Connections * 3;

		if(numberOf13Connectionsx3 > connected13.length){
			// if we've run out of space, expand the connected13 array
			int largerConnectArray[] = new int[numberOf13Connectionsx3 + CONNECT_INCREMENT_SIZE*3];
			System.arraycopy(connected13, 0, largerConnectArray, 0, connected13.length);
			connected13 = largerConnectArray;
		}

		connected13[numberOf13Connectionsx3 - 3] = atom1;
		connected13[numberOf13Connectionsx3 - 2] = atom2;
		connected13[numberOf13Connectionsx3 - 1] = atom3;
	}

	// Adds an entry in the connected14 table for the connection
	//  between atom1, atom2, atom3, and atom4
	private void connect14(int atom1, int atom2, int atom3, int atom4){
		if (atom1 > atom4){
			int temp = atom4;
			atom4 = atom1;
			atom1 = temp;
			temp = atom3;
			atom3 = atom2;
			atom2 = temp;
		}
		// check for duplicate entries
		for(int i=0; i<numberOf14Connections; i++){
			int ix4 = i * 4;
			if (connected14[ix4] == atom1)
				if (connected14[ix4 + 1] == atom2)
					if ((connected14[ix4 + 2] == atom3) &&
						(connected14[ix4 + 3] == atom4))
						return;
		}
		numberOf14Connections++;
		int numberOf14Connectionsx4 = numberOf14Connections * 4;

		if(numberOf14Connectionsx4 > connected14.length){
			// if we've run out of space, expand the connected14 array
			int largerConnectArray[] = new int[numberOf14Connectionsx4 + CONNECT_INCREMENT_SIZE*4];
			System.arraycopy(connected14, 0, largerConnectArray, 0, connected14.length);
			connected14 = largerConnectArray;
		}

		connected14[numberOf14Connectionsx4 - 4] = atom1;
		connected14[numberOf14Connectionsx4 - 3] = atom2;
		connected14[numberOf14Connectionsx4 - 2] = atom3;
		connected14[numberOf14Connectionsx4 - 1] = atom4;
	}
	
	//Determines if the two atoms are 1-2 connected (directly)
	public boolean are12connected(int atom1, int atom2){
		
		for (int i=0; i<numberOf12Connections; i++){
			int ix2 = i*2;
			if (connected12[ix2]==atom1){
				if (connected12[ix2+1]==atom2)
					return true;
			}
			else if (connected12[ix2]==atom2){
				if (connected12[ix2+1]==atom1)
					return true;
			}
		}
		return false;
	}
	
	//Determines if the two atoms are 1-3 connected (1 atom between)
	public boolean are13connected(int atom1, int atom2){
		
		for (int i=0; i<numberOf13Connections; i++){
			int ix3 = i*3;
			if (connected13[ix3]==atom1){
				if (connected13[ix3+2]==atom2)
					return true;
			}
			else if (connected13[ix3]==atom2){
				if (connected13[ix3+2]==atom1)
					return true;
			}
		}
		return false;
	}
	
	//Determines if the two atoms are non-bonded: uses the nonBonded[] array;
	//NOTE: nonBonded[] is only valid if connectivityValid is true;
	//NOTE: nonBonded[] contains all pairs of atoms that are more than 3 bonds apart
	public boolean areNonBonded(int atom1, int atom2){
		
		int i=0;
		while ( i<(nonBonded.length-1) ){
			if ( (nonBonded[i]==atom1) && (nonBonded[i+1]==atom2) )
				return true;
			else if ( (nonBonded[i]==atom2) && (nonBonded[i+1]==atom1) )
				return true;
			
			i = i+2;
		}
		return false;
	}

	// This function returns the moleculeAtomNumber for the atom named atomName
	//  in the residue moleculeResidueNumber (molecule based numbering)
	public int getMoleculeAtomNumber(String atomName, int moleculeResidueNumber){
		for ( int i = 0; i < residue[ moleculeResidueNumber ].numberOfAtoms; i++ )
			if ( residue[ moleculeResidueNumber ].atom[ i ].name.equals( atomName ) )
				return( residue[ moleculeResidueNumber ].atom[ i ].moleculeAtomNumber );
		return -1;
	}
	
	//Return the coordinates from actualCoordinates[] for atom atomNum
	public float [] getActualCoord(int atomNum){
		float c[] = new float[3];
		c[0] = actualCoordinates[atomNum*3];
		c[1] = actualCoordinates[atomNum*3+1];
		c[2] = actualCoordinates[atomNum*3+2];
		return c;
	}
	
	
	//Determines if a bond exists between each pair of atoms in this molecule as a function of the distance (atom.coord[] and not
	//		actualCoordinates[]) and atom types;
	//Clears the bond information for all atoms in the molecule first;
	// TODO: for now determining bonds by distance is ok, but we eventually might want to read the CONECT terms if they exist
	public void determineBonds(){
		
		for ( int i = 0; i < numberOfAtoms; i++ ){ //reset the bond information
			atom[i].bond = null;
			atom[i].numberOfBonds = 0;
		}

		for ( int i = 0; i < numberOfAtoms; i++ ){
			for ( int j = i + 1 ; j < numberOfAtoms ; j++ ){
				double dist = atom[ i ].distance( atom[ j ] );
				if ((atom[i].elementType.equalsIgnoreCase("H") || atom[j].elementType.equalsIgnoreCase("H")) &&
						(atom[i].elementType.equalsIgnoreCase("S") || atom[j].elementType.equalsIgnoreCase("S"))) {
					if (dist < 1.5) {
						if(atom[i].strandNumber == atom[j].strandNumber){
							atom[i].addBond(j);
							atom[j].addBond(i);
						}
					}
				}
				else if ((atom[i].elementType.equalsIgnoreCase("H") || atom[j].elementType.equalsIgnoreCase("H")) &&
						!(atom[i].elementType.equalsIgnoreCase("S") || atom[j].elementType.equalsIgnoreCase("S"))) {
					if (dist < 1.2) {
						if(atom[i].strandNumber == atom[j].strandNumber){
							atom[i].addBond(j);
							atom[j].addBond(i);
						}
					}
				}
				else if (dist < 1.7) {
					if(atom[i].strandNumber == atom[j].strandNumber){
						atom[i].addBond(j);
						atom[j].addBond(i);
					}
				}
				else if (dist < 1.92) {
					if (atom[i].elementType.equalsIgnoreCase("S") ||
							 atom[j].elementType.equalsIgnoreCase("S")) {
						if(atom[i].strandNumber == atom[j].strandNumber){
							atom[i].addBond(j);
							atom[j].addBond(i);
						}
					}
					else if (atom[i].elementType.equalsIgnoreCase("Cl") ||
							 atom[j].elementType.equalsIgnoreCase("Cl")) {
						if(atom[i].strandNumber == atom[j].strandNumber){
							atom[i].addBond(j);
							atom[j].addBond(i);
						}
					}
				}
				else if (dist < 2.4) {
					if (atom[i].elementType.equalsIgnoreCase("S") &&
						 atom[j].elementType.equalsIgnoreCase("S")) {
						if(atom[i].strandNumber == atom[j].strandNumber){
							atom[i].addBond(j);
							atom[j].addBond(i);
						}
					}
				}
			}
		}
	}
	
	//Returns the molecule-relative residue number of the residue with PDB residue number pdbResNum
	public int mapPDBresNumToMolResNum(int pdbResNum){
		for (int i=0; i<numberOfResidues; i++){
			if (residue[i].getResNumber()==pdbResNum){
				return residue[i].moleculeResidueNumber;
			}
		}
		return -1;
	}

	////////////////////////////////////////////////////////////
	// Used by the SimpleMinimizer
	
	//Backs up the coordinates of the atoms in the atom[] array
	public void backupAtomCoord (){		
		backupCoordinates = new float [actualCoordinates.length];
		
		for (int i=0; i<numberOfAtoms; i++)
			backupCoord(atom[i]);
	}
	
	private void backupCoord(Atom updateAtom){
		int atomNumberx3 = updateAtom.moleculeAtomNumber * 3;
			
		backupCoordinates[atomNumberx3] = updateAtom.coord[0];	
		backupCoordinates[atomNumberx3 + 1] = updateAtom.coord[1];	
		backupCoordinates[atomNumberx3 + 2] = updateAtom.coord[2];	
	}
	
	//This function copies the coordinates of the actualCoordinates array
	//  into the coordinates for the atom[] array. All atoms are
	//  updated
	public void updateAtomCoord(){
		for(int i=0;i<numberOfAtoms;i++)
			updateCoord(atom[i]);
	}
	
	private void updateCoord(Atom updateAtom){
		int atomNumberx3 = updateAtom.moleculeAtomNumber * 3;
		
		updateAtom.coord[0] = actualCoordinates[atomNumberx3];	
		updateAtom.coord[1] = actualCoordinates[atomNumberx3 + 1];	
		updateAtom.coord[2] = actualCoordinates[atomNumberx3 + 2];
	}
	
	//Restore the backed-up coordinates for the atoms in the atom[] array
	public void restoreAtomCoord(){
		for(int i=0;i<numberOfAtoms;i++)
			restoreCoord(atom[i]);
	}
	
	private void restoreCoord(Atom updateAtom){				
		int atomNumberx3 = updateAtom.moleculeAtomNumber * 3;
		
		updateAtom.coord[0] = backupCoordinates[atomNumberx3];	
		updateAtom.coord[1] = backupCoordinates[atomNumberx3 + 1];	
		updateAtom.coord[2] = backupCoordinates[atomNumberx3 + 2];
	}
	////////////////////////////////////////////////////////////

	public float [] getActualCoords(){
		float retCoord[] = new float[actualCoordinates.length];
		System.arraycopy(actualCoordinates, 0, retCoord, 0, actualCoordinates.length);
		return retCoord;
	}
	
	public void setActualCoords(float ac[]){
		System.arraycopy(ac, 0, actualCoordinates, 0, actualCoordinates.length);
	}
}
