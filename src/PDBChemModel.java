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
// PDBChemModel.java
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
 * Written by Ryan Lilien (2000-2004); some modifications by Ivelin Georgiev (2004-2009)
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

import java.io.InputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;

/**
 * Reads a pdb file and creates a molecule object.
 */
class PDBChemModel {

  // Reads a pdb file and creates a molecule object
	PDBChemModel (Molecule m, InputStream is) throws Exception {
		
		int modelAtomNumber = 0;
		int residueNumber = 0;
		int strandNumber = 0;
		String atomName = "ZZZ", residueName = "ZZZ", strandName = "0"; 
		String lastResidueName = "ZZZ", fullResidueName = "ZZZZZZZZZ";
		String elementType = "  ";
		String curLine = null;
		String tmpStg = null;
		int tmpInt;
		char[] tmpChr = new char[15];
		float	x = 0f, y = 0f, z = 0f;
		Atom	newAtom;
		boolean newStrandPending = true;  // Will add the first strand with the first atom
		Residue newResidue = null;

		BufferedReader bufread = new BufferedReader(new InputStreamReader(is));

		curLine = bufread.readLine();
		if (curLine == null) { // stop if we've reached EOF
			bufread.close();
		}

	  	while(curLine != null) {
	  		// First pad line to 80 characters
	  		tmpInt = curLine.length();
	  		for (int i=0; i < (80-tmpInt); i++)
	  			curLine += " ";
	
	   		if ((curLine.regionMatches(true,0,"ATOM  ",0,6)) || (curLine.regionMatches(true,0,"HETATM",0,6))) {
					
	   			// Is an ATOM line
				tmpStg = curLine.substring(6,11);  // Snag atom serial number
				tmpStg = tmpStg.trim();
				modelAtomNumber = (new Integer(tmpStg)).intValue();
	   			atomName = curLine.substring(12,16);  // Snag atom name
	   			atomName = atomName.trim();
	   			residueName = curLine.substring(17,20);  // Snag short residue name
	   			residueName = residueName.trim();
	   			fullResidueName = curLine.substring(17,26);  // Snag full residue atom name
	   			fullResidueName = fullResidueName.trim();
	
				if (!(fullResidueName.equals(lastResidueName))) {
					if (newResidue != null) {
						if (newStrandPending)
							m.addResidue(strandNumber-1, newResidue);
						else
							m.addResidue(strandNumber, newResidue);
					}
					newResidue = new Residue();
					newResidue.name = residueName;
					newResidue.fullName = fullResidueName;
					lastResidueName = fullResidueName;
					residueNumber++;
				}
	
				if (newStrandPending) {
					strandName = curLine.substring(21,22);  // Snag strand name
					m.addStrand(strandName);
				 	newStrandPending = false;
				}
				tmpStg = curLine.substring(30,38);  // Snag x coord
				x = (float) new Double(tmpStg).doubleValue();
				tmpStg = curLine.substring(38,46);  // Snag y coord
				y = (float) new Double(tmpStg).doubleValue();
				tmpStg = curLine.substring(46,54);  // Snag z coord
				z = (float) new Double(tmpStg).doubleValue();
	
				elementType = curLine.substring(76,78);  // Snag atom elementType
				elementType = elementType.trim();
				// If we can't get element type from substring(76,78) snag
				//  the first character of the atom name
				if (elementType.equalsIgnoreCase(""))
					elementType = getEleType(curLine.substring(12,15));
				newAtom = new Atom(atomName,x,y,z);
				newAtom.modelAtomNumber = modelAtomNumber;
				newAtom.strandNumber = strandNumber;
				newAtom.elementType = elementType;
				newResidue.addAtom(newAtom);
			} // end ATOM line
			else if (curLine.regionMatches(true,0,"TER   ",0,6)) {
				// Is the end of a strand
				lastResidueName = "ZZZ";
				residueNumber = 0;
				strandNumber++;
				newStrandPending = true;
			} // end TER line
			else   // is a line we skip
				;
			curLine = bufread.readLine();  // attempt to read next line
		}  // end while (curLine != null)
	
		if (newStrandPending)
			m.addResidue(strandNumber-1,newResidue);
		else
			m.addResidue(strandNumber,newResidue);
		bufread.close();  // close the buffer
		
		//Determine the bonds between the atoms in the molecule
		m.determineBonds();
		
		// Assign the molecule relative atom numbers
		m.updateMoleculeAtomNumbers();
	}
	
	// This function pulls the element type from
	//  the atom name
	private String getEleType(String str){
		
		int start=0, end=-1;
		int i=0;
		while( (str.charAt(i)==' ') || ((str.charAt(i)>='0') && (str.charAt(i)<='9')) ) {
			i++;
		}
		start = i;
		end = i++;
		if (i<str.length())
			if((str.charAt(i)>='a') && (str.charAt(i)<='z'))
				end = i;
		return(str.substring(start,end+1));	
	}
}
