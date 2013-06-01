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
// AminoAcidTemplates.java
//
//  Version:           0.1
//
//
// authors:
//    initials    name            organization                email
//   ---------   --------------  ------------------------    ------------------------------
//     RHL        Ryan Lilien     Dartmouth College           ryan.lilien@dartmouth.edu
//
////////////////////////////////////////////////////////////////////////////////////////////

/** 
 * Writtten by Ryan Lilien (2001-2004)
 *
 */

import java.io.InputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.util.*;
import java.io.Serializable;

/**
 * 
 * This class reads from three data files specifying amino acid residues, N-terminal amino acid residues,
 * and C-terminal amino acid residues. Information read includes element type, atom type, limited connectivity,
 * and partial charge. By matching these amino acid templates to actual residues in a molecule, the corresponding
 * template atom types and partial charges can be assigned to the matched residues.
 *
 */
public class AminoAcidTemplates implements Serializable {

	String aaFilename = "all_amino94.in";
	String aaNTFilename = "all_aminont94.in";
	String aaCTFilename = "all_aminoct94.in";	

	Residue aaResidues[];   // arrays of residues
	Residue aaNTResidues[];
	Residue aaCTResidues[];
	int numAAs, numAANTs, numAACTs;

	// Constructor which reads amino acid information from the template files
	AminoAcidTemplates() throws Exception {

		// HANDLE THE NORMAL AAs	
		FileInputStream is = new FileInputStream( aaFilename );
		BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
		String curLine = null, tmpName = null;
		int dumPresent = 0;
		numAAs = numAANTs = numAACTs = 0;
		aaResidues = new Residue[30];
		aaNTResidues = new Residue[30];
		aaCTResidues = new Residue[30];

		// Skip over first 2 lines of header info
		curLine = bufread.readLine();
		curLine = bufread.readLine();

  		while( curLine != null ) {
			// Skip over first line which is the long amino acid name
			curLine = bufread.readLine();
			if (curLine.length() >= 4)
				if (curLine.substring(0,4).equalsIgnoreCase("stop")) {
					curLine = bufread.readLine();
					continue;
				}
			// Skip blank line
			curLine = bufread.readLine();
			// The next line contains the 3 letter amino acid name
			curLine = bufread.readLine();
			tmpName = getToken(curLine,1);
			Residue newRes = new Residue();
			newRes.name = tmpName;
			// Skip next 2 lines
			curLine = bufread.readLine();
			curLine = bufread.readLine();
			// Now we're into the section with atoms
			curLine = bufread.readLine();
			// Skip the dummy atoms
			dumPresent = 0;
			while (getToken(curLine,2).equalsIgnoreCase("DUMM")) {
				dumPresent++;
				curLine = bufread.readLine();
			}
			dumPresent++; // to adjust for 0-based
			while (!getToken(curLine,2).equals("")) {
				Atom at = new Atom();
				tmpName = getToken(curLine,2);
				at.changeType(tmpName);
				at.forceFieldType = getToken(curLine,3);
				at.charge = (float) (new Float(getToken(curLine,11)).floatValue());
				at.isBBatom = at.setIsBBatom();
				at.addBond(((new Integer(getToken(curLine,5))).intValue())-dumPresent);
				newRes.addAtom(at);  // add atom
				curLine = bufread.readLine();  // read next line
			}
			// Eventually we might want to be able to handle the improper
			//  torsions listed here

			// Add the residue to the array
			aaResidues[numAAs++] = newRes;
			// Read until the end of the residue
			boolean atDone = false;
			if (curLine.length() >= 4)
				atDone = curLine.substring(0,4).equalsIgnoreCase("done");
			else
				atDone = false;
			while (!atDone) {
				curLine = bufread.readLine();
				if (curLine.length() >= 4)
					atDone = curLine.substring(0,4).equalsIgnoreCase("done");
				}
		}
		bufread.close();


		// HANDLE THE N-terminal AAs
		is = new FileInputStream( aaNTFilename );
		bufread = new BufferedReader(new InputStreamReader(is));
		curLine = null;
		tmpName = null;
		
		// Skip over first 2 lines of header info
		curLine = bufread.readLine();
		curLine = bufread.readLine();

  		while( curLine != null ) {
			// Skip over first line which is the long amino acid name
			curLine = bufread.readLine();
			if (curLine.length() >= 4)
				if (curLine.substring(0,4).equalsIgnoreCase("stop")) {
					curLine = bufread.readLine();
					continue;
				}
			// Skip blank line
			curLine = bufread.readLine();
			// The next line contains the 3 letter amino acid name
			curLine = bufread.readLine();
			tmpName = getToken(curLine,1);
			Residue newRes = new Residue();
			newRes.name = tmpName;
			// Skip next 2 lines
			curLine = bufread.readLine();
			curLine = bufread.readLine();
			// Now we're into the section with atoms
			curLine = bufread.readLine();
			// Skip the dummy atoms
			dumPresent = 0;
			while (getToken(curLine,2).equalsIgnoreCase("DUMM")) {
				dumPresent++;
				curLine = bufread.readLine();
			}
			dumPresent++; // to adjust for 0-based
			while (!getToken(curLine,2).equals("")) {
				Atom at = new Atom();
				tmpName = getToken(curLine,2);
				at.name = tmpName;
				tmpName = tmpName.substring(0,1);		
				at.elementType = tmpName;
				at.forceFieldType = getToken(curLine,3);
				at.charge = (float) (new Float(getToken(curLine,11)).floatValue());
				at.isBBatom = at.setIsBBatom();
				at.addBond(((new Integer(getToken(curLine,5))).intValue())-dumPresent);
				newRes.addAtom(at);  // add atom
				curLine = bufread.readLine();  // read next line
			}
			// Eventually we might want to be able to handle the improper
			//  torsions listed here

			// Add the residue to the array
			aaNTResidues[numAANTs++] = newRes;
			// Read until the end of the residue
			boolean atDone = false;
			if (curLine.length() >= 4)
				atDone = curLine.substring(0,4).equalsIgnoreCase("done");
			else
				atDone = false;
			while (!atDone) {
				curLine = bufread.readLine();
				if (curLine.length() >= 4)
					atDone = curLine.substring(0,4).equalsIgnoreCase("done");
				}
		}
		bufread.close();

		// HANDLE THE C-terminal AAs
		is = new FileInputStream( aaCTFilename );
		bufread = new BufferedReader(new InputStreamReader(is));
		curLine = null;
		tmpName = null;
		
		// Skip over first 2 lines of header info
		curLine = bufread.readLine();
		curLine = bufread.readLine();

  		while( curLine != null ) {
			// Skip over first line which is the long amino acid name
			curLine = bufread.readLine();
			if (curLine.length() >= 4)
				if (curLine.substring(0,4).equalsIgnoreCase("stop")) {
					curLine = bufread.readLine();
					continue;
				}
			// Skip blank line
			curLine = bufread.readLine();
			// The next line contains the 3 letter amino acid name
			curLine = bufread.readLine();
			tmpName = getToken(curLine,1);
			Residue newRes = new Residue();
			newRes.name = tmpName;
			// Skip next 2 lines
			curLine = bufread.readLine();
			curLine = bufread.readLine();
			// Now we're into the section with atoms
			curLine = bufread.readLine();
			// Skip the dummy atoms
			dumPresent = 0;
			while (getToken(curLine,2).equalsIgnoreCase("DUMM")) {
				dumPresent++;
				curLine = bufread.readLine();
			}
			dumPresent++; // to adjust for 0-based
			while (!getToken(curLine,2).equals("")) {
				Atom at = new Atom();
				tmpName = getToken(curLine,2);
				at.name = tmpName;
				tmpName = tmpName.substring(0,1);		
				at.elementType = tmpName;
				at.forceFieldType = getToken(curLine,3);
				at.charge = (float) (new Float(getToken(curLine,11)).floatValue());
				at.isBBatom = at.setIsBBatom();
				at.addBond(((new Integer(getToken(curLine,5))).intValue())-dumPresent);
				newRes.addAtom(at);  // add atom
				curLine = bufread.readLine();  // read next line
			}
			// Eventually we might want to be able to handle the improper
			//  torsions listed here

			// Add the residue to the array
			aaCTResidues[numAACTs++] = newRes;
			// Read until the end of the residue
			boolean atDone = false;
			if (curLine.length() >= 4)
				atDone = curLine.substring(0,4).equalsIgnoreCase("done");
			else
				atDone = false;
			while (!atDone) {
				curLine = bufread.readLine();
				if (curLine.length() >= 4)
					atDone = curLine.substring(0,4).equalsIgnoreCase("done");
				}
		}
		bufread.close();
		
}

	/******************************/
	// This function returns the xth token in string s
	private String getToken(String s, int x) {
	
		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");
		
		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
			  st.nextToken();
			else {
				return(new String(""));
			}
		}

		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken

}
