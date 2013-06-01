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
//	EEF1.java
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

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.Serializable;
import java.util.StringTokenizer;

/**
 * Manages the EEF1 solvation parameters;
 * 		EEF1 handles only natural amino acids, so solvation energies are computed for proteins and ligands (natural amino acids) only, and not cofactors;
 * 		Some important notes about assumptions of the EEF1 model are given in the comments to getSolvationParameters()
 * 
*/
public class EEF1 implements Serializable {
	
	Molecule m = null;
	
	//Variables to store the EEF1 solvation paramaters;
	//These values are specific to eef1parm.dat, and may have to be modified for different
	//		versions of the parameter file
	final int numAtTypesEEF1 = 17;
	String groupEEF1names[] = new String[numAtTypesEEF1];
	String atTypeEEF1names[] = new String[numAtTypesEEF1];
	float atEEF1Vol[] = new float[numAtTypesEEF1];
	float dGiRef[] = new float[numAtTypesEEF1];
	float dGiFree[] = new float[numAtTypesEEF1];
	float dHiRef[] = new float[numAtTypesEEF1];
	float dCpiRef[] = new float[numAtTypesEEF1];
	float lambdai[] = new float[numAtTypesEEF1];
	float vdWri[] = new float[numAtTypesEEF1];	

	
	//constructor
	EEF1(Molecule mol){
		m = mol;
	}
	
	//Reads the solvation parameters for EEF1 from the file eef1parm.dat;
	//If a different versio of the paramater file is used, some changes to the
	//		current function may be necessary
	public void readEEF1parm() throws Exception {
		
		FileInputStream is = new FileInputStream( "eef1parm.dat" );
		BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
		String curLine = null;
		int tmpInt = 0;

		//Skip over the first line of header info
		curLine = bufread.readLine();
		
		//Read each parameter line, which is in the format:
		//	Group type   Atom type   Volume   DGi(ref)   DGi(free)   DHi(ref)   DCpi(ref)   lambda_i   vdWri
		curLine = bufread.readLine();		
		tmpInt = 0; // temporary integer
		// Until we're at a blank line (or until we've read numAtomTypes)
		while (!(curLine==null)) {
			groupEEF1names[tmpInt] = getToken(curLine,1);  // snag group name
			atTypeEEF1names[tmpInt] = getToken(curLine,2);  // snag atom type
			atEEF1Vol[tmpInt] = (new Float(getToken(curLine,3))).floatValue();
			dGiRef[tmpInt] = (new Float(getToken(curLine,4))).floatValue();
			dGiFree[tmpInt] = (new Float(getToken(curLine,5))).floatValue();
			dHiRef[tmpInt] = (new Float(getToken(curLine,6))).floatValue();
			dCpiRef[tmpInt] = (new Float(getToken(curLine,7))).floatValue();
			lambdai[tmpInt] = (new Float(getToken(curLine,8))).floatValue();
			vdWri[tmpInt] = (new Float(getToken(curLine,9))).floatValue();
			tmpInt++;
			curLine = bufread.readLine();
		}
		
		bufread.close();
	}
	
	//Returns the dG(ref), dG(free), atom volume, and lambda parameters for 
	//		the atom specified by moleculeAtomNumber atom1;
	//The mapping is based on the group types in eef1parm.dat: different versions
	//		of the parameter file may require modifications in this function;
	//The assumption is made that the temperature in Kelvin is 298.15, since this is
	//		the temperature for which dGi(free) and dGi(ref) are derived in eef1parm.dat
	//Note that some of the C and S groups are *extended* heavy atoms (heavy atoms which include
	//		the bonded hydrogen atoms), so in the original EEF1 and CHARMM implementation, the
	//		corresponding vdW radii are increased accordingly. We thus use the vdW radii from the
	//		extended Param 19 (see Appendix B, Neria & Karplus, 1996) parameters for the solvation energy computation only;
	//The groups in eef1 are based on the natural amino acids. Additional groups may have to be included, 
	//		and the current groups may have to be modified if new non-natural amino acids or other
	//		atom types or molecules (such as AMP or ATP) are included 
	public boolean getSolvationParameters(int atom1,double dGref[],double dGfree[],
			double atVolume[],double lambda[], double vdWr[]){
		
		Atom at1 = m.atom[atom1];
		
		int solvGroupIndex = getSolvGroupIndex(at1);
		
		//System.out.print(m.residue[at1.moleculeResidueNumber].fullName+" "+at1.name);
		
		if (solvGroupIndex==-1) //group not found
			return false;
		else {
			dGref[0] = dGiRef[solvGroupIndex];
			dGfree[0] = dGiFree[solvGroupIndex];
			atVolume[0] = atEEF1Vol[solvGroupIndex];
			lambda[0] = lambdai[solvGroupIndex];
			vdWr[0] = vdWri[solvGroupIndex];
			
			//System.out.println("           "+groupEEF1names[solvGroupIndex]);
			
			return true;
		}		
	}
	
	//Determines the solvation group index for the given atom;
	//This function assumes that the predfefined groups are as given in eef1parm.dat;
	//		Modifcations to this functions will be required if different group types are used
	private int getSolvGroupIndex(Atom at1){
		
		String elementType = at1.elementType;
		String AAname = m.residue[at1.moleculeResidueNumber].name; //the AA name to which this atom belongs
		boolean aromatic = isAromatic(at1);
		int numBoundH = getNumBoundH(at1);
		
		if (elementType.equalsIgnoreCase("C")) {
			
			if ( at1.name.equalsIgnoreCase("CG") && 
					( AAname.equalsIgnoreCase("TYR") || AAname.equalsIgnoreCase("PHE") || AAname.equalsIgnoreCase("HIS") ) )
				return findSolvGroup("CR"); //CG of (TYR or PHE or HIS)
			
			else if ( ( at1.name.equalsIgnoreCase("CG") || at1.name.equalsIgnoreCase("CD2") || at1.name.equalsIgnoreCase("CE2") )
					&& AAname.equalsIgnoreCase("TRP") )
				return findSolvGroup("CR"); // (CG or CD2 or CE2) of TRP
			
			else if ( at1.name.equalsIgnoreCase("CZ") && AAname.equalsIgnoreCase("ARG") )
				return findSolvGroup("CR"); //CZ of PHE
			
			else if ( (aromatic) && (numBoundH==1) ) // extended aromatic C with 1 H, so CR1E group
				return findSolvGroup("CR1E");
			
			else if ( (!aromatic) && (numBoundH==1) ) // extended aliphatic C with 1 H, so CH1E group
				return findSolvGroup("CH1E");
			
			else if ( (!aromatic) && (numBoundH==2) ) // extended aliphatic C with 2 H, so CH2E group
				return findSolvGroup("CH2E");
			
			else if ( (!aromatic) && (numBoundH==3) ) // extended aliphatic C with 3 H, so CH3E group
				return findSolvGroup("CH3E");
			
			else // default C group
				return findSolvGroup("C");
		}
		else if (elementType.equalsIgnoreCase("N")) {
			
			if (AAname.equalsIgnoreCase("PRO")) // backbone N of PRO, so N group
				return findSolvGroup("N");
			
			else if (numBoundH==1) // any N with 1 H
				return findSolvGroup("NH1");
			
			else if (AAname.equalsIgnoreCase("ARG") && (numBoundH==2)) // guanidinium N in (side-chain of) ARG, so NC2 group
				return findSolvGroup("NC2");
			
			else if ( (aromatic) && (numBoundH==0) ) // aromatic (side-chain) N with 0 H, so NR group
				return findSolvGroup("NR");
			
			else if (numBoundH==2) // any (side-chain) N with 2 H, so NH2 group
				return findSolvGroup("NH2");
			
			else if (numBoundH==3) // any N with 3 H, so NH3 group
				return findSolvGroup("NH3");
			
			else //unknown nitrogen group
				return -1;
		}
		else if (elementType.equalsIgnoreCase("O")) {
			
			if (numBoundH==1) // hydroxyl oxygen, so OH1 group
				return findSolvGroup("OH1");
			
			else if (isCarboxylO(at1)) // carboxyl oxygen, so OC group
				return findSolvGroup("OC");
			
			else //all other oxygens (carbonyl oxygen)
				return findSolvGroup("O");
		}
		else if (elementType.equalsIgnoreCase("S")) {
			
			if (numBoundH==1) // extended S with 1 H, so SH1E group
				return findSolvGroup("SH1E");
			
			else // any other S, so S group
				return findSolvGroup("S");
		}
		else //group not found
			return -1;
	}
	
	//Find the solvation group index for the group name grName
	private int findSolvGroup(String grName){
		for (int i=0; i<groupEEF1names.length; i++){
			if (groupEEF1names[i].equalsIgnoreCase(grName))
				return i;
		}
		return -1; //group name not found
	}
	
	//Determines if ths given atom is a carboxyl oxygen (there should be only one carbon bound to the
	//		given oxygen and exactly two oxygens (including the given one) should be bound to that carbon)
	private boolean isCarboxylO(Atom at1){
		
		if (!at1.elementType.equalsIgnoreCase("O")) //not O
			return false;
		else {
			
			int numBoundC = 0;
			boolean isCarboxylO = false;
			if (at1.bond!=null){ //check bonds
				for (int i=0; i<at1.bond.length; i++){		
					
					Atom at2 = m.atom[at1.bond[i]];
					
					if (at2.elementType.equalsIgnoreCase("C")){ //found a bound C
						
						numBoundC++;
						if (numBoundC>1) //more than 1 C bound
							return false;
						else
							isCarboxylO = isCOhelper(at2);
					}
				}
			}
			return isCarboxylO;			
		}
	}
	
	//Used by isCarboxylO(); see comments above
	private boolean isCOhelper(Atom at2){		
		
		int numBoundO = 0;
		if (at2.bond!=null){ //check bonds
			for (int j=0; j<at2.bond.length; j++){
				if (m.atom[at2.bond[j]].elementType.equalsIgnoreCase("O"))
					numBoundO++;
			}
		}
		
		if (numBoundO==2)
			return true;
		else
			return false;
	}
	
	//Determines if the given heavy atom is aromatic
	private boolean isAromatic(Atom at1){
		
		String AAname = m.residue[at1.moleculeResidueNumber].name; //the AA name of the residue to which this atom belongs
		
		boolean isHis = (AAname.equalsIgnoreCase("HIS") || AAname.equalsIgnoreCase("HIP") || AAname.equalsIgnoreCase("HID") || AAname.equalsIgnoreCase("HIE"));
		
		if ( !(AAname.equalsIgnoreCase("PHE") || AAname.equalsIgnoreCase("TYR") 
					|| AAname.equalsIgnoreCase("TRP") || isHis) ) //not PHE and not TYR and not TRP or HIS, so not aromatic
				return false;
		else {//PHE or TYR or TRP or HIS, so check if aromatic atom
			
			if (!at1.elementType.equalsIgnoreCase("C")){ //not a carbon
				if ( at1.elementType.equalsIgnoreCase("N") && 
						( !at1.name.equalsIgnoreCase("N") && (AAname.equalsIgnoreCase("TRP") || isHis) ) ) //N in the ring of TRP or HIS
					return true;
				else
					return false;
			}
			else {//a carbon, so check if CA or CB or C
				String atomName = at1.name;
				if (atomName.equalsIgnoreCase("CA") || atomName.equalsIgnoreCase("CB") 
						|| atomName.equalsIgnoreCase("C")) //CA or CB or C, so cannot be aromatic
					return false;
				else //aromatic
					return true;
			}
		}
	}
	
	//Get the number of H bound to the given atom
	private int getNumBoundH(Atom at1){
		
		int numBoundH = 0;
		for (int i=0; i<at1.bond.length; i++){
			if (m.atom[at1.bond[i]].elementType.equalsIgnoreCase("H"))
				numBoundH++;
		}
		
		return numBoundH;
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
