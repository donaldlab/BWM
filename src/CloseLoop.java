/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.1 beta
	Copyright (C) 2001-2012 Bruce Donald Lab, Duke University

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

	<signature of Bruce Donald>, Mar 1, 2012
	Bruce Donald, Professor of Computer Science
*/

///////////////////////////////////////////////////////////////////////////////////////////////
//	CloseLoop.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  JDJ        Jonathan Jou    	  Duke University             jj@cs.duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/* This perturbation involves finding alternate ways to close loops
 * A plain "loop closure adjustment" does this by finding alternate sets of dihedrals
 * for a three-residue segment without altering the chain elsewhere
 * This is a discrete type of perturbation

 * There are also some specialized versions:
 * 1. A helix or sheet reduction is a loop closure adjustment containing one or more non-loop residues
 * (which will be removed from their helix or sheet)
 * 2. A helix or sheet expansion is formed by adding a loop residue to an adjacent helix or sheet
 * and then trying to close the loop
 * 3. A partial structure switch uses experimental data (maybe its own class??)
 *
 */

public class CloseLoop
{

    
    int solutions = 0;
	int[][] atomList;
	double[][][] solutionDihedrals; /*[solution][residue][phi/psi]*/
	double [][] initialDihedrals = new double[3][2];
	Molecule m;
	public KSParser debugWriter = null;
	boolean CToNClosure = false;
	
	int[] resAffected;

    
    public CloseLoop(String t, Molecule molec, int[] reslist, int[] NAnchorAtoms, int[] CAnchorAtoms){
		m = molec;
		resAffected = reslist;
		atomList = new int[3][];
		atomList[0] = NAnchorAtoms;
		atomList[2] = CAnchorAtoms;
		if(NAnchorAtoms != null && CAnchorAtoms == null)
			CToNClosure = true;
		for(int resNumber = 0; resNumber < resAffected.length; resNumber++)
		{
			int[] phiAtoms = m.getResiduePhiDihedralAtoms(resAffected[resNumber]);
			if(phiAtoms[0] < 0 || phiAtoms[1] < 0 || phiAtoms[2] < 0 || phiAtoms[3] < 0) return;
			int[] psiAtoms = m.getResiduePsiDihedralAtoms(resAffected[resNumber]);

			if(psiAtoms[0] < 0 || psiAtoms[1] < 0 || psiAtoms[2] < 0 || psiAtoms[3] < 0) return;
			initialDihedrals[resNumber][0] = m.getTorsion(phiAtoms[0], phiAtoms[1], phiAtoms[2], phiAtoms[3]);
			initialDihedrals[resNumber][1] = m.getTorsion(psiAtoms[0], psiAtoms[1], psiAtoms[2], psiAtoms[3]);
			
		}
	}


    public CloseLoop(String t, Molecule molec, int resList[]){
        m=molec;
        resAffected=resList;
    }

	private void populateAtomList()
	{
		if(atomList == null)
			atomList = new int[3][];
		for(int i = 0; i < resAffected.length; i++)
		{
			if(atomList[i] != null)
				continue;
			atomList[i] = new int[3];
			atomList[i][0] = m.residue[resAffected[i]].getAtomNameToMolnum("N");
			atomList[i][1] = m.residue[resAffected[i]].getAtomNameToMolnum("CA");
			atomList[i][2] = m.residue[resAffected[i]].getAtomNameToMolnum("C");
		}
	}
	
	

	private double[][][] storeDihedrals(int numSolutions, float r_soln_n[][][], float r_soln_a[][][], float r_soln_c[][][])
	{
		double[][][] dihedrals = new double[numSolutions][][];
		int[][] atomList2 = new int[3][3];
		float[][] backupCoordinates = new float[9][3];
		for(int i = 0; i < resAffected.length; i++)
		{
			atomList2[i][0] = m.residue[resAffected[i]].getAtomNameToMolnum("N");
			atomList2[i][1] = m.residue[resAffected[i]].getAtomNameToMolnum("CA");
			atomList2[i][2] = m.residue[resAffected[i]].getAtomNameToMolnum("C");
		}
		for(int s = 0; s < numSolutions; s++)
		{
			System.out.println("Begin Solution "+s);
	        for(int i = 0; i < 3; i++)
	        {
	        	for(int j = 0; j < 3; j++)
	        	{
	        		float[] solutionArray = r_soln_n[s][j];
	        		switch(i)
	        		{
		        		case 0:
		        			solutionArray = r_soln_n[s][j];
		        			break;
		        		case 1:
		        			solutionArray = r_soln_a[s][j];
		        			break;
		        		case 2:
		        			solutionArray = r_soln_c[s][j];
		        			break;
	        		}
	        		for(int k = 0; k < 3; k++)
	        		{
	        			backupCoordinates[i*3+j][k] = m.actualCoordinates[atomList2[j][i]*3+k];
	        			m.actualCoordinates[atomList2[j][i]*3+k] = solutionArray[k];
	        		}
	        	}
    			m.resolveCoordinates();

	        }
	        if(debugWriter != null)
	        	debugWriter.saveMolecule(m, "output-2-13-"+s+"-backbone.pdb", 0);
	        
	        
	        dihedrals[s] = getDihedrals();

	        
	        /* Restore coordinates */
	        for(int i = 0; i < 3; i++)
	        {
	        	for(int j = 0; j < 3; j++)
	        	{
	        		for(int k = 0; k < 3; k++) 
	        		{
	        			m.actualCoordinates[atomList2[j][i]*3 + k] = backupCoordinates[i*3 + j][k];
	        		}
	        	}
	        }
			m.resolveCoordinates();
	        if(dihedrals[s] == null) 
	        	return null;
		}

        return dihedrals;
		
	}
	
	private double[][] getDihedrals()
	{
		double[][] dihedrals = new double[3][2];
		for(int resNumber = 0; resNumber < resAffected.length; resNumber++)
		{
			int[] phiAtoms = m.getResiduePhiDihedralAtoms(resAffected[resNumber]);
			if(phiAtoms[0] < 0 || phiAtoms[1] < 0 || phiAtoms[2] < 0 || phiAtoms[3] < 0) return null;
			int[] psiAtoms = m.getResiduePsiDihedralAtoms(resAffected[resNumber]);

			if(psiAtoms[0] < 0 || psiAtoms[1] < 0 || psiAtoms[2] < 0 || psiAtoms[3] < 0) return null;
			dihedrals[resNumber][0] = m.getTorsion(phiAtoms[0], phiAtoms[1], phiAtoms[2], phiAtoms[3]);
			//System.out.println("Phi: "+dihedrals[resNumber][0]);
			dihedrals[resNumber][1] = m.getTorsion(psiAtoms[0], psiAtoms[1], psiAtoms[2], psiAtoms[3]);
			//System.out.println("Psi atoms: "+psiAtoms[0]+", "+psiAtoms[1]+", "+psiAtoms[2]+", "+psiAtoms[3]);
			//System.out.println("Psi: "+dihedrals[resNumber][1]+", from "+ m.getTorsion(psiAtoms[0], psiAtoms[1], psiAtoms[2], psiAtoms[3]));
			
			
		}
		return dihedrals;
	}

	public void translateBackBone(int s, float r_soln_n[][][], float r_soln_a[][][], float r_soln_c[][][])
	{
		int[][] atomList2 = new int[3][3];
		float[][] backupCoordinates = new float[9][3];
		for(int i = 0; i < resAffected.length; i++)
		{
			atomList2[i][0] = m.residue[resAffected[i]].getAtomNameToMolnum("N");
			atomList2[i][1] = m.residue[resAffected[i]].getAtomNameToMolnum("CA");
			atomList2[i][2] = m.residue[resAffected[i]].getAtomNameToMolnum("C");
		}

	        for(int i = 0; i < 3; i++)
	        {
	        	for(int j = 0; j < 3; j++)
	        	{
	        		float[] solutionArray = r_soln_n[s][j];
	        		switch(i)
	        		{
		        		case 0:
		        			solutionArray = r_soln_n[s][j];
		        			break;
		        		case 1:
		        			solutionArray = r_soln_a[s][j];
		        			break;
		        		case 2:
		        			solutionArray = r_soln_c[s][j];
		        			break;
	        		}
	        		for(int k = 0; k < 3; k++)
	        		{
	        			backupCoordinates[i*3+j][k] = m.actualCoordinates[atomList2[j][i]*3+k];
	        			m.actualCoordinates[atomList2[j][i]*3+k] = solutionArray[k];
	        		}
	        	}
    			m.resolveCoordinates();

	        }
	}
	
	
	public void restore()
	{
		applyDihedrals(initialDihedrals, CToNClosure);
	}
	
	private void applyDihedrals(double[][] dihedrals, boolean CToN)
	{
		int allAtoms[] = concatenateArrays(m.residue[resAffected[0]].getAtomList(false, true, true, true),
				m.residue[resAffected[1]].getAtomList(true,  true,  true, true), m.residue[resAffected[2]].getAtomList(true,  true,  true,  true));
		if(CToN)
			applyReverseDihedral(dihedrals, CToN, allAtoms);
		else
			applyForwardDihedral(dihedrals, CToN, allAtoms);
	}


	private void applyForwardDihedral(double[][] dihedrals, boolean CToN,
			int[] allAtoms) {
		for(int resNumber = 0; resNumber < 3; resNumber++)
		{
			/* TODO: Fetch atom list */
			int[] phiAtoms = m.getResiduePhiDihedralAtoms(resAffected[resNumber]);
			int[] psiAtoms = m.getResiduePsiDihedralAtoms(resAffected[resNumber]);
			for(int angleIndex = 0; angleIndex < 2; angleIndex++)
			{
				int[] atoms = phiAtoms;
				if(angleIndex>0) atoms = psiAtoms;
				int[] remainingAtoms = allAtoms;
				remainingAtoms = getAtoms(resNumber, angleIndex, CToN);

				double lastPsi= m.getTorsion(atoms[0], atoms[1], atoms[2], atoms[3]);
				
				//System.out.println((resNumber*2 + angleIndex)+": setting to "+dihedrals[resNumber][angleIndex]+" from "+lastPsi);
				m.setTorsion(atoms[0], atoms[1], atoms[2], atoms[3], dihedrals[resNumber][angleIndex], remainingAtoms, remainingAtoms.length);
			}
		}
	}
	
	private void applyReverseDihedral(double[][] dihedrals, boolean CToN,
			int[] allAtoms) {
		for(int resNumber = 2; resNumber >= 0; resNumber--)
		{
			/* TODO: Fetch atom list */
			int[] phiAtoms = m.getResiduePhiDihedralAtoms(resAffected[resNumber]);
			int[] psiAtoms = m.getResiduePsiDihedralAtoms(resAffected[resNumber]);
			for(int angleIndex = 1; angleIndex >= 0; angleIndex--)
			{
				int[] atoms = phiAtoms;
				if(angleIndex>0) atoms = psiAtoms;
				int[] remainingAtoms = allAtoms;
				remainingAtoms = getAtoms(resNumber, angleIndex, CToN);

				double lastPsi= m.getTorsion(atoms[0], atoms[1], atoms[2], atoms[3]);
				
				//System.out.println((resNumber*2 + angleIndex)+": setting to "+dihedrals[resNumber][angleIndex]+" from "+lastPsi);
				m.setTorsion(atoms[3], atoms[2], atoms[1], atoms[0], dihedrals[resNumber][angleIndex], remainingAtoms, remainingAtoms.length);
				if(debugWriter != null)
					debugWriter.saveMolecule(m, "pSolution-"+(resNumber*2+angleIndex)+"-"+".pdb", 0);
				
			}
		}
	}
	
	
	private int[] getAtoms(int resNumber, int angleIndex, boolean CToN)
	{
		int[] remainingAtoms = null;
		if(CToN)
		{
			switch(resNumber*2 + angleIndex)
			{
			case 5:
				remainingAtoms = concatenateArrays(m.residue[resAffected[0]].getAtomList(true, true, true, true),
						m.residue[resAffected[1]].getAtomList(true,  true,  true, true), m.residue[resAffected[2]].getAtomList(true,  true,  true,  false));
				break;
			case 4:
				remainingAtoms = concatenateArrays(m.residue[resAffected[0]].getAtomList(true, true, true, true),
						m.residue[resAffected[1]].getAtomList(true,  true,  true, true), m.residue[resAffected[2]].getAtomList(true,  false,  false,  false));
				break;
			case 3:
				remainingAtoms = concatenateArrays(m.residue[resAffected[0]].getAtomList(true, true, true, true),
						m.residue[resAffected[1]].getAtomList(true,  true,  true, false));
				break;
			case 2:
				remainingAtoms = concatenateArrays(m.residue[resAffected[0]].getAtomList(true, true, true, true),
						m.residue[resAffected[1]].getAtomList(true,  false,  false, false));
				break;
			case 1:
				remainingAtoms = concatenateArrays(m.residue[resAffected[0]].getAtomList(true, true, true, false));
				break;
			case 0:
				remainingAtoms = concatenateArrays(m.residue[resAffected[0]].getAtomList(true, false, false, false));
				break;
			}
		}
		else
			switch(resNumber*2 + angleIndex)
			{
				case 0:
					remainingAtoms = concatenateArrays(m.residue[resAffected[0]].getAtomList(false, true, true, true),
							m.residue[resAffected[1]].getAtomList(true,  true,  true, true), m.residue[resAffected[2]].getAtomList(true,  true,  true,  true));
					break;
				case 1:
					remainingAtoms = concatenateArrays(m.residue[resAffected[0]].getAtomList(false, false, false, true),
							m.residue[resAffected[1]].getAtomList(true,  true,  true, true), m.residue[resAffected[2]].getAtomList(true,  true,  true,  true));
					break;
				case 2:
					remainingAtoms = concatenateArrays(
							m.residue[resAffected[1]].getAtomList(false,  true,  true, true), m.residue[resAffected[2]].getAtomList(true,  true,  true,  true));
					break;
				case 3:
					remainingAtoms = concatenateArrays(
							m.residue[resAffected[1]].getAtomList(false,  false,  false, true), m.residue[resAffected[2]].getAtomList(true,  true,  true,  true));
					break;
				case 4:
					remainingAtoms = concatenateArrays(
							m.residue[resAffected[2]].getAtomList(false,  true,  true,  true));
					break;
				case 5:
					remainingAtoms = concatenateArrays(
							m.residue[resAffected[2]].getAtomList(false,  false,  false,  true));
					break;
			}
			return remainingAtoms;
	}
	
	public void applySolution(int solutionNumber)
	{
		applyDihedrals(solutionDihedrals[solutionNumber], CToNClosure);
		
	}

    public int CalculateSolutions(){//Solve the loop closure equations given the current state of predecessor perturbations
        //Store conformational change info and return the number of solutions

        //Get necessary lengths, angles, dihedrals

        int tripepRes[] = null;//Residues in the tripeptide to be closed

        tripepRes = resAffected;


        TripeptideClosure tc = new TripeptideClosure(m, tripepRes);

        float r_soln_n[][][] = new float[16][3][3];
        float r_soln_a[][][] = new float[16][3][3];
        float r_soln_c[][][] = new float[16][3][3];

		populateAtomList();

        float NAnchorN[] = m.getActualCoord( atomList[0][0] );//Coordinates of the first residue's N
        float NAnchorCA[] = m.getActualCoord( atomList[0][1] );//Its CA
        float NAnchorC[] = m.getActualCoord( atomList[0][2] );

        float midN[] = m.getActualCoord( atomList[1][0] );//Starting coordinates of the middle N
        float midCA[] = m.getActualCoord( atomList[1][1] );
        float midC[] = m.getActualCoord( atomList[1][2] );

        float CAnchorN[] = m.getActualCoord( atomList[2][0] );
        float CAnchorCA[] = m.getActualCoord( atomList[2][1] );
        float CAnchorC[] = m.getActualCoord( atomList[2][2] );
        
        float firstN[] = m.getActualCoord(m.residue[resAffected[0]].getAtomNameToMolnum("N"));
		float firstCA[] =  m.getActualCoord(m.residue[resAffected[0]].getAtomNameToMolnum("CA"));
		float firstC[] =  m.getActualCoord(m.residue[resAffected[0]].getAtomNameToMolnum("C"));
		
        float lastN[] = m.getActualCoord(m.residue[resAffected[2]].getAtomNameToMolnum("N"));
		float lastCA[] =  m.getActualCoord(m.residue[resAffected[2]].getAtomNameToMolnum("CA"));
		float lastC[] =  m.getActualCoord(m.residue[resAffected[2]].getAtomNameToMolnum("C"));

		float lastO[] =  m.getActualCoord(m.residue[resAffected[2]].getAtomNameToMolnum("O"));


        int numSoln = tc.solve_3pep_poly( NAnchorN, NAnchorCA, CAnchorCA, CAnchorC, r_soln_n, r_soln_a, r_soln_c);
        solutions = numSoln;
        
        solutionDihedrals = storeDihedrals(numSoln, r_soln_n, r_soln_a, r_soln_c);
        if(solutionDihedrals == null)
        {
        	System.out.println("invalid residue positions; no anchors.");
        	numSoln = 0;
        }
        
        //This is to make room for the unperturbed state at matrices[0]


            if(numSoln == 0){//At least the unperturbed state should have been generated
                System.out.println("No solutions found. Error applying tripeptide closure to residues " + resAffected[0] + " through " + resAffected[2]);
            }


        return numSoln;

    }


    
    private int[] concatenateArrays(int[]... toMerge){
        //Concatenates integer arrays (e.g. atom lists)

        int fullSize = 0;

        for ( int[] arr : toMerge )
            fullSize += arr.length;

        int[] ans = new int[fullSize];

        int cursor = 0;

        for( int[] arr : toMerge ){
            System.arraycopy(arr, 0, ans, cursor, arr.length);
            cursor += arr.length;
        }

        return ans;
    }






}
