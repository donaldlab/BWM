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
//	MinDEEIntervals.java
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
 * Computes the single and pair interval terms in the MinDEE/BD/BRDEE criteria. 
 * This class is not used for traditional DEE.
 */
public class MinDEEIntervals {
	
	//two pairwise energy matrices: one for the min energies and one for the max
	private float pairwiseMinEnergyMatrix [][][][][][] = null;
	private float pairwiseMaxEnergyMatrix [][][][][][] = null;

	//eliminated rotamers at position i, for all positions
	private boolean eliminatedRotAtPos [] = null;
	
	//number of residues under consideration
	private int numSiteResidues;
	
	//for each residue, number of possible amino acids
	private int numTotalRot;
	
	//number of possible rotamers for the ligand
	int numLigRot;
	
	//offset of the given rotamer in the total rotamer set (?152?)
	int rotIndOffset[];
	
	//the number of AA types allowed for each AS residue
	int numAAtypes[] = null;
	
	//the rotamer library
	RotamerLibrary rl = null;
	
	//The system rotamer handler
	StrandRotamers sysLR = null;
	
	//The mapping from AS position to actual residue numbers
	int residueMap[] = null;
	
	//the single and pair interval terms in the MinDEE criterion
	double indIntMinDEE[] = null;
	double pairIntMinDEE[] = null;
	
	//determines if the interval terms should be scaled
	boolean scaleInt = false;
	
	//the distance between the closest side-chain atoms of each pair of active site residues (the ligand is not considered here) 
	float asResDist[][] = null;
	
	//the max scale
	float maxScale = 1.0f;
	
	//the max distance at which an interaction between residues is counted
	final float maxDist = 9.0f;
	
	//the current ligand amino acid index
	int ligAANum = -1;

	//constructor
	MinDEEIntervals(float arpMatrix[][][][][][], float arpMatrixMax[][][][][][], int numResInActiveSite, 
			int numTotalRotamers, int numLigRotamers, int rotamerIndexOffset[], int resMap[],
			StrandRotamers systemLRot, boolean prunedRotAtRes[], boolean scInt, float dist[][], float maxSc, RotamerLibrary rlP, StrandRotamers ligROT) {
		
		pairwiseMinEnergyMatrix = arpMatrix;
		pairwiseMaxEnergyMatrix = arpMatrixMax;
		
		eliminatedRotAtPos = prunedRotAtRes;
		rotIndOffset = rotamerIndexOffset;		
		residueMap = resMap;
		sysLR = systemLRot;
		rl = rlP;
		scaleInt = scInt;
		asResDist = dist;
		maxScale = maxSc;
		
		numSiteResidues = numResInActiveSite;		
		numTotalRot = numTotalRotamers;				// ?152?
		numLigRot = numLigRotamers;					// 0 if no ligand
		if (numLigRot>0)
			ligAANum = ligROT.getIndexOfNthAllowable(0,0);
	
		numAAtypes = new int[numSiteResidues];
		for (int i=0; i<numAAtypes.length; i++) //the number of AAs allowed for each AS residue
			numAAtypes[i] = sysLR.getNumAllowable(residueMap[i]);
		
		if (!scaleInt) //no scaling performed
			maxScale = 1.0f;
	}
	
	public double [] getIndInt(){
		return indIntMinDEE;
	}
	
	public double [] getPairInt(){
		return pairIntMinDEE;
	}
	
	//Determines the scaling factor for the interval term;
	//		Uses a distance-dependent function
	private float getScaleFactor(int pos1, int pos2){
		
		if (!scaleInt)
			return 1.0f;
		else {
			if ((pos1>=numSiteResidues)||(pos2>=numSiteResidues)) //one of the residues is the ligand, so scale with maxScale
				return maxScale;
			else { //both residues are active site residues
				if (asResDist[pos1][pos2] > maxDist)
					return 0.0f;
				else {
					double k = 2.0;
					float sf = maxScale * (float)Math.pow(( 1 - (float)Math.pow( (asResDist[pos1][pos2]/maxDist),(1/k) )),k);
					return sf;
				}
			}
		}
	}
	
	//Compute the two interval terms in the summation of the MinDEE criteria
	public void compMinDEEIntervals(){
		
		int numRes = numSiteResidues;		
		if (numLigRot!=0) //ligand is present
			numRes++;
		
		indIntMinDEE = new double[numRes];
		pairIntMinDEE = new double[numRes];
		
		for (int curPos=0; curPos<numRes; curPos++){
			
			indIntMinDEE[curPos] = SumMaxIndInt(curPos);
			pairIntMinDEE[curPos] = SumSumMaxPairInt(curPos);
		}
	}
	
	//Computes the single interval for the given residue position
	private double SumMaxIndInt (int withoutPos){
		
		double sum = 0;
		
		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<numSiteResidues; curPos++){
			
			if (curPos != withoutPos){
				
				sum += getScaleFactor(curPos,withoutPos)*MaxIndInt(curPos);
			}
		}
		
		if (numLigRot!=0){ //there is a ligand
			//get the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position j here for which to add
			if (withoutPos!=numSiteResidues) //if we are not currently checking ligand rotamers for pruning
				sum += maxScale*LigandMaxIndInt();
		}
	
		return sum;
	}
	
	//Called by SumMaxIndInt(.)
	private double MaxIndInt (int atPos){
		
		double maxEInt = -999999.0;//this is OK, since E intervals are always positive
		double curEInt;
		
		int numRotForAAatPos;
		
		for (int AA=0; AA<numAAtypes[atPos]; AA++){
			
			int curAA = sysLR.getIndexOfNthAllowable(residueMap[atPos],AA);
			
			numRotForAAatPos = rl.getNumRotForAAtype(curAA);
			if (numRotForAAatPos==0)	//ala or gly
				numRotForAAatPos = 1;
			
			for (int curRot=0; curRot<numRotForAAatPos; curRot++){
			
				curEInt = IndInt(atPos, curAA, curRot);
				if (curEInt > maxEInt){
					maxEInt = curEInt;
				}
			}
		}
		
		return maxEInt;
	}

	//Called by MaxIndInt(.)
	private double IndInt (int atPos, int atAA, int atRot){
		
		int index1 = atPos*numTotalRot + rotIndOffset[atAA] + atRot;
		
		if (!eliminatedRotAtPos[index1]){
		
			double maxE = pairwiseMaxEnergyMatrix[atPos][atAA][atRot][atPos][0][0];
			double minE = pairwiseMinEnergyMatrix[atPos][atAA][atRot][atPos][0][0];
			
			double maxShell = pairwiseMaxEnergyMatrix[atPos][atAA][atRot][atPos][0][1];
			double minShell = pairwiseMinEnergyMatrix[atPos][atAA][atRot][atPos][0][1];
			
			return ((maxE+maxShell) - (minE+minShell));
		}
		else //the interval is always non-negative, so returning 0.0 if pruned is correct
			return 0.0;
	}
	
	//Called by SumMaxIndInt(.)
	private double LigandMaxIndInt(){
		
		double maxEInt = -999999.0;//this is OK, since E intervals are always positive
		double curEInt;
		
		for (int curLigPos=0; curLigPos<numLigRot; curLigPos++){
			
			curEInt = LigandIndInt(curLigPos);
			
			if (curEInt > maxEInt){
				maxEInt = curEInt;
			}
		}
		
		return maxEInt;
	}

	//Called by LigandMaxIndInt(.)
	private double LigandIndInt (int ligRot){
		
		//s at j (the ligand residue)
		int index1 = numSiteResidues*numTotalRot + ligRot;
		
		if (!eliminatedRotAtPos[index1]){
		
			double maxE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][ligRot][numSiteResidues][0][0];
			double minE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][ligRot][numSiteResidues][0][0];
			
			double maxShell = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][ligRot][numSiteResidues][0][1];
			double minShell = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][ligRot][numSiteResidues][0][1];
		
			return ((maxE+maxShell) - (minE+minShell));
		}
		else //the interval is always non-negative, so returning 0.0 if pruned is correct
			return 0.0;
	}
	//////////////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////////////
	//Computes the pair interval for the given residue position
	private double SumSumMaxPairInt(int withoutPos){
		
		double sum = 0;
		
		//get the contribution from the active site residue rotamers
		for (int curPos1=0; curPos1<numSiteResidues; curPos1++){
			if (curPos1 != withoutPos){
				for (int curPos2=0; curPos2<curPos1; curPos2++){
					if (curPos2 != withoutPos){
						
						float sf = 1.0f;
						if (scaleInt){
							float sf1 = getScaleFactor(curPos1,withoutPos);
							float sf2 = getScaleFactor(curPos2,withoutPos);
							sf = (sf1+sf2)/2.0f;
						}
					
						sum += sf*MaxPairInt(curPos1,curPos2);
					}
				}
			}
		}
		
		if (numLigRot!=0){ //there is a ligand
			//get the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position k here for which to add;
			//the range of j is the number of active site residues
			if (withoutPos!=numSiteResidues){ //if we are not currently checking ligand rotamers for pruning
				for (int curPos=0; curPos<numSiteResidues; curPos++){
					if (curPos != withoutPos){
						
						sum += maxScale*LigandMaxPairInt(curPos);
					}
				}
			}
		}
		
		return sum;
	}
	
	//Called by SumSumMaxPairInt(.)
	private double MaxPairInt (int atPos1, int atPos2){
		
		double maxEInt = -999999.0;//this is OK, since E intervals are always positive
		double curEInt;
		
		int numRotForAAatPos1;
		
		for (int AA1=0; AA1<numAAtypes[atPos1]; AA1++){
			
			int curAA1 = sysLR.getIndexOfNthAllowable(residueMap[atPos1],AA1);
			
			numRotForAAatPos1 = rl.getNumRotForAAtype(curAA1);
			if (numRotForAAatPos1==0)	//ala or gly
				numRotForAAatPos1 = 1;
		
			for (int curRot1=0; curRot1<numRotForAAatPos1; curRot1++){
				
				int numRotForAAatPos2;
				
				for (int AA2=0; AA2<numAAtypes[atPos2]; AA2++){
					
					int curAA2 = sysLR.getIndexOfNthAllowable(residueMap[atPos2],AA2);;
					
					numRotForAAatPos2 = rl.getNumRotForAAtype(curAA2);
					if (numRotForAAatPos2==0)	//ala or gly
						numRotForAAatPos2 = 1;
					
					for (int curRot2=0; curRot2<numRotForAAatPos2; curRot2++){
					
						curEInt = PairInt(atPos1, curAA1, curRot1, atPos2, curAA2, curRot2);
						if (curEInt > maxEInt){
							maxEInt = curEInt;
						}
					}
				}
			}
		}
		
		return maxEInt;
	}
	
	//Called by MaxPairInt(.)
	private double PairInt (int atPos1, int atAA1, int atRot1, int atPos2, int atAA2, int atRot2){
		
		//There is a displacement: column 0 and row 0 have special entries, 
		//so pairwise energies start from row 1, column 1
		int index1 = atPos1*numTotalRot + rotIndOffset[atAA1] + atRot1;//u at k
		int index2 = atPos2*numTotalRot + rotIndOffset[atAA2] + atRot2;//s at j
		
		if ((!eliminatedRotAtPos[index1])&&(!eliminatedRotAtPos[index2])){
		
			double maxE = pairwiseMaxEnergyMatrix[atPos1][atAA1][atRot1][atPos2][atAA2][atRot2];
			double minE = pairwiseMinEnergyMatrix[atPos1][atAA1][atRot1][atPos2][atAA2][atRot2];
		
			return (maxE - minE);
		}
		else //the interval is always non-negative, so returning 0.0 if pruned is correct
			return 0.0;
	}
	
	//Called by SumSumMaxPairInt(.)
	private double LigandMaxPairInt (int atPos){
		
		double maxEInt = -999999.0;//this is OK, since E intervals are always positive
		double curEInt;
		
		int numRotForAAatPos;
		
		for (int AA=0; AA<numAAtypes[atPos]; AA++){
			
			int curAA = sysLR.getIndexOfNthAllowable(residueMap[atPos],AA);
			
			numRotForAAatPos = rl.getNumRotForAAtype(curAA);
			if (numRotForAAatPos==0)	//ala or gly
				numRotForAAatPos = 1;
		
			for (int curRot=0; curRot<numRotForAAatPos; curRot++){
				
				for (int curLigPos=0; curLigPos<numLigRot; curLigPos++){
					
					curEInt = LigandPairInt(atPos, curAA, curRot, curLigPos);
					if (curEInt > maxEInt){
						maxEInt = curEInt;
					}
				}
			}
		}
		
		return maxEInt;
	}
	
	//Called by LigandMaxPairInt(.)
	private double LigandPairInt (int atPos, int atAA, int atRot, int ligRot){
		
		//There is a displacement: colum 0 and row 0 have special entries, 
		//so pairwise energies start from row 1, column 1
		int index1 = numSiteResidues*numTotalRot + ligRot;//u at k (the ligand residue)
		int index2 = atPos*numTotalRot + rotIndOffset[atAA] + atRot;//s at j
		
		if ((!eliminatedRotAtPos[index1])&&(!eliminatedRotAtPos[index2])){

			double maxE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][ligRot][atPos][atAA][atRot];
			double minE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][ligRot][atPos][atAA][atRot];
	
			return (maxE - minE);
		}
		else //the interval is always non-negative, so returning 0.0 if pruned is correct
			return 0.0;
	}
}
