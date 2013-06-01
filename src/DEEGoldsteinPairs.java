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
//	DEEGoldsteinPairs.java
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
 * Performs DEE Goldstein pairs rotamer pruning
 * 
*/
public class DEEGoldsteinPairs {

	//two pairwise energy matrices: one for the min energies and one for the max
	private float pairwiseMinEnergyMatrix[][][][][][] = null;
	private float pairwiseMaxEnergyMatrix[][][][][][] = null;

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
	
	//number of rotamers for the current AA type at the given residue
	//int numRotForAAtypeAtRes[];
	
	//this value depends on the particular value specified in the pairwise energy matrices;
	//		in KSParser, this value is 99999.0;
	//entries with this particular value will not be examined, as they are not allowed;
	//note that when computing E intervals, if a steric is not allowed, (maxE-minE)=0,
	//		so no comparison with stericE is necessary there
	private float bigE = (float)Math.pow(10,38);
	
	private float curEw = 0.0f;	//the max allowable difference from the GMEC (checkSum<=curEw should not be pruned)
	
	//the minimum difference in the checkSum when a rotamer cannot be pruned
	private double minDiff = -(float)Math.pow(10,30);
	
	//int count = 0;
	
	//the rotamer library
	RotamerLibrary rl = null;
	
	//The system rotamer handler
	StrandRotamers sysLR = null;
	
	//The mapping from AS position to actual residue numbers
	int residueMap[] = null;
	
	//the number of runs
	int numRuns = 1;
	
	//determines if energy minimization is performed: either traditional-DEE or MinDEE is used
	boolean doMinimize = false;
	
	//the single and pair interval terms in the DE pairs MinDEE criterion
	double indIntMinDEE2Pos[][] = null;
	double pairIntMinDEE2Pos[][] = null;
	
	//determines if magic bullet or full pairs is used
	boolean magicBullet = false;
	
	//split flags for all rotamer pairs;
	//only the dead-ending pairs that contain i_r (the rotamer to be pruned) can be 
	// 		discarded from the summation; the pairs with i_t and the pairs in the interval terms (if MinDEE/BD/BRDEE)
	// 		must still be a part of the summation.
	boolean splitFlags[][] = null;
	
	//determines if split flags are used
	boolean useFlags = false;
	
	//determines which two residues are in the current pair (only for the distributed DEE)
	boolean resInPair[] = null;
	
	//determines if distributed DEE is performed
	boolean distrDEE = false;
	
	//determines if backbone minimization is performed
	boolean minimizeBB = false;
	
	//the template interval energy (0.0 if fixed backbone)
	float templateInt = 0.0f;
	
	//the max scaling factor for the interval terms
	float maxScale = 1.0f;
	
	//the current ligand amino acid index
	int ligAANum = -1;

	//constructor
	DEEGoldsteinPairs(float arpMatrix[][][][][][], float arpMatrixMax[][][][][][], int numResInActiveSite, 
			int numTotalRotamers, int numLigRotamers,
			int rotamerIndexOffset[], int resMap[], float initEw, 
			StrandRotamers systemLRot, boolean prunedRotAtRes[], boolean residueMut[],
			boolean doMin, boolean spFlags[][], boolean useSF, boolean mb, boolean dDEE, boolean minBB, boolean scaleInt, float maxSc, 
			RotamerLibrary rlP, StrandRotamers ligROT) {
		
		doMinimize = doMin;
		
		pairwiseMinEnergyMatrix = arpMatrix;
		if (doMinimize) //max matrix is different
			pairwiseMaxEnergyMatrix = arpMatrixMax;
		else //no minimization, so the same matrix
			pairwiseMaxEnergyMatrix = pairwiseMinEnergyMatrix;
		
		splitFlags = spFlags;
		eliminatedRotAtPos = prunedRotAtRes;
		rotIndOffset = rotamerIndexOffset;		
		residueMap = resMap;
		sysLR = systemLRot;
		rl = rlP;
		useFlags = useSF;
		resInPair = residueMut;
		distrDEE = dDEE;
		magicBullet = mb;
		minimizeBB = minBB;
		
		numSiteResidues = numResInActiveSite;		// tested with 9
		numTotalRot = numTotalRotamers;				// ?152?
		numLigRot = numLigRotamers;					// 0 if no ligand
		if (numLigRot>0)
			ligAANum = ligROT.getIndexOfNthAllowable(0,0);
	
		numAAtypes = new int[numSiteResidues];
		for (int i=0; i<numAAtypes.length; i++) //the number of AAs allowed for each AS residue
			numAAtypes[i] = sysLR.getNumAllowable(residueMap[i]);
		
		curEw = initEw;
		
		numRuns = 1;
		
		templateInt = 0.0f;
		if (minimizeBB) //backbone minimization, so we need the template interval energy (otherwise, templateInt will be 0.0)			
			templateInt = pairwiseMaxEnergyMatrix[pairwiseMaxEnergyMatrix.length-1][0][0][0][0][0] - pairwiseMinEnergyMatrix[pairwiseMinEnergyMatrix.length-1][0][0][0][0][0];
		
		if (!scaleInt) //no scaling of the interval terms performed
			maxScale = 1.0f;
		else
			maxScale = maxSc;
	}
	
	//return the split flags for all rotamer pairs
	public boolean[][] getSplitFlags(){
		return splitFlags;
	}

	//Compute the conformations that can be eliminated
	//Return a boolean matrix in which an element is true if
	//the corresponding rotamer pair can be eliminated, and false otherwise
	public void ComputeEliminatedRotConf(){
		
		if (doMinimize){ //compute the DE pairs MinDEE interval terms
			
			int numRes = numSiteResidues;		
			if (numLigRot!=0) //ligand is present
				numRes++;
			
			indIntMinDEE2Pos = new double[numRes][numRes];
			pairIntMinDEE2Pos = new double[numRes][numRes];
			
			for (int posNum1=0; posNum1<numRes; posNum1++){
				for (int posNum2=posNum1+1; posNum2<numRes; posNum2++){
		
					indIntMinDEE2Pos[posNum1][posNum2] = SumMaxIndInt(posNum1,posNum2);				//formula term 3
					pairIntMinDEE2Pos[posNum1][posNum2] = SumSumMaxPairInt(posNum1,posNum2);		//formula term 4
					indIntMinDEE2Pos[posNum2][posNum1] = indIntMinDEE2Pos[posNum1][posNum2];		//formula term 3
					pairIntMinDEE2Pos[posNum2][posNum1] = pairIntMinDEE2Pos[posNum1][posNum2];		//formula term 4
				}
			}
		}
			
		//Check for pairs pruning
		int numRotForCurAAatPos1;
		
		int prunedCurRun = 0;
		boolean done = false;
		numRuns = 1;
		
		while (!done){
			
			prunedCurRun = 0;
			
			System.out.println("Current run: "+numRuns);
		
			//Compute for the AS residues first
			for (int curPos1=0; curPos1<numSiteResidues; curPos1++){
				
				if ((magicBullet)||(!distrDEE)||(resInPair[curPos1])){ //mb-pairs or not distrDEE or cur res is in distr pair
				
					System.out.print("Starting AS residue "+curPos1);
					System.out.print("..");
					
					for (int curPos2=curPos1+1; curPos2<numSiteResidues; curPos2++){
						
						if ((magicBullet)||(!distrDEE)||(resInPair[curPos2])){ //mb-pairs or not distrDEE or cur res is in distr pair
							
							for (int AA1=0; AA1<numAAtypes[curPos1]; AA1++){
								
								int curAA1 = sysLR.getIndexOfNthAllowable(residueMap[curPos1],AA1);
								
								//find how many rotamers are allowed for the current AA type at the given residue;
								//note that ala and gly have 0 possible rotamers
								numRotForCurAAatPos1 = rl.getNumRotForAAtype(curAA1);
								if (numRotForCurAAatPos1==0)	//ala or gly
									numRotForCurAAatPos1 = 1;
								
								for(int curRot1=0; curRot1<numRotForCurAAatPos1; curRot1++){
									
									int index_1 = curPos1*numTotalRot + rotIndOffset[curAA1] + curRot1;
								
									if (!eliminatedRotAtPos[index_1]){//not already pruned
										
										for (int AA2=0; AA2<numAAtypes[curPos2]; AA2++){
											int curAA2 = sysLR.getIndexOfNthAllowable(residueMap[curPos2],AA2);
											int numRotForCurAAatPos2 = rl.getNumRotForAAtype(curAA2);
											if (numRotForCurAAatPos2==0)	//ala or gly
												numRotForCurAAatPos2 = 1;
											
											for(int curRot2=0; curRot2<numRotForCurAAatPos2; curRot2++){
											
												int index_2 = curPos2*numTotalRot + rotIndOffset[curAA2] + curRot2;
												
												if (!eliminatedRotAtPos[index_2]){//not already pruned
												
													if (!splitFlags[index_1][index_2]){ //rotamer pair not already pruned
											
														if (CanEliminate(curPos1, curAA1, curRot1, curPos2, curAA2, curRot2)){
															splitFlags[index_1][index_2] = true;
															splitFlags[index_2][index_1] = true;
															
															prunedCurRun++;
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
					System.out.println("done");
				}
			}
			
			//If there is a ligand, compute DEE for the lig rotamers as well
			if (numLigRot!=0){
				
				if ((magicBullet)||(!distrDEE)||(resInPair[numSiteResidues])){ //mb-pairs or not distrDEE or ligand is in distr pair
					
					System.out.print("Starting ligand run");
					System.out.print("..");
					
					for (int curPos2=0; curPos2<numSiteResidues; curPos2++){ //curPos2 is AS residue, always != ligand
						
						if ((magicBullet)||(!distrDEE)||(resInPair[curPos2])){ //mb-pairs or not distrDEE or cur res is in distr pair
							
							for (int curRotLig=0; curRotLig<numLigRot; curRotLig++){
								int indexLig = numSiteResidues*numTotalRot+curRotLig;
								
								if (!eliminatedRotAtPos[indexLig]){//not already pruned
									
									for (int AA2=0; AA2<numAAtypes[curPos2]; AA2++){
										int curAA2 = sysLR.getIndexOfNthAllowable(residueMap[curPos2],AA2);
										int numRotForCurAAatPos2 = rl.getNumRotForAAtype(curAA2);
										if (numRotForCurAAatPos2==0)	//ala or gly
											numRotForCurAAatPos2 = 1;
										
										for(int curRot2=0; curRot2<numRotForCurAAatPos2; curRot2++){
										
											int index_2 = curPos2*numTotalRot + rotIndOffset[curAA2] + curRot2;
											
											if (!eliminatedRotAtPos[index_2]){//not already pruned
											
												if (!splitFlags[indexLig][index_2]){ //rotamer pair not already pruned
								
													if (CanEliminateLig(curRotLig,curPos2,curAA2,curRot2)){
														splitFlags[indexLig][index_2] = true;
														splitFlags[index_2][indexLig] = true;
								
														prunedCurRun++;
													}
												}
											}
										}
									}
								}
							}
						}
					}
					System.out.println("done");
				}
			}
			
			System.out.println("Number of pairs pruned this run: "+prunedCurRun);
			System.out.println("DEE: The minimum difference is "+minDiff);
			System.out.println();
			
			if (prunedCurRun==0) //no rotamers pruned this run, so done
				done = true;
			else 
				numRuns++;
			
			if ((!magicBullet)&&(!distrDEE)) //full non-distributed pairs, so do not repeat (computationally-expensive)
				done = true;
		}
	}
	
	//Called only by ComputeEliminatedRotConf(.)
	private boolean CanEliminate (int posNum1, int AANumAtPos1, int rotNumAtPos1, int posNum2, 
			int AANumAtPos2, int rotNumAtPos2){
		
		double minIndVoxelE, maxIndVoxelE;
		double minShellResE, maxShellResE;
		double minPairE, maxPairE;
		double indVoxelInterval, pairVoxelInterval;
		double minDiffPairVoxelE;
		
		int index_r1, index_r2, index_t1, index_t2;
		
		double checkSum;
		
		//In the energy matrix, column 0 gives the individual energies for each r at i;
		//skip row 0, as the individual energies start from row 1 (and are all in column 0)
		index_r1 = posNum1*numTotalRot + rotIndOffset[AANumAtPos1] + rotNumAtPos1;
		index_r2 = posNum2*numTotalRot + rotIndOffset[AANumAtPos2] + rotNumAtPos2;
		minIndVoxelE = pairwiseMinEnergyMatrix[posNum1][AANumAtPos1][rotNumAtPos1][posNum1][0][0] + pairwiseMinEnergyMatrix[posNum2][AANumAtPos2][rotNumAtPos2][posNum2][0][0]; //formula term 1
		minShellResE = pairwiseMinEnergyMatrix[posNum1][AANumAtPos1][rotNumAtPos1][posNum1][0][1] + pairwiseMinEnergyMatrix[posNum2][AANumAtPos2][rotNumAtPos2][posNum2][0][1];
		minPairE = pairwiseMinEnergyMatrix[posNum1][AANumAtPos1][rotNumAtPos1][posNum2][AANumAtPos2][rotNumAtPos2];
		
		
		//if ((minIndVoxelE<=stericEThreshIntra)&&(minShellResE<=stericEThreshPair)){//check only if not an unallowed steric
		if ((!eliminatedRotAtPos[index_r1])&&(!eliminatedRotAtPos[index_r2])
				&&(!splitFlags[index_r1][index_r2])){ //not pruned
		
			if (doMinimize){ //MinDEE, so compute the interval terms
				indVoxelInterval = indIntMinDEE2Pos[posNum1][posNum2];							//formula term 3
				pairVoxelInterval = pairIntMinDEE2Pos[posNum1][posNum2];						//formula term 4
			}
			else { //traditional-DEE, so no interval terms
				indVoxelInterval = 0.0;
				pairVoxelInterval = 0.0;
			}
			
			//For the particular position, compare the energy performance (one by one)
			//of the remaining rotamer possibilities to that of the given rotamer:
			//given r at i, compare it to all t at i for pruning
			int numRotForAAatPos1;
			
			for (int AA1=0; AA1<numAAtypes[posNum1]; AA1++){
				
				int altAA1 = sysLR.getIndexOfNthAllowable(residueMap[posNum1],AA1);
				
				numRotForAAatPos1 = rl.getNumRotForAAtype(altAA1);
				if (numRotForAAatPos1==0)	//ala or gly
					numRotForAAatPos1 = 1;
				
				for (int altRot1=0; altRot1<numRotForAAatPos1; altRot1++){
							
					//if t and r are not actually the same rotamer of the same AA
					if (!((altAA1==AANumAtPos1)&&(altRot1==rotNumAtPos1))){
						
						index_t1 = posNum1*numTotalRot + rotIndOffset[altAA1] + altRot1;
						
						if ((!eliminatedRotAtPos[index_t1])){ //not pruned
							
							int numRotForAAatPos2;
							
							for (int AA2=0; AA2<numAAtypes[posNum2]; AA2++){
								
								int altAA2 = sysLR.getIndexOfNthAllowable(residueMap[posNum2],AA2);
								
								numRotForAAatPos2 = rl.getNumRotForAAtype(altAA2);
								if (numRotForAAatPos2==0)	//ala or gly
									numRotForAAatPos2 = 1;
								
								for (int altRot2=0; altRot2<numRotForAAatPos2; altRot2++){
											
									//if t and r are not actually the same rotamer of the same AA
									if (!((altAA2==AANumAtPos2)&&(altRot2==rotNumAtPos2))){
										
										index_t2 = posNum2*numTotalRot + rotIndOffset[altAA2] + altRot2;
										
										if ((!eliminatedRotAtPos[index_t2])){ //not pruned 
							
											maxIndVoxelE = pairwiseMaxEnergyMatrix[posNum1][altAA1][altRot1][posNum1][0][0] + pairwiseMaxEnergyMatrix[posNum2][altAA2][altRot2][posNum2][0][0];		//formula term 2
											maxShellResE = pairwiseMaxEnergyMatrix[posNum1][altAA1][altRot1][posNum1][0][1] + pairwiseMaxEnergyMatrix[posNum2][altAA2][altRot2][posNum2][0][1];	
											maxPairE = pairwiseMaxEnergyMatrix[posNum1][altAA1][altRot1][posNum2][altAA2][altRot2];
											
											minDiffPairVoxelE = SumMinDiffPVE(posNum1, AANumAtPos1, rotNumAtPos1, altAA1, altRot1,
													posNum2, AANumAtPos2, rotNumAtPos2, altAA2, altRot2);	//formula term 5
											
											checkSum = -templateInt + (minIndVoxelE + minShellResE + minPairE) - (maxIndVoxelE + maxShellResE + maxPairE)
														- indVoxelInterval - pairVoxelInterval + minDiffPairVoxelE;
											
											if (checkSum > curEw){
												//System.out.println(index_r+" "+index_t+" "+checkSum+" "+minIndVoxelE+" "+minShellResE+" "+maxIndVoxelE+" "+maxShellResE+" "+indVoxelInterval+" "+pairVoxelInterval+" "+minDiffPairVoxelE);
					
												return true;}//this rotamer can be pruned/eliminated
											else {
												minDiff = Math.max(minDiff,checkSum);
												
												if (magicBullet) //magic bullet pairs, so no further checks
													return false;
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
		
		//We have tried all of the other rotamers at the current position and none
		//of them is able to prune the given rotamer, so we return false
		return false;
	}
	
	////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	
	//Called only by CanEliminate(.)
	private double SumMinDiffPVE (int atPos1, int withAA1, int withRot1, int altAA1, int altRot1,
			int atPos2, int withAA2, int withRot2, int altAA2, int altRot2){
		
		double sum = 0;
		
		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<numSiteResidues; curPos++){			
				
			if ((curPos != atPos1)&&(curPos!=atPos2))
			
				sum += IndMinDiffPVE(atPos1, withAA1, withRot1, altAA1, altRot1, 
						atPos2, withAA2, withRot2, altAA2, altRot2, curPos);
		}
		
		if (numLigRot!=0){ //there is a ligand
			//get the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position j here for which to add
			sum += LigandIndMinDiffPVE (atPos1, withAA1, withRot1, altAA1, altRot1, atPos2, withAA2, withRot2, altAA2, altRot2);
		}

		return sum;
	}
	
	//Called by SumMaxMaxPVE(.)
	private double IndMinDiffPVE (int firstPos, int firstAA1, int firstRot1, int firstAltAA1, int firstAltRot1, 
			int secondPos, int secondAA1, int secondRot1, int secondAltAA1, int secondAltRot1, int thirdPos){
		
		double minE = bigE;
		double curEmin, curEmax;
		
		int index_r1, index_r2, index_t1, index_t2, index2;
		int numRotForAAatPos;
		
		//r at i
		index_r1 = firstPos*numTotalRot + rotIndOffset[firstAA1] + firstRot1;
		index_r2 = secondPos*numTotalRot + rotIndOffset[secondAA1] + secondRot1;
		
		//t at i
		index_t1 = firstPos*numTotalRot + rotIndOffset[firstAltAA1] + firstAltRot1;
		index_t2 = secondPos*numTotalRot + rotIndOffset[secondAltAA1] + secondAltRot1;
		
		boolean found = false;
		
		if (((!eliminatedRotAtPos[index_r1])&&(!eliminatedRotAtPos[index_r2]))&&
				((!eliminatedRotAtPos[index_t1])&&(!eliminatedRotAtPos[index_t2]))){ //not pruned
		
			for (int AA=0; AA<numAAtypes[thirdPos]; AA++){
				
				int curAA = sysLR.getIndexOfNthAllowable(residueMap[thirdPos],AA);
				
				numRotForAAatPos = rl.getNumRotForAAtype(curAA);
				if (numRotForAAatPos==0)	//ala or gly
					numRotForAAatPos = 1;
				
				for (int curRot=0; curRot<numRotForAAatPos; curRot++){
						
					//There is a displacement: column 0 and row 0 have special entries, 
					//so pairwise energies start from row 1, column 1
					
					//s at j
					index2 = thirdPos*numTotalRot + rotIndOffset[curAA] + curRot;
					
					if ((!eliminatedRotAtPos[index2])){ //not pruned 
						
						if ((!useFlags)||((!splitFlags[index_r1][index2])&&(!splitFlags[index_r2][index2]))){ //not using split flags or not flagged
		
							curEmin = pairwiseMinEnergyMatrix[firstPos][firstAA1][firstRot1][thirdPos][curAA][curRot] + pairwiseMinEnergyMatrix[secondPos][secondAA1][secondRot1][thirdPos][curAA][curRot];
							curEmax = pairwiseMaxEnergyMatrix[firstPos][firstAltAA1][firstAltRot1][thirdPos][curAA][curRot] + pairwiseMaxEnergyMatrix[secondPos][secondAltAA1][secondAltRot1][thirdPos][curAA][curRot];
							
							if ((curEmin-curEmax) < minE)
								minE = curEmin-curEmax;
							
							found = true;
						}
					}					
				}
			}
			
			//if(minE==bigE)//make it contribute nothing to the sum
			//	minE = 0;
		}
		
		if (!found) //no possible pairs found
			minE = 0.0; //contributes nothing to the sum
		
		return minE;
	}
	
	//Called by SumMaxMaxPVE(.)
	private double LigandIndMinDiffPVE (int firstPos, int firstAA1, int firstRot1, int firstAltAA1, int firstAltRot1, 
			int secondPos, int secondAA1, int secondRot1, int secondAltAA1, int secondAltRot1){
		
		double minE = bigE;
		double curEmin, curEmax;
		
		int index_r1, index_r2, index_t1, index_t2, index2;
		
		//r at i
		index_r1 = firstPos*numTotalRot + rotIndOffset[firstAA1] + firstRot1;
		index_r2 = secondPos*numTotalRot + rotIndOffset[secondAA1] + secondRot1;
		
		//t at i
		index_t1 = firstPos*numTotalRot + rotIndOffset[firstAltAA1] + firstAltRot1;
		index_t2 = secondPos*numTotalRot + rotIndOffset[secondAltAA1] + secondAltRot1;
		
		boolean found = false;
		
		if (((!eliminatedRotAtPos[index_r1])&&(!eliminatedRotAtPos[index_r2]))&&
				((!eliminatedRotAtPos[index_t1])&&(!eliminatedRotAtPos[index_t2]))){ //not pruned
			
			for (int curLigPos=0; curLigPos<numLigRot; curLigPos++){
				
				//s at j (the ligand residue)
				index2 = numSiteResidues*numTotalRot + curLigPos;
				
				if ((!eliminatedRotAtPos[index2])){ //not pruned 
					
					if ((!useFlags)||((!splitFlags[index_r1][index2])&&(!splitFlags[index_r2][index2]))){ //not using split flags or not flagged
				
						curEmin = pairwiseMinEnergyMatrix[firstPos][firstAA1][firstRot1][numSiteResidues][ligAANum][curLigPos] + pairwiseMinEnergyMatrix[secondPos][secondAA1][secondRot1][numSiteResidues][ligAANum][curLigPos];
						curEmax = pairwiseMaxEnergyMatrix[firstPos][firstAltAA1][firstAltRot1][numSiteResidues][ligAANum][curLigPos] + pairwiseMaxEnergyMatrix[secondPos][secondAltAA1][secondAltRot1][numSiteResidues][ligAANum][curLigPos];
					
						if ((curEmin-curEmax) < minE)
							minE = curEmin-curEmax;
						
						found = true;
					}
				}
			}
			
			//if(minE==bigE)//make it contribute nothing to the sum
			//	minE = 0;
		}
		
		if (!found)
			minE = 0.0; //contributes nothing to the sum
		
		return minE;
	}
	//////////////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////////////		
	//Called only by CanEliminate() adn CanEliminateLig()
	private double SumMaxIndInt (int withoutPos1, int withoutPos2){
		
		double sum = 0;
		
		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<numSiteResidues; curPos++){
			
			if ((curPos != withoutPos1)&&(curPos != withoutPos2)){
				
				sum += maxScale*MaxIndInt(curPos);
			}
		}
		
		if (numLigRot!=0){ //there is a ligand
			//get the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position j here for which to add
			if ((withoutPos1!=numSiteResidues)&&(withoutPos2!=numSiteResidues)) //if we are not currently checking ligand rotamers for pruning
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
		
		//s at j
		int index1 = atPos*numTotalRot + rotIndOffset[atAA] + atRot;
		
		if ((!eliminatedRotAtPos[index1])){ //not pruned 
		
			double maxE = pairwiseMaxEnergyMatrix[atPos][atAA][atRot][atPos][0][0];
			double minE = pairwiseMinEnergyMatrix[atPos][atAA][atRot][atPos][0][0];
			
			double maxShell = pairwiseMaxEnergyMatrix[atPos][atAA][atRot][atPos][0][1];
			double minShell = pairwiseMinEnergyMatrix[atPos][atAA][atRot][atPos][0][1];
			
			//if ((maxE<=stericEThreshIntra)&&(minE<=stericEThreshIntra)
			//		&&(maxShell<=stericEThreshPair)&&(minShell<=stericEThreshPair))
				return ((maxE+maxShell) - (minE+minShell));
		}
		else
			return 0.0;//contributes nothing
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
		
		if ((!eliminatedRotAtPos[index1])){ //not pruned 
		
			double maxE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][ligRot][numSiteResidues][0][0];
			double minE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][ligRot][numSiteResidues][0][0];
			
			double maxShell = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][ligRot][numSiteResidues][0][1];
			double minShell = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][ligRot][numSiteResidues][0][1];
			
			//if ((maxE<=stericEThreshIntra)&&(minE<=stericEThreshIntra)
			//		&&(maxShell<=stericEThreshPair)&&(minShell<=stericEThreshPair))
				return ((maxE+maxShell) - (minE+minShell));
		}
		else
			return 0.0;//contributes nothing
	}
	//////////////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////////////
	//Called only by CanEliminate() and CanEliminateLig()
	private double SumSumMaxPairInt(int withoutPos1, int withoutPos2){
		
		double sum = 0;
		
		//get the contribution from the active site residue rotamers
		for (int curPos1=0; curPos1<numSiteResidues; curPos1++){
			if ((curPos1 != withoutPos1)&&(curPos1 != withoutPos2)){
				for (int curPos2=0; curPos2<curPos1; curPos2++){
					if ((curPos2 != withoutPos1)&&(curPos2 != withoutPos2)){
					
						sum += maxScale*MaxPairInt(curPos1,curPos2);
					}
				}
			}
		}
		
		if (numLigRot!=0){ //there is a ligand
			//get the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position k here for which to add;
			//the range of j is the number of active site residues
			if ((withoutPos1!=numSiteResidues)&&(withoutPos2!=numSiteResidues)){ //if we are not currently checking ligand rotamers for pruning
				for (int curPos=0; curPos<numSiteResidues; curPos++){
					if ((curPos != withoutPos1)&&(curPos != withoutPos2)){
						
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
		
		//There is a displacement: colum 0 and row 0 have special entries, 
		//so pairwise energies start from row 1, column 1
		int index1 = atPos1*numTotalRot + rotIndOffset[atAA1] + atRot1;//u at k
		int index2 = atPos2*numTotalRot + rotIndOffset[atAA2] + atRot2;//s at j
		
		if (((!eliminatedRotAtPos[index1]))&&
				((!eliminatedRotAtPos[index2]))){ //not pruned 
		
			double maxE = pairwiseMaxEnergyMatrix[atPos1][atAA1][atRot1][atPos2][atAA2][atRot2];
			double minE = pairwiseMinEnergyMatrix[atPos1][atAA1][atRot1][atPos2][atAA2][atRot2];
		
			//if ((maxE<=stericEThreshPair)&&(minE<=stericEThreshPair))
				return (maxE - minE);
		}
		else
			return 0.0;//contributes nothing
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
		
		if (((!eliminatedRotAtPos[index1]))&&
				((!eliminatedRotAtPos[index2]))){ //not pruned 
		
			double maxE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][ligRot][atPos][atAA][atRot];
			double minE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][ligRot][atPos][atAA][atRot];
		
			//if ((maxE<=stericEThreshPair)&&(minE<=stericEThreshPair))
				return (maxE - minE);
		}
		else
			return 0.0;//contributes nothing
	}
	///////////////////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////////////////
	//Same as CanEliminate(), just checks the ligand rotamers for pruning
	//Called by ComputeEliminatedRotConf()
	private boolean CanEliminateLig (int curLigRot,int posNum2, int AANumAtPos2, int rotNumAtPos2){
		
		double minIndVoxelE, maxIndVoxelE;
		double minShellResE, maxShellResE;
		double minPairE, maxPairE;
		double indVoxelInterval, pairVoxelInterval;
		double minDiffPairVoxelE;
		
		int index_r1, index_r2, index_t1, index_t2;
		
		double checkSum;
		
		//In the energy matrix, column 0 gives the individual energies for each r at i;
		//skip row 0, as the individual energies start from row 1 (and are all in column 0)
		index_r1 = numSiteResidues*numTotalRot + curLigRot;
		index_r2 = posNum2*numTotalRot + rotIndOffset[AANumAtPos2] + rotNumAtPos2;
		minIndVoxelE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][curLigRot][numSiteResidues][0][0] + pairwiseMinEnergyMatrix[posNum2][AANumAtPos2][rotNumAtPos2][posNum2][0][0]; //formula term 1
		minShellResE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][curLigRot][numSiteResidues][0][1] + pairwiseMinEnergyMatrix[posNum2][AANumAtPos2][rotNumAtPos2][posNum2][0][1];
		minPairE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][curLigRot][posNum2][AANumAtPos2][rotNumAtPos2];
		
		//if ((minIndVoxelE<=stericEThreshIntra)&&(minShellResE<=stericEThreshPair)){//check only if not an unallowed steric
		if ((!eliminatedRotAtPos[index_r1])&&(!eliminatedRotAtPos[index_r2])
				&&(!splitFlags[index_r1][index_r2])){ //not pruned
		
			if (doMinimize){ //MinDEE, so compute the interval terms
				indVoxelInterval = indIntMinDEE2Pos[numSiteResidues][posNum2];							//formula term 3
				pairVoxelInterval = pairIntMinDEE2Pos[numSiteResidues][posNum2];						//formula term 4
			}
			else { //traditional-DEE, so no interval terms
				indVoxelInterval = 0.0;
				pairVoxelInterval = 0.0;
			}			
			
			//For the particular position, compare the energy performance (one by one)
			//of the remaining rotamer possibilities to that of the given rotamer:
			//given r at i, compare it to all t at i for pruning
			for (int altRot=0; altRot<numLigRot; altRot++){
							
				//if t and r are not actually the same lig rotamer
				if (curLigRot!=altRot){
					
					//at this point, we know what r at i and t at i are
					
					index_t1 = numSiteResidues*numTotalRot + altRot;
					
					if ((!eliminatedRotAtPos[index_t1])){ //not pruned 
						
						int numRotForAAatPos2;
						
						for (int AA2=0; AA2<numAAtypes[posNum2]; AA2++){
							
							int altAA2 = sysLR.getIndexOfNthAllowable(residueMap[posNum2],AA2);
							
							numRotForAAatPos2 = rl.getNumRotForAAtype(altAA2);
							if (numRotForAAatPos2==0)	//ala or gly
								numRotForAAatPos2 = 1;
							
							for (int altRot2=0; altRot2<numRotForAAatPos2; altRot2++){
										
								//if t and r are not actually the same rotamer of the same AA
								if (!((altAA2==AANumAtPos2)&&(altRot2==rotNumAtPos2))){
									
									index_t2 = posNum2*numTotalRot + rotIndOffset[altAA2] + altRot2;
									
									if ((!eliminatedRotAtPos[index_t2])){ //not pruned 
						
										maxIndVoxelE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][altRot][numSiteResidues][0][0] + pairwiseMaxEnergyMatrix[posNum2][altAA2][altRot2][posNum2][0][0]; //formula term 2
										maxShellResE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][altRot][numSiteResidues][0][1] + pairwiseMaxEnergyMatrix[posNum2][altAA2][altRot2][posNum2][0][1];	
										maxPairE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][altRot][posNum2][altAA2][altRot2];				
					
										minDiffPairVoxelE = SumMinDiffPVELig(curLigRot, altRot, posNum2, AANumAtPos2, rotNumAtPos2, altAA2, altRot2);	//formula term 5
										
										checkSum = -templateInt + (minIndVoxelE + minShellResE + minPairE) - (maxIndVoxelE + maxShellResE + maxPairE)
													- indVoxelInterval - pairVoxelInterval + minDiffPairVoxelE;
										
										if (checkSum > curEw){
											//System.out.println(checkSum+" "+minIndVoxelE+" "+minShellResE+" "+maxIndVoxelE+" "+maxShellResE+" "+indVoxelInterval+" "+pairVoxelInterval+" "+minDiffPairVoxelE);
											
											return true;}//this rotamer can be pruned/eliminated
										else {
											minDiff = Math.max(minDiff,checkSum);
											
											if (magicBullet) //magic bullet pairs, so no further checks
												return false;
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		//We have tried all of the other rotamers at the current position and none
		//of them is able to prune the given rotamer, so we return false
		return false;
	}
	
	//Same as SumMinDiffPVE(), just checks the ligand rotamers for pruning;
	//Called by CanEliminateLig()
	private double SumMinDiffPVELig (int withRot1, int altRot1, int atPos2, int withAA2, int withRot2, int altAA2, int altRot2){
		
		double sum = 0;
		
		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<numSiteResidues; curPos++){	
			
			if ((curPos!=atPos2))
			
				sum += IndMinDiffPVELig(withRot1, altRot1, atPos2, withAA2, withRot2, altAA2, altRot2, curPos);
		}
		
		return sum;
	}
	
	//Same as IndMinDiffPVE(), just checks the ligand rotamers for pruning
	//Called by SumMinDiffPVELig()
	private double IndMinDiffPVELig (int firstRot1, int firstAltRot1, int secondPos, int secondAA1, int secondRot1, 
			int secondAltAA1, int secondAltRot1, int thirdPos){
		
		double minE = bigE;
		double curEmin, curEmax;
		
		int index_r1, index_r2, index_t1, index_t2, index2;
		int numRotForAAatPos;
		
		//r at i
		index_r1 = numSiteResidues*numTotalRot + firstRot1;
		index_r2 = secondPos*numTotalRot + rotIndOffset[secondAA1] + secondRot1;
		
		//t at i
		index_t1 = numSiteResidues*numTotalRot + firstAltRot1;
		index_t2 = secondPos*numTotalRot + rotIndOffset[secondAltAA1] + secondAltRot1;
		
		boolean found = false;
		
		if (((!eliminatedRotAtPos[index_r1])&&(!eliminatedRotAtPos[index_r2]))&&
				((!eliminatedRotAtPos[index_t1])&&(!eliminatedRotAtPos[index_t2]))){ //not pruned
			
			for (int AA=0; AA<numAAtypes[thirdPos]; AA++){
				
				int curAA = sysLR.getIndexOfNthAllowable(residueMap[thirdPos],AA);
				
				numRotForAAatPos = rl.getNumRotForAAtype(curAA);
				if (numRotForAAatPos==0)	//ala or gly
					numRotForAAatPos = 1;
				
				for (int curRot=0; curRot<numRotForAAatPos; curRot++){
						
					//There is a displacement: column 0 and row 0 have special entries, 
					//so pairwise energies start from row 1, column 1
					
					//s at j
					index2 = thirdPos*numTotalRot + rotIndOffset[curAA] + curRot;
					
					if ((!eliminatedRotAtPos[index2])){ //not pruned 
						
						if ((!useFlags)||((!splitFlags[index_r1][index2])&&(!splitFlags[index_r2][index2]))){ //not using split flags or not flagged
							
							curEmin = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][firstRot1][thirdPos][curAA][curRot] + pairwiseMinEnergyMatrix[secondPos][secondAA1][secondRot1][thirdPos][curAA][curRot];
							curEmax = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][firstAltRot1][thirdPos][curAA][curRot] + pairwiseMaxEnergyMatrix[secondPos][secondAltAA1][secondAltRot1][thirdPos][curAA][curRot];
							
							if ((curEmin-curEmax) < minE)
								minE = curEmin-curEmax;
							
							found = true;
						}
					}					
				}
			}
			
			//if(minE==bigE)//make it contribute nothing to the sum
			//	minE = 0;
		}
		
		if (!found) //no possible pairs found
			minE = 0.0; //contributes nothing to the sum
		
		return minE;
	}
}
