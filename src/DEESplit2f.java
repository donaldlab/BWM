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
//	DEESplit2f.java
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
 * Performs full split-DEE (conformational splitting) with 2 plit positions
 * 
*/
public class DEESplit2f {

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
	//		in KSParser, this value is 10^38;
	//entries with this particular value will not be examined, as they are not allowed;
	//note that when computing E intervals, if a steric is not allowed, (maxE-minE)=0,
	//		so no comparison with stericE is necessary there
	private float bigE = (float)Math.pow(10,38);
	
	//steric energy that determines incompatibility of a rotamer with the template
	float stericE = bigE;
	
	//size of the pairwise energy matrix
	private int PEMsize;
	
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
	
	int numRotForRes[] = null;
	
	//the number of split positions
	int numSplits = 0;
	
	//the single and pair interval terms in the MinDEE criterion
	double indIntMinDEE[] = null;
	double pairIntMinDEE[] = null;
	
	//split flags for all rotamer pairs;
	//only the dead-ending pairs that contain i_r (the rotamer to be pruned) can be 
	// 		discarded from the summation; the pairs with i_t and the pairs in the interval terms (if MinDEE/BD/BRDEE)
	// 		must still be a part of the summation.
	boolean splitFlags[][] = null;
	
	//determines if split flags are used
	boolean useFlags = false;
	
	//the number of runs
	int numRuns = 1;
	
	//determines if energy minimization is performed: either traditional-DEE or MinDEE is used
	boolean doMinimize = false;
	
	//determines which residue is checked (only for the distributed DEE)
	boolean resInMut[] = null;
	
	//determines if distributed DEE is performed
	boolean distrDEE = false;
	
	//determines if backbone minimization is performed
	boolean minimizeBB = false;
	
	//the template interval energy (0.0 if fixed backbone)
	float templateInt = 0.0f;
	
	//the current ligand amino acid index
	int ligAANum = -1;

	//constructor
	DEESplit2f(float arpMatrix[][][][][][], float arpMatrixMax[][][][][][], int numResInActiveSite, 
			int numTotalRotamers, int numLigRotamers, int rotamerIndexOffset[], int resMap[], float initEw, 
			StrandRotamers systemLRot, boolean prunedRotAtRes[], boolean residueMut[], 
			boolean doMin, double indInt[], double pairInt[], boolean spFlags[][], boolean useSF, boolean dDEE, boolean minBB, 
			RotamerLibrary rlP, StrandRotamers ligROT) {
		
		doMinimize = doMin;
		
		//size of the pairwise energy matrix: (km+1)x(km+1)
		PEMsize = arpMatrix.length;
		
		pairwiseMinEnergyMatrix = arpMatrix;
		if (doMinimize) //max matrix is different
			pairwiseMaxEnergyMatrix = arpMatrixMax;
		else //no minimization, so the same matrix
			pairwiseMaxEnergyMatrix = pairwiseMinEnergyMatrix;
		
		eliminatedRotAtPos = prunedRotAtRes;
		splitFlags = spFlags;
		rotIndOffset = rotamerIndexOffset;		
		residueMap = resMap;
		indIntMinDEE = indInt;
		pairIntMinDEE = pairInt;
		sysLR = systemLRot;
		rl = rlP;
		useFlags = useSF;
		resInMut = residueMut;
		distrDEE = dDEE;
		minimizeBB = minBB;
		
		numSiteResidues = numResInActiveSite;		// tested with 9
		numTotalRot = numTotalRotamers;				// ?152?
		numLigRot = numLigRotamers;					// 0 if no ligand
		if (numLigRot>0)
			ligAANum = ligROT.getIndexOfNthAllowable(0,0);
		
		numAAtypes = new int[numSiteResidues];
		for (int i=0; i<numAAtypes.length; i++) //the number of AAs allowed for each AS residue
			numAAtypes[i] = sysLR.getNumAllowable(residueMap[i]);
		
		compNumRotForRes(); //compute the number of rotamers for each residue position
		
		curEw = initEw;
		
		numRuns = 1;
		
		templateInt = 0.0f;
		if (minimizeBB) //backbone minimization, so we need the template interval energy (otherwise, templateInt will be 0.0)			
			templateInt = pairwiseMaxEnergyMatrix[pairwiseMaxEnergyMatrix.length-1][0][0][0][0][0] - pairwiseMinEnergyMatrix[pairwiseMinEnergyMatrix.length-1][0][0][0][0][0];
	}
	
	//return the split flags for all rotamer pairs
	public boolean[][] getSplitFlags(){
		return splitFlags;
	}

	//Compute the conformations that can be eliminated
	//Return a boolean matrix in which an element is true if
	//the corresponding r at i can be eliminated, and false otherwise
	public boolean [] ComputeEliminatedRotConf(){
			
		int numRotForCurAAatPos;
		
		int prunedCurRun = 0;
		boolean done = false;
		numRuns = 1;
		
		while (!done){
			
			prunedCurRun = 0;
			
			System.out.println("Current run: "+numRuns);
		
			//Compute for the AS residues first
			for (int curPos=0; curPos<numSiteResidues; curPos++){
				
				if ((!distrDEE)||(resInMut[curPos])) {//not distributed DEE or the residue to be checked for distrDEE
				
					System.out.print("Starting AS residue "+curPos);
					
					for (int AA=0; AA<numAAtypes[curPos]; AA++){
						
						System.out.print(".");
						
						int curAA = sysLR.getIndexOfNthAllowable(residueMap[curPos],AA);
						
						//find how many rotamers are allowed for the current AA type at the given residue;
						//note that ala and gly have 0 possible rotamers
						numRotForCurAAatPos = rl.getNumRotForAAtype(curAA);
						if (numRotForCurAAatPos==0)	//ala or gly
							numRotForCurAAatPos = 1;
						
						for(int curRot=0; curRot<numRotForCurAAatPos; curRot++){
						
							if (!eliminatedRotAtPos[curPos*numTotalRot + rotIndOffset[curAA] + curRot]){//not already pruned
								
								if (CanEliminate(curPos, curAA, curRot)){
									eliminatedRotAtPos[curPos*numTotalRot + rotIndOffset[curAA] + curRot] = true;
									//System.out.println(curEc);
									
									prunedCurRun++;
								}
								else
									eliminatedRotAtPos[curPos*numTotalRot + rotIndOffset[curAA] + curRot] = false;
							}
						}
					}
					System.out.println("done");
				}
			}//System.out.println("Ec: "+Ec);
			
			//If there is a ligand, compute MinDEE for the lig rotamers as well
			if (numLigRot!=0){
				
				if ((!distrDEE)||(resInMut[numSiteResidues])) {//not distributed DEE or ligand to be checked for distrDEE
					
					System.out.print("Starting ligand run");
					System.out.print("..");
					for (int curRot=0; curRot<numLigRot; curRot++){
						if (!eliminatedRotAtPos[numSiteResidues*numTotalRot+curRot]){//not already pruned
							
							if (CanEliminateLig(curRot)){
								eliminatedRotAtPos[numSiteResidues*numTotalRot+curRot] = true;
								
								prunedCurRun++;
							}
							else
								eliminatedRotAtPos[numSiteResidues*numTotalRot+curRot] = false;
						}
					}
					System.out.println("done");
				}
			}
			
			System.out.println("Number of rotamers pruned this run: "+prunedCurRun);
			System.out.println("DEE: The minimum difference is "+minDiff);
			System.out.println();
			
			if (prunedCurRun==0) //no rotamers pruned this run, so done
				done = true;
			else 
				numRuns++;
		}
		
		return eliminatedRotAtPos;
	}
	
	//Called only by ComputeEliminatedRotConf(.)
	/*
	 * The logic is as follows:
	 * 	for every AS residue
	 * 		for every AA type at this residue
	 * 			for every possible rotamer at the given AA type
	 * 				check against every other possible rotamer for all AA types at the same residue;
	 * 				eliminate if provable
	 * 
	 * That is, for each residue: out of the 152 rotamer possibilities, choose 1
	 * and compare it to the other 151 until elimination can be proven or there
	 * are no more comparisons left. Repeat this for all 152 rotamers and all AS residues
	 * 
	 */
	private boolean CanEliminate (int posNum, int AANumAtPos, int rotNumAtPos){
		
		double minIndVoxelE, maxIndVoxelE;
		double minShellResE, maxShellResE;
		double indVoxelInterval, pairVoxelInterval;
		double minDiffPairVoxelE;
		
		int index_r, index_t;
		
		double checkSum;
		
		//In the energy matrix, column 0 gives the individual energies for each r at i;
		//skip row 0, as the individual energies start from row 1 (and are all in column 0)
		index_r = posNum*numTotalRot + rotIndOffset[AANumAtPos] + rotNumAtPos;
		
		if ((!eliminatedRotAtPos[index_r])){ //not already pruned
			
			minIndVoxelE = pairwiseMinEnergyMatrix[posNum][AANumAtPos][rotNumAtPos][posNum][0][0]; 	//formula term 1
			minShellResE = pairwiseMinEnergyMatrix[posNum][AANumAtPos][rotNumAtPos][posNum][0][1];
			
			if ((minIndVoxelE + minShellResE)>=stericE) //rotamer incompatible with template, so prune
				return true;				
		
			if (doMinimize){ //MinDEE, so compute the interval terms
				indVoxelInterval = indIntMinDEE[posNum];							//formula term 3
				pairVoxelInterval = pairIntMinDEE[posNum];							//formula term 4
			}
			else { //traditional-DEE, so no interval terms
				indVoxelInterval = 0.0;
				pairVoxelInterval = 0.0;
			}
			
			
			//int splitPos[] = chooseSplitPos(posNum,majorSplitPos,index_r); //choose the split position
			int splitPos[] = new int[2];
			for (int i1=0; i1<numSiteResidues; i1++){
				
				boolean found = true;
				/*for (int i=0; i<majorSplitPos.length; i++){
					if (i1==majorSplitPos[i])
						found = false;
				}*/
				if ((i1!=posNum)&&(found)){
					for (int i2=i1+1; i2<numSiteResidues; i2++){
						
						found = true;
						/*for (int i=0; i<majorSplitPos.length; i++){
							if (i2==majorSplitPos[i])
								found = false;
						}*/
						if ((i2!=posNum)&&(found)){
							
							splitPos[0] = i1;
							splitPos[1] = i2;
			
							//System.out.println("Split positions: "+splitPos[0]+" "+splitPos[1]);
							
							int numPartitions = numRotForRes[splitPos[0]]*numRotForRes[splitPos[1]];
							boolean partitionPruned[] = new boolean[numPartitions]; //prune r at i for each partition
							for (int partP=0; partP<numPartitions; partP++){
								partitionPruned[partP] = false;
							}
							
							//local split flags (for the current partition)
							//sf1 is first indexed over the first residue; sf2 is first indexed over the second residue;
							//	two arrays are used for easier mapping
							boolean sf1[][] = new boolean[numRotForRes[splitPos[0]]][numRotForRes[splitPos[1]]];
							for (int partP1=0; partP1<sf1.length; partP1++){
								for (int partP2=0; partP2<sf1[0].length; partP2++){
									sf1[partP1][partP2] = false;
								}
							}
							
							boolean sf2[][] = new boolean[numRotForRes[splitPos[1]]][numRotForRes[splitPos[0]]];
							for (int partP1=0; partP1<sf2.length; partP1++){
								for (int partP2=0; partP2<sf2[0].length; partP2++){
									sf2[partP1][partP2] = false;
								}
							}
							
							//get the mapping between rotamer number for the split residues and the index of that
							//	rotamer into the *boolean* matrix (1 is NOT added, in contrast to the PEM indices)
							int indexMap1[] = getIndexMap(splitPos[0]);
							int indexMap2[] = getIndexMap(splitPos[1]);
							
							
							//For the particular position, compare the energy performance (one by one)
							//of the remaining rotamer possibilities to that of the given rotamer:
							//given r at i, compare it to all t at i for pruning
							int numRotForAAatPos;
							
							for (int AA=0; AA<numAAtypes[posNum]; AA++){
								
								int altAA = sysLR.getIndexOfNthAllowable(residueMap[posNum],AA);
								
								numRotForAAatPos = rl.getNumRotForAAtype(altAA);
								if (numRotForAAatPos==0)	//ala or gly
									numRotForAAatPos = 1;
								
								for (int altRot=0; altRot<numRotForAAatPos; altRot++){
											
									//if t and r are not actually the same rotamer of the same AA
									if (!((altAA==AANumAtPos)&&(altRot==rotNumAtPos))){
										
										//at this point, we know what r at i and t at i are
										
										index_t = posNum*numTotalRot + rotIndOffset[altAA] + altRot;
										
										maxIndVoxelE = pairwiseMaxEnergyMatrix[posNum][altAA][altRot][posNum][0][0];		//formula term 2
										maxShellResE = pairwiseMaxEnergyMatrix[posNum][altAA][altRot][posNum][0][1];
										
										//if ((maxIndVoxelE<=stericEThreshIntra)&&(maxShellResE<=stericEThreshPair)){//check only if not an unallowed steric
										if ((!eliminatedRotAtPos[index_t])){ //not pruned 	
										
											minDiffPairVoxelE = SumMinDiffPVE(posNum, AANumAtPos, rotNumAtPos, altAA, altRot, splitPos[0], splitPos[1]);	//formula term 5
											
											
											int partCount = 0; //the visited partitions count
											int countFirst = 0; //the current rotamer for the first residue
											for (int spAA1=0; spAA1<numAAtypes[splitPos[0]]; spAA1++){//for each AA at the first splitting residue
												
												int partAA1 = sysLR.getIndexOfNthAllowable(residueMap[splitPos[0]],spAA1);
												int numRotForPartAA1 = rl.getNumRotForAAtype(partAA1);
												if (numRotForPartAA1==0)	//ala or gly
													numRotForPartAA1 = 1;
											
												for (int partRot1=0; partRot1<numRotForPartAA1; partRot1++){//for each rot for the given partitioning AA
													
													countFirst++;
													
													int index_hv1 = splitPos[0]*numTotalRot + rotIndOffset[partAA1] + partRot1; //the index of the partitioning rotamer
													
													if ((!eliminatedRotAtPos[index_hv1])&&((!useFlags)||(!splitFlags[index_r][index_hv1]))){ //not pruned 	
														
														double splitPosDiffE1 = pairwiseMinEnergyMatrix[posNum][AANumAtPos][rotNumAtPos][splitPos[0]][partAA1][partRot1] - pairwiseMaxEnergyMatrix[posNum][altAA][altRot][splitPos[0]][partAA1][partRot1]; //formula term 6
														
														int countSecond = 0; //the current rotamer for the second residue
														for (int spAA2=0; spAA2<numAAtypes[splitPos[1]]; spAA2++){//for each AA at the second splitting residue
															
															int partAA2 = sysLR.getIndexOfNthAllowable(residueMap[splitPos[1]],spAA2);
															int numRotForPartAA2 = rl.getNumRotForAAtype(partAA2);
															if (numRotForPartAA2==0)	//ala or gly
																numRotForPartAA2 = 1;
														
															for (int partRot2=0; partRot2<numRotForPartAA2; partRot2++){//for each rot for the given partitioning AA
																
																countSecond++;
																
																if (!partitionPruned[partCount]){ //only if r at i not pruned for this partition yet
																	
																	int index_hv2 = splitPos[1]*numTotalRot + rotIndOffset[partAA2] + partRot2; //the index of the partitioning rotamer
																	
																	if ((!eliminatedRotAtPos[index_hv2])&&((!useFlags)||(!splitFlags[index_r][index_hv2]))){ //not pruned 	
																		
																		double splitPosDiffE2 = pairwiseMinEnergyMatrix[posNum][AANumAtPos][rotNumAtPos][splitPos[1]][partAA2][partRot2] - pairwiseMaxEnergyMatrix[posNum][altAA][altRot][splitPos[1]][partAA2][partRot2]; //formula term 6
														
																		checkSum = -templateInt + (minIndVoxelE + minShellResE) - (maxIndVoxelE + maxShellResE)
																					- indVoxelInterval - pairVoxelInterval + minDiffPairVoxelE + splitPosDiffE1 + splitPosDiffE2;
																		
																		if (checkSum > curEw){
																			//System.out.println(index_r+" "+index_t+" "+checkSum+" "+minIndVoxelE+" "+minShellResE+" "+maxIndVoxelE+" "+maxShellResE+" "+indVoxelInterval+" "+pairVoxelInterval+" "+minDiffPairVoxelE);
																			partitionPruned[partCount] = true; //r at i can be pruned for this partition
																			
																			if (useFlags){
																				sf1[countFirst-1][countSecond-1] = true;//incremented above, so subtract 1
																				sf2[countSecond-1][countFirst-1] = true;
																			}
																		}										
																			
																		else {
																			minDiff = Math.max(minDiff,checkSum);
																		}
																	}
																	else //v at h is pruned, so we set the pruning of r at i for this partition to true
																		partitionPruned[partCount] = true;
																}
																partCount++; //the next partition
															}
														}
													}
													else {//v at h is pruned, so we set the pruning of r at i for this partition to true
														//all of the partitions that involve hv1 are pruned; also re-synchronize the count
														for (int synchCount=0; synchCount<numRotForRes[splitPos[1]]; synchCount++){
															partitionPruned[partCount] = true;
															partCount++;
														}
													}
												}
											}
											
											if (useFlags){
												//check for global split flags
												for (int partP1=0; partP1<sf1.length; partP1++){
													int sfCount = 0;
													for (int partP2=0; partP2<sf1[0].length; partP2++){
														if (sf1[partP1][partP2])
															sfCount++;
													}
													if (sfCount==sf1[0].length){ //r at i is pruned for all subpartitions of the first residue
														splitFlags[index_r][indexMap1[partP1]] = true; //flag as a dead-ending pair
														splitFlags[indexMap1[partP1]][index_r] = true;
													}
												}
												for (int partP1=0; partP1<sf2.length; partP1++){
													int sfCount = 0;
													for (int partP2=0; partP2<sf2[0].length; partP2++){
														if (sf2[partP1][partP2])
															sfCount++;
													}
													if (sfCount==sf2[0].length){ //r at i is pruned for all subpartitions of the second residue
														splitFlags[index_r][indexMap2[partP1]] = true; //flag as a dead-ending pair
														splitFlags[indexMap2[partP1]][index_r] = true;
													}
												}
											}
											
											//after checking all partitions, return if r at i is pruned for all of the partitions;
											//	otherwise, go to the next competitor t at i
											boolean canPrune = true;
											for (int curPartCheck=0; curPartCheck<partitionPruned.length; curPartCheck++){
												if (!partitionPruned[curPartCheck]){
													canPrune = false;
													break;
												}								
											}
											if (canPrune)
												return true;
										}
									}
								}
							}
						}
					}
				}
			}
		}
		else //aready pruned
			return true;
		
		//We have tried all of the other rotamers at the current position and none
		//of them is able to prune the given rotamer, so we return false
		return false;
	}
	
	//Get the mapping between rotamer indices (into the pruning matrix) and the number of the
	//	current rotamer for the giveen residue; assumes sysLR in rs is valid (all allowables for the AS residues)
	private int [] getIndexMap(int curRes){
		
		int indexMap[] = new int[numRotForRes[curRes]];
		int indNum = 0;
		for (int AA=0; AA<sysLR.getNumAllowable(residueMap[curRes]); AA++){ //for each AA for the given AS residue
			int curAA = sysLR.getIndexOfNthAllowable(residueMap[curRes],AA);
			int numRotForAA = rl.getNumRotForAAtype(curAA);
			if (numRotForAA==0) //GLY or ALA
				numRotForAA = 1;
			
			for (int curRot=0; curRot<numRotForAA; curRot++){ //for each rot for the given AA
				indexMap[indNum] = curRes*numTotalRot + rotIndOffset[curAA] + curRot;
				indNum++;
			}
		}
		return indexMap;
	}
	
	//Compute the number of rotamers for each residue position (assign to numRotForRes[])
	private void compNumRotForRes(){
		
		boolean ligPresent = (numLigRot==0); //ligand present
		int treeLevels = numSiteResidues;
		if (ligPresent)
			treeLevels++;
		
		numRotForRes = new int[treeLevels];
		
		int curNumRot = 0;
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			if ((ligPresent)&&(curLevel==(treeLevels-1))){ //the ligand level
				curNumRot = numLigRot;
			}
			else { //AS residue				
				curNumRot = 0;
				for (int i=0; i<sysLR.getNumAllowable(residueMap[curLevel]); i++){ //add the rot for all allowable AA at this residue
					int newRot = rl.getNumRotForAAtype(sysLR.getIndexOfNthAllowable(residueMap[curLevel],i));
					if (newRot==0) //GLY or ALA
						newRot = 1;
					curNumRot += newRot; 
				}
			}
			numRotForRes[curLevel] = curNumRot;
		}
	}
	
	////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	
	//Called only by CanEliminate(.)
	private double SumMinDiffPVE (int atPos, int withAA1, int withRot1, int withAA2, int withRot2, int splitPos1, int splitPos2){
		
		double sum = 0;
		
		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<numSiteResidues; curPos++){			
				
			if ((curPos != atPos)&&(curPos!=splitPos1)&&(curPos!=splitPos2)) // j!=i and j!=k
			
				sum += IndMinDiffPVE(atPos, withAA1, withRot1, withAA2, withRot2, curPos);
		}
		
		if (numLigRot!=0){ //there is a ligand
			//get the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position j here for which to add
			sum += LigandIndMinDiffPVE (atPos, withAA1, withRot1, withAA2, withRot2);
		}

		return sum;
	}
	
	//Called by SumMaxMaxPVE(.)
	private double IndMinDiffPVE (int firstPos, int firstAA1, int firstRot1, int firstAA2, int firstRot2, int secondPos){
		
		double minE = bigE;
		double curEmin, curEmax;
		
		int index1, index2, index3;
		int numRotForAAatPos;
		
		//r at i
		index1 = firstPos*numTotalRot + rotIndOffset[firstAA1] + firstRot1;
		
		//t at i
		index3 = firstPos*numTotalRot + rotIndOffset[firstAA2] + firstRot2;
		
		boolean found = false;
		
		if (((!eliminatedRotAtPos[index1]))&&
				((!eliminatedRotAtPos[index3]))){ //not pruned 
		
			for (int AA=0; AA<numAAtypes[secondPos]; AA++){
				
				int curAA = sysLR.getIndexOfNthAllowable(residueMap[secondPos],AA);;
				
				numRotForAAatPos = rl.getNumRotForAAtype(curAA);
				if (numRotForAAatPos==0)	//ala or gly
					numRotForAAatPos = 1;
				
				for (int curRot=0; curRot<numRotForAAatPos; curRot++){
						
					//There is a displacement: column 0 and row 0 have special entries, 
					//so pairwise energies start from row 1, column 1
					
					//s at j
					index2 = secondPos*numTotalRot + rotIndOffset[curAA] + curRot;
					
					if ((!eliminatedRotAtPos[index2])){ //not pruned 
						
						if ((!useFlags)||(!splitFlags[index1][index2])){ //not using split flags or not flagged
		
							curEmin = pairwiseMinEnergyMatrix[firstPos][firstAA1][firstRot1][secondPos][curAA][curRot];
							curEmax = pairwiseMaxEnergyMatrix[firstPos][firstAA2][firstRot2][secondPos][curAA][curRot];
							//if (/*(curEmin<=stericEThreshPair)&&*/(curEmax<=stericEThreshPair)){//check only if not an unallowed steric
								if ((curEmin-curEmax) < minE)
									minE = curEmin-curEmax;
							//}
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
	private double LigandIndMinDiffPVE (int firstPos, int firstAA1, int firstRot1, int firstAA2, int firstRot2){
		
		double minE = bigE;
		double curEmin, curEmax;
		
		int index1, index2, index3;
		
		//r at i
		index1 = firstPos*numTotalRot + rotIndOffset[firstAA1] + firstRot1;
		
		//t at i
		index3 = firstPos*numTotalRot + rotIndOffset[firstAA2] + firstRot2;
		
		boolean found = false;
		
		if (((!eliminatedRotAtPos[index1]))&&
				((!eliminatedRotAtPos[index3]))){ //not pruned 
			
			for (int curLigPos=0; curLigPos<numLigRot; curLigPos++){
				
				//s at j (the ligand residue)
				index2 = numSiteResidues*numTotalRot + curLigPos;
				
				if ((!eliminatedRotAtPos[index2])){ //not pruned 
					
					if ((!useFlags)||(!splitFlags[index1][index2])){ //not using split flags or not flagged
				
						curEmin = pairwiseMinEnergyMatrix[firstPos][firstAA1][firstRot1][numSiteResidues][ligAANum][curLigPos];
						curEmax = pairwiseMaxEnergyMatrix[firstPos][firstAA2][firstRot2][numSiteResidues][ligAANum][curLigPos];
						//if (/*(curEmin<=stericEThreshPair)&&*/(curEmax<=stericEThreshPair)){//check only if not an unallowed steric
							if ((curEmin-curEmax) < minE)
								minE = curEmin-curEmax;
						//}
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
	
	///////////////////////////////////////////////////////////////////////////////////////
	//Same as CanEliminate(), just checks the ligand rotamers for pruning
	//Called by ComputeEliminatedRotConf()
	private boolean CanEliminateLig (int curLigRot){
		
		double minIndVoxelE, maxIndVoxelE;
		double minShellResE, maxShellResE;
		double indVoxelInterval, pairVoxelInterval;
		double minDiffPairVoxelE;
		
		int index_r, index_t;
		
		double checkSum;
		
		//In the energy matrix, column 0 gives the individual energies for each r at i;
		//skip row 0, as the individual energies start from row 1 (and are all in column 0)
		index_r = numSiteResidues*numTotalRot + curLigRot;
		
		if ((!eliminatedRotAtPos[index_r])){ //not already pruned
			
			minIndVoxelE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][curLigRot][numSiteResidues][0][0]; //formula term 1
			minShellResE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][curLigRot][numSiteResidues][0][1];
			
			if ((minIndVoxelE + minShellResE)>=stericE) //rotamer incompatible with template, so prune
				return true;		
		
			if (doMinimize){ //MinDEE, so compute the interval terms
				indVoxelInterval = indIntMinDEE[numSiteResidues];							//formula term 3
				pairVoxelInterval = pairIntMinDEE[numSiteResidues];							//formula term 4
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
					
					index_t = numSiteResidues*numTotalRot + altRot;
					
					maxIndVoxelE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][altRot][numSiteResidues][0][0]; //formula term 2
					maxShellResE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][altRot][numSiteResidues][0][1];
					
					//if ((maxIndVoxelE<=stericEThreshIntra)&&(maxShellResE<=stericEThreshPair)){//check only if not an unallowed steric
					if ((!eliminatedRotAtPos[index_t])){ //not pruned 
					
						minDiffPairVoxelE = SumMinDiffPVELig(curLigRot, altRot);	//formula term 5
						
						checkSum = -templateInt + (minIndVoxelE + minShellResE) - (maxIndVoxelE + maxShellResE)
									- indVoxelInterval - pairVoxelInterval + minDiffPairVoxelE;
						
						if (checkSum > curEw){
							//System.out.println(checkSum+" "+minIndVoxelE+" "+minShellResE+" "+maxIndVoxelE+" "+maxShellResE+" "+indVoxelInterval+" "+pairVoxelInterval+" "+minDiffPairVoxelE);
												
							return true;}//this rotamer can be pruned/eliminated
						else {
							minDiff = Math.max(minDiff,checkSum);
						}
					}
				}
			}
		}
		else //aready pruned
			return true;
		
		//We have tried all of the other rotamers at the current position and none
		//of them is able to prune the given rotamer, so we return false
		return false;
	}
	
	//Same as SumMinDiffPVE(), just checks the ligand rotamers for pruning;
	//Called by CanEliminateLig()
	private double SumMinDiffPVELig (int withRot1, int withRot2){
		
		double sum = 0;
		
		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<numSiteResidues; curPos++){			
			
				sum += IndMinDiffPVELig(withRot1, withRot2, curPos);
		}
		
		return sum;
	}
	
	//Same as IndMinDiffPVE(), just checks the ligand rotamers for pruning
	//Called by SumMinDiffPVELig()
	private double IndMinDiffPVELig (int firstRot1, int firstRot2, int secondPos){
		
		double minE = bigE;
		double curEmin, curEmax;
		
		int index1, index2, index3;
		int numRotForAAatPos;
		
		//r at i
		index1 = numSiteResidues*numTotalRot + firstRot1;
		
		//t at i
		index3 = numSiteResidues*numTotalRot + firstRot2;
		
		boolean found = false;
		
		if (((!eliminatedRotAtPos[index1]))&&
				((!eliminatedRotAtPos[index3]))){ //not pruned 
			
			for (int AA=0; AA<numAAtypes[secondPos]; AA++){
				
				int curAA = sysLR.getIndexOfNthAllowable(residueMap[secondPos],AA);
				
				numRotForAAatPos = rl.getNumRotForAAtype(curAA);
				if (numRotForAAatPos==0)	//ala or gly
					numRotForAAatPos = 1;
				
				for (int curRot=0; curRot<numRotForAAatPos; curRot++){
						
					//There is a displacement: column 0 and row 0 have special entries, 
					//so pairwise energies start from row 1, column 1
					
					//s at j
					index2 = secondPos*numTotalRot + rotIndOffset[curAA] + curRot;
					
					if ((!eliminatedRotAtPos[index2])){ //not pruned 
						
						if ((!useFlags)||(!splitFlags[index1][index2])){ //not using split flags or not flagged
		
							curEmin = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][firstRot1][secondPos][curAA][curRot];
							curEmax = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][firstRot2][secondPos][curAA][curRot];
							//if (/*(curEmin<=stericEThreshPair)&&*/(curEmax<=stericEThreshPair)){//check only if not an unallowed steric
								if ((curEmin-curEmax) < minE)
									minE = curEmin-curEmax;
							//}
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
