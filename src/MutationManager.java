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
// MutationManager.java
//
//  Version:           1.0
//
//
// authors:
//    initials    name                 organization                email
//   ---------   --------------      ------------------------     ------------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//
////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Written by Ryan Lilien (2002-2004) and Ivelin Georgiev (2004-2009)
 *
 */

import java.io.*;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.util.Hashtable;
import java.math.*;
import mpi.*;

/**
 * The MutationManager class maintains a list of mutations to be tested, maintains
 *  their scores, prints a log file, and generally manages the mutations to test.
 */
public class MutationManager
{
	
    //the algorithm options that define what pruning criteria will be applied
	//	NOTE!!! These must be the same as in KSParser.java
    final int optSimple = 1;
    final int optBounds = 2;
    final int optSplit = 3;
    final int optPairs = 4;

	// Information needed by all mutations
    CommucObj cObjArray[] = null;
	int residueMap[] = null;
	String resDefault[] = null;
	int rotamerIndexOffset[] = null;
	String ligType = null;
	boolean ligPresent = false;
	int numMutations = 0;
	String arpFilenameMin = null;
	String arpFilenameMax = null;
	boolean resMutatable[][] = null;
	int algOption = 0;
	int numSplits = 0;
	String AAallowed[] = null;
	String minDEEfile = null;
	float initEw = 0.0f;
	float pruningE = (float)Math.pow(10,38);
	double gamma = 0.01;  // The gamma used in inter-mutation pruning
	float epsilon = 0.03f;  // The epsilon used in intra-mutation pruning
	float stericThresh = -10000.0f;
	float softStericThresh = -10000.0f;
	int numInAS = 0;
	int numTotalRotamers = 152;
	int numResAllowed = 0;
	boolean computeEVEnergy = true;
	boolean doMinimization = true;
	boolean minimizeBB = false;
	boolean doBackrubs = false;
	String backrubFile = null;
	boolean repeatSearch = true;
	boolean calculateVolumes = true;
	BigDecimal bestScore = new BigDecimal("0.0");
	BigDecimal q_L = new BigDecimal("0.0");
	ParamSet sParams = null;
	boolean approxMinGMEC = false;
	float lambda = (float)Math.pow(10,38);
	double stericE = Math.pow(10,38);
	boolean distDepDielect = true;
	double dielectConst = 1.0;
	boolean doDihedE = false;
	boolean doSolvationE = false;
	double solvScale = 1.0;
	double vdwMult = 1.0;
	boolean scaleInt = false;
	float maxIntScale = 1.0f;
	boolean useEref = false;
	float eRef[] = null;
	boolean entropyComp = false; //this *must* be false for the pairwise matrix energy computation
	float asasE[][][][] = null;
	boolean compASdist = false;
	boolean asDist[][] = null;
	float dist = 0.0f;

	PrintStream logPS = null;
	OneMutation mutArray[] = null;
	
	//Variables specific to PEM computation	
	float pairEMatrixMin[][][][][][] = null;
	float pairEMatrixMax[][][][][][] = null;
	float curMaxE = -(float)Math.pow(10,30);
	int numLigRotamers = 0;
	
	float pairEMatrixMinEntropy[][] = null;
	
	//Variables specific to distributed DACS and distributed DEE computations
	boolean prunedRot[] = null;
	String rotFile = null;
	boolean useSF = false;
	boolean splitFlags[][] = null;
	String sfFile = null;
	boolean distrDACS = false;
	boolean distrDEE = false;
	int numSpPos = -1;
	int msp[] = null;
	int typeDEE = -1;
	int initDepth = -1;
	int subDepth = -1;
	int diffFact = -1;
	double minRatioDiff = 0.0;
	BigInteger numInitUnprunedConf = null;
	String outputPruneInfo = null;
	String outputConfInfo = null;
	
	boolean PEMcomp = false; //true if PEM computation is performed; false if mut search is performed
	
	
	// Generic constructor
	MutationManager(String logName, OneMutation mArray[], boolean PEMcomputation) {

		if (logName!=null){ //open log file for writing
			try {
				FileOutputStream fileOutputStream = new FileOutputStream(logName);
				BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
				logPS = new PrintStream( bufferedOutputStream );
			}
			catch (Exception ex) {
				System.out.println("ERROR: An exception occured while opening log file");
			}
		}

		mutArray = mArray;
		cObjArray = new CommucObj[mutArray.length];
		
		PEMcomp = PEMcomputation;
		
		distrDACS = false; //if this is a distributed DACS computation, the flag will be set with setDistrDACS()
		distrDEE = false; //if this is a distributed DEE computation, the flag will be set with setDistrDEE()
	}
	
	// Returns the next mutation packaged in a communication object
	public synchronized CommucObj getNextComObj(int curMutIndex) {
		
		CommucObj cObj = new CommucObj();
		cObj.residueMap = residueMap;
		cObj.resDefault = resDefault;
		cObj.rotamerIndexOffset = rotamerIndexOffset;
		cObj.ligPresent = ligPresent;
		cObj.ligType = ligType;
		cObj.distrDEE = distrDEE;
		cObj.distrDACS = distrDACS;
		cObj.arpFilenameMin = arpFilenameMin;
		cObj.arpFilenameMax = arpFilenameMax;
		cObj.params = sParams;
		cObj.stericThresh = stericThresh;
		cObj.softStericThresh = softStericThresh;
		cObj.numInAS = numInAS;
		cObj.numTotalRotamers = numTotalRotamers;
		cObj.numResAllowed = numResAllowed;
		cObj.computeEVEnergy = computeEVEnergy;
		cObj.doMinimization = doMinimization;
		cObj.minimizeBB = minimizeBB;
		cObj.doBackrubs = doBackrubs;
		cObj.backrubFile = backrubFile;
		cObj.calculateVolumes = calculateVolumes;
		cObj.distDepDielect = distDepDielect;
		cObj.dielectConst = dielectConst;
		cObj.doDihedE = doDihedE;
		cObj.doSolvationE = doSolvationE;
		cObj.solvScale = solvScale;
		cObj.vdwMult = vdwMult;
		cObj.PEMcomp = PEMcomp;
		
		if (PEMcomp) {//PEM computation
			
			if (!entropyComp){
				cObj.numLigRotamers = numLigRotamers;
				cObj.flagMutType = mutArray[curMutIndex].flagMutType;
				cObj.resMut = new int[numInAS];
				cObj.currentMutation = new String[numInAS];
				for(int i=0;i<numInAS;i++) {
					cObj.resMut[i] = mutArray[curMutIndex].resMut[i];
					cObj.currentMutation[i] = resDefault[i];
				}
			}
			else { //entropy energy computation run
				cObj.entropyComp = entropyComp;
				if (compASdist){ //AS-AS distance computation
					cObj.compASdist = compASdist;
					cObj.asDist = new boolean[numInAS];
					cObj.dist = dist;
				}
				else { //AS-AS or SHL-AS energy matrix computation
					cObj.flagMutType = mutArray[curMutIndex].flagMutType;
					if (cObj.flagMutType.equalsIgnoreCase("AS-AS")){ //AS-AS run
						cObj.numInAS = 2;
						cObj.residueMap = new int[2];
						cObj.residueMap[0] = mutArray[curMutIndex].resMut[0];
						cObj.residueMap[1] = mutArray[curMutIndex].resMut[1];
					}
					else if (cObj.flagMutType.equalsIgnoreCase("INTRA")){ //INTRA run
						cObj.numInAS = 1;
						cObj.residueMap = new int[1];
						cObj.residueMap[0] = curMutIndex;
					}
					else {
						System.out.println("ERROR: only AS-AS and INTRA runs allowed for the pairwise entropy matrix precomputation.");
						System.exit(1);
					}
				}
			}
		}			
		else {//mutation search
			
			if ((!distrDEE)&&(!distrDACS)){ //mutation search run, not (distributed DACS or distributed DEE)
			
				cObj.q_L = q_L;
				cObj.numMutations = numMutations;
				cObj.repeatSearch = repeatSearch;
				cObj.initEw = initEw;
				cObj.scaleInt = scaleInt;
				cObj.maxIntScale = maxIntScale;
				cObj.pruningE = pruningE;
				cObj.stericE = stericE;
				cObj.gamma = gamma;
				cObj.epsilon = epsilon;
				cObj.currentMutation = new String[numInAS];
				for(int i=0;i<numInAS;i++) {
					cObj.currentMutation[i] = mutArray[curMutIndex].resTypes[i];
				}
				cObj.bestScore = bestScore;
			}
			else {//distributed DACS or distributed DEE
				
				cObj.initEw = initEw;
				cObj.scaleInt = scaleInt;
				cObj.maxIntScale = maxIntScale;
				cObj.prunedRot = prunedRot;
				cObj.useSF = useSF;
				cObj.sfFileIn = sfFile;
				cObj.numLigRotamers = numLigRotamers;
				cObj.useEref = useEref;
				cObj.eRef = eRef;
				
				if (distrDACS){ //distributed DACS
					cObj.numMutations = numMutations;
					cObj.rotFileIn = rotFile;
					cObj.pruningE = pruningE;
					cObj.approxMinGMEC = approxMinGMEC;
					cObj.lambda = lambda;
					cObj.algOption = algOption;
					cObj.initDepth = initDepth;
					cObj.subDepth = subDepth;
					cObj.diffFact = diffFact;
					cObj.minRatioDiff = minRatioDiff;
					cObj.msp = msp;
					cObj.numInitUnprunedConf = numInitUnprunedConf;
					cObj.currentMutation = new String[numInAS];
					cObj.outputPruneInfo = outputPruneInfo;
					cObj.outputConfInfo = outputConfInfo;
					cObj.partIndex = new int[initDepth];
					for (int i=0; i<initDepth; i++)
						cObj.partIndex[i] = mutArray[curMutIndex].resMut[i];
					for (int i=0; i<numInAS; i++){
						cObj.currentMutation[i] = resDefault[i];
					}
					cObj.bestScore = bestScore;
				}
				else { //distributed DEE
					cObj.resMut = new int[mutArray[curMutIndex].resMut.length];
					for (int i=0; i<cObj.resMut.length; i++){
						cObj.resMut[i] = mutArray[curMutIndex].resMut[i];
					}
					cObj.AAallowed = new String[numInAS];
					cObj.currentMutation = new String[numInAS];
					for (int i=0; i<numInAS; i++){
						cObj.currentMutation[i] = resDefault[i];
						cObj.AAallowed[i] = AAallowed[i];
					}
					cObj.numSpPos = numSpPos;
					cObj.typeDEE = typeDEE;
				}
			}
		}
		
		cObj.mutationNumber = curMutIndex;
		curMutIndex++;
		
		return(cObj);
	}
	

	// Output a finished mutation to the results file
	public synchronized void processFinishedMutation(CommucObj cObj) {
		
		if (PEMcomp){ //energy matrix computation
			int countNewEntries = cObj.compEE.length;
			//Update the E matrices computed so far with the new computations supplied
			//	by the current cObj
			if (!entropyComp){ //PEM computation
				for (int i=0; i<countNewEntries; i++){
					pairEMatrixMin[cObj.compEE[i].i1][cObj.compEE[i].i2][cObj.compEE[i].i3][cObj.compEE[i].i4][cObj.compEE[i].i5][cObj.compEE[i].i6] = cObj.compEE[i].minE;
					if (doMinimization)
						pairEMatrixMax[cObj.compEE[i].i1][cObj.compEE[i].i2][cObj.compEE[i].i3][cObj.compEE[i].i4][cObj.compEE[i].i5][cObj.compEE[i].i6] = cObj.compEE[i].maxE;
				}
			}
			else { //entropy E matrix computation
				if (compASdist){ //AS-AS distance computation
					asDist[cObj.mutationNumber] = cObj.asDist;
				}
				else {
					if (cObj.flagMutType.equalsIgnoreCase("INTRA")){
						for (int i=0; i<countNewEntries; i++){
							int index1 = 1 + rotamerIndexOffset[cObj.compEE[i].i2] + cObj.compEE[i].i3;
							pairEMatrixMinEntropy[cObj.mutationNumber*numTotalRotamers+index1][0] = cObj.compEE[i].minE;
						}
					}
					else { //AS-AS run
						for (int i=0; i<countNewEntries; i++){
							int ind1 = -1;
							int ind2 = -1;
							int index1 = cObj.compEE[i].i1*numTotalRotamers + rotamerIndexOffset[cObj.compEE[i].i2] + cObj.compEE[i].i3;
							int index2 = cObj.compEE[i].i4*numTotalRotamers + rotamerIndexOffset[cObj.compEE[i].i5] + cObj.compEE[i].i6;
							if (index1<numTotalRotamers){
								ind1 = index1;
								ind2 = index2-numTotalRotamers;
							}
							else {
								ind1 = index2;
								ind2 = index1-numTotalRotamers;
							}
							asasE[cObj.residueMap[0]][cObj.residueMap[1]][ind1][ind2] = cObj.compEE[i].minE;
						}
					}
				}
				
				//Output mutation information to results file (for resume)
				System.out.println("MutNUM: "+cObj.mutationNumber+" produced "+countNewEntries+" new entries.");
				logPS.print("Completed mutation "+cObj.mutationNumber);
				logPS.print(" SlaveNum "+cObj.slaveNum);
				logPS.print(" Time "+(cObj.elapsedTime/60.0));
				if (!entropyComp){
					for(int i=0;i<cObj.numInAS;i++)
						logPS.print(" "+cObj.resMut[i]);
					logPS.print(" "+cObj.flagMutType);
				}
				logPS.println();
				logPS.flush();
			}			
		}
		else { //mutation search
			if ((!distrDEE)&&(!distrDACS)){ //Hybrid MinDEE-K*, not (distributed DACS or distributed DEE)
				System.out.println("MutNUM: "+cObj.mutationNumber);
				logPS.print("Completed mutation "+cObj.mutationNumber);
				BigDecimal score = new BigDecimal("0.0");
				if (cObj.q_E.compareTo(new BigDecimal("0.0")) != 0)
					score = cObj.q_EL.divide(cObj.q_E.multiply(cObj.q_L),4);
				logPS.print(" Score "+score);
				logPS.print(" Volume "+mutArray[cObj.mutationNumber].vol);
				logPS.print(" SlaveNum "+cObj.slaveNum);
				logPS.print(" Time "+(cObj.q_E_Time/60.0)+" "+(cObj.q_EL_Time/60.0));
				logPS.print(" InitBest "+cObj.bestScore);
				BigDecimal bs = cObj.bestScore;
				if (score.compareTo(cObj.bestScore) >0)
					bs = score;
				logPS.print(" FinalBest "+bs);
				for(int i=0;i<cObj.numInAS;i++)
					logPS.print(" "+cObj.currentMutation[i]);
				logPS.print(" EConfInfo "+cObj.E_searchNumConfsEvaluated+" "+cObj.E_searchNumPrunedMinDEE+" "
						+cObj.E_searchNumConfsPrunedByS+" "+cObj.E_searchNumConfsLeft);
				logPS.print(" ELConfInfo "+cObj.EL_searchNumConfsEvaluated+" "+cObj.EL_searchNumPrunedMinDEE+" "
						+cObj.EL_searchNumConfsPrunedByS+" "+cObj.EL_searchNumConfsLeft);
				logPS.print(" MinEMinimized "+cObj.bestUnBoundEMin+" "+cObj.bestBoundEMin);
				logPS.print(" MinEUnMinimized "+cObj.bestUnBoundE+" "+cObj.bestBoundE);
				logPS.print(" Partial_q_E "+cObj.q_E+" Partial_q_EL "+cObj.q_EL);
				logPS.print(" E_total "+cObj.E_searchNumConfsTotal+" EL_total "+cObj.EL_searchNumConfsTotal);	
				logPS.print(" ESecondEw "+cObj.E_repeatEw);
				logPS.print(" ELSecondEw "+cObj.EL_repeatEw);
				logPS.print(" E_allPruned "+cObj.E_allPruned+" EL_allPruned "+cObj.EL_allPruned);
				logPS.print(" q_L "+q_L);
				logPS.println();
				if (score.compareTo(bestScore) >0){
					//logPS.println("BestScoreChange "+bestScore+" to "+score);
					bestScore = score;
				}
				logPS.flush();	
			}
			else if (distrDACS){ //distributed DACS
				bestScore = bestScore.min(cObj.bestScore);
				pruningE = bestScore.floatValue();
				logPS.print("Completed mutation "+cObj.mutationNumber);
				logPS.print(" Score "+cObj.bestScore);
				logPS.print(" BestScore "+bestScore);
				logPS.print(" PartitionIndices");
				for (int i=0; i<initDepth; i++)
					logPS.print(" "+cObj.partIndex[i]);
				logPS.print(" Time "+(cObj.elapsedTime/60.0));
				logPS.println();
				logPS.flush();
				System.out.println("Partition "+cObj.mutationNumber+" done; best energy: "+cObj.bestScore);
			}
			else {//distributed DEE

				if (typeDEE==optPairs){ //pairs DEE
					boolean tmpSF[][] = null;
					tmpSF = (boolean [][])readFromFile(tmpSF,cObj.sfFileOut);
					
					for (int i1=0; i1<tmpSF.length; i1++){
						for (int i2=0; i2<tmpSF[0].length; i2++){
							if (tmpSF[i1][i2]) //get the DE pairs from this computation
								splitFlags[i1][i2] = true; 
						}
					}
					deleteFile(cObj.sfFileOut);
					outputObject(splitFlags,sfFile); //must be output after each result read
				}
				else { //singles DEE
					for (int i=0; i<prunedRot.length; i++){
						if (cObj.prunedRot[i]) //get the pruned rotamers from this computation
							prunedRot[i] = true;
					}
					outputObject(prunedRot,rotFile); //must be output after each result read
				}
			}
		}
	}
	
	private synchronized Object readFromFile(Object inObj, String inFile){
		try{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(inFile));
			inObj = in.readObject();
			in.close();
		}
		catch (Exception e){
			System.out.println(e.toString());
			System.out.println("ERROR: An exception occurred while reading from object file");
			System.exit(0);
		}
		return inObj;
	}
	
	private void outputObject(Object outObj, String outFile){
		FileOutputStream fout = null;
		FileChannel ch = null;
		FileLock lock = null;
		try{
			fout = new FileOutputStream(outFile);
			ch = fout.getChannel();
			lock = ch.lock();
			ObjectOutputStream out = new ObjectOutputStream(fout);
			out.writeObject(outObj);
			//out.close();
		}
		catch (Exception e){
			System.out.println(e.toString());
			System.out.println("ERROR: An exception occurred while writing object file");
			System.exit(0);
		}
		finally {
			try {
				lock.release();
			}
			catch (Exception e){
				System.out.println(e.toString());
				System.out.println("ERROR: unable to release lock on file "+outFile);
				System.exit(1);
			}
		}
	}
	
	private void deleteFile(String file) {

		Runtime rt = Runtime.getRuntime();
		try {
			String tmpScr = new String("tmp_del"+System.currentTimeMillis());
			FileOutputStream fout = new FileOutputStream(tmpScr);
			DataOutputStream textOut = new DataOutputStream(fout);
			textOut.writeBytes("#!/bin/sh\n");
			String s = "rm -f "+file;
			textOut.writeBytes(s);
			textOut.close();
			String cmd = new String(tmpScr);
			System.out.println(cmd);
			rt.exec("startscript "+cmd);
				// startscript, changes the file permissions to 744, executes
				//  the file, and then deletes the file
		}
		catch(Exception ex) {
			System.out.println("Exception: runtime");
			System.out.println(ex.getMessage());
		}	
	}
	
	public void closeLog() {
		if (logPS!=null){
			logPS.flush();
			logPS.close();
		}
	}
	
	public void setResidueMap(int rm[]) {
		residueMap = rm;
	}
	public void setResDefault(String rd[]) {
		resDefault = rd;
	}
	public void setRotamerIndexOffset(int rio[]) {
		rotamerIndexOffset = rio;
	}
	public void setLigType(String lt) {
		ligType = lt;
	}
	public void setLigPresent(boolean lp) {
		ligPresent = lp;
	}
	public void setNumMutations(int nm){
		numMutations = nm;
	}
	public void setResMutatable (int resMut[][]){
		resMutatable = new boolean[resMut.length][resMut[0].length];
		for (int i=0; i<resMutatable.length; i++){
			for (int j=0; j<resMutatable[0].length; j++){
				if (resMut[i][j]==1)
					resMutatable[i][j] = true;
				else
					resMutatable[i][j] = false;
			}
		}
	}
	public void setAAallowed(String aal[]){
		AAallowed = aal;
	}
	public void setarpFilenameMin(String afnm) {
		arpFilenameMin = afnm;
	}
	public void setarpFilenameMax(String afnm) {
		arpFilenameMax = afnm;
	}
	public void setAlgOption(int ao){
		algOption = ao;
	}
	public void setNumSplits(int ns){
		numSplits = ns;
	}
	public void setMinDEEFileName(String mdf) {
		minDEEfile = mdf;
	}
	public void setInitEw(float iew){
		initEw = iew;
	}
	public void setPruningE(float pe){
		pruningE = pe;
	}
	public float getPruningE(){
		return pruningE;
	}
	public void setGamma(double g) {
		gamma = g;
	}
	public void setEpsilon(float g) {
		epsilon = g;
	}
	public void setEpsilon(double g) {
		epsilon = (new Double(g)).floatValue();
	}
	public void setParams(ParamSet theParams) {
		sParams = theParams;
	}
	public void setStericThresh(float st) {
		stericThresh = st;
	}
	public void setSoftStericThresh(float st){
		softStericThresh = st;
	}
	public void setNumInAS(int nas) {
		numInAS = nas;
	}
	public void numTotalRotamers(int ntr) {
		numTotalRotamers = ntr;
	}
	public void numResAllowed(int nra) {
		numResAllowed = nra;
	}
	public void setComputeEVEnergy(boolean ceve) {
		computeEVEnergy = ceve;
	}
	public void setDoMinimization(boolean dm) {
		doMinimization = dm;
	}
	public void setMinimizeBB(boolean mbb){
		minimizeBB = mbb;
	}
	public void setDoBackrubs(boolean br){
		doBackrubs = br;
	}
	public void setBackrubFile(String brf){
		backrubFile = brf;
	}
	public void setRepeatSearch(boolean rs){
		repeatSearch = rs;
	}
	public void setCalculateVolumes(boolean cv) {
		calculateVolumes = cv;
	}
	public void setnumLigRotamers(int nlr) {
		numLigRotamers = nlr;
	}
	public void setPairEMatrixMin(float pemMin[][][][][][]){
		pairEMatrixMin = pemMin;
	}
	public void setPairEMatrixMax(float pemMax[][][][][][]){
		pairEMatrixMax = pemMax;
	}
	public void setPrunedRot(boolean pr[]){
		prunedRot = pr;
	}
	public void setRotFile(String rf){
		rotFile = rf;
	}
	public void setUseSF(boolean usf){
		useSF = usf;
	}
	public void setSpFlags(boolean spFlags[][]){
		splitFlags = spFlags;
	}
	public void setSfFile(String sff){
		sfFile = sff;
	}
	public void setDistrDACS(boolean dDACS){
		distrDACS = dDACS;
	}
	public boolean getDistrDACS(){
		return distrDACS;
	}
	public void setDistrDEE(boolean dDEE){
		distrDEE = dDEE;
	}
	public void setBestScore(BigDecimal bs){
		bestScore = bs;
	}
	public void setNumSpPos(int spp){
		numSpPos = spp;
	}
	public void setMSP(int m[]){
		msp = m;
	}
	public void setTypeDEE(int t){
		typeDEE = t;
	}
	public void setInitDepth(int id){
		initDepth = id;
	}
	public void setSubDepth(int sd){
		subDepth = sd;
	}
	public void setDiffFact(int df){
		diffFact = df;
	}
	public void setMinRatioDiff(double mrd){
		minRatioDiff = mrd;
	}
	public void setNumInitUnprunedConf(BigInteger niuc){
		numInitUnprunedConf = niuc;
	}
	public void setOutputPruneInfo(String opi){
		outputPruneInfo = opi;
	}
	public void setOutputConfInfo(String oci){
		outputConfInfo = oci;
	}
	public void setApproxMinGMEC(boolean amg){
		approxMinGMEC = amg;
	}
	public void setLambda(float l){
		lambda = l;
	}
	public void setStericE(double se){
		stericE = se;
	}
	public void setDistDepDielect(boolean ddd){
		distDepDielect = ddd;
	}
	public void setDielectConst(double dc){
		dielectConst = dc;
	}
	public void setDoDihedE(boolean dde){
		doDihedE = dde;
	}
	public void setDoSolvationE(boolean dse){
		doSolvationE = dse;
	}
	public void setSolvScale(double ss){
		solvScale = ss;
	}
	public void setVdwMult(double vm){
		vdwMult = vm;
	}
	public void setScaleInt(boolean si){
		scaleInt = si;
	}
	public void setMaxIntScale(float is){
		maxIntScale = is;
	}
	public void setUseEref(boolean uer){
		useEref = uer;
	}
	public void setEref(float er[]){
		eRef = er;
	}
	public void setLigPartFn(BigDecimal ql){
		q_L = ql;
	}
	public void setEntropyComp(boolean ec){
		entropyComp = ec;
	}
	public float [][][][][][] getMinEmatrix(){
		return pairEMatrixMin;
	}
	public float [][][][][][] getMaxEmatrix(){
		return pairEMatrixMax;
	}
	public void setPairEntropyMatrix(float aae[][][][]){
		asasE = aae;
	}
	public float [][][][] getPairEntropyEmatrix(){
		return asasE;
	}
	public void setIntraEntropyMatrixMin(float pemMin[][]){
		pairEMatrixMinEntropy = pemMin;
	}
	public void setASdistMatrix(boolean ad[][]){
		asDist = ad;
	}
	public void setASdist(float d){
		dist = d;
	}
	public void setCompASdist(boolean ad){
		compASdist = ad;
	}
	public boolean [][] getASdistMatrix(){
		return asDist;
	}
	public float [][] getMinEmatrixEntropy(){
		return pairEMatrixMinEntropy;
	}
}
