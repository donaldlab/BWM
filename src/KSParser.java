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
// KSParser.java
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

/**
 * Written by Ryan Lilien (2001-2004) and Ivelin Georgiev (2004-2009)
 * 
 */

import java.io.*;
import java.nio.channels.*;
import java.lang.Runtime;
import java.util.*;
import java.lang.Integer;
import java.math.*;
import java.lang.management.*; // added by Swati for computing CPU, USer and System time

import BranchDecomposition.BranchTree;
// import com.neva.*;   // Not compatible with linux

import mpi.MPI;
import mpi.MPIException;
import mpi.Status;

/**
 * 
 * The main class that sets up and handles the basic OSPREY computation and related functions.
 * 
 * The OSPREY functions include:
 * 		doDEE - perform DEE/A* redesign (this includes MinDEE, BD, and BRDEE);
 * 		genStructDEE - generate structures for a selected set of the top doDEE conformations;
 * 		precomputeBackrubs - precompute a list of allowed backrubs for each flexible residue position (used by BRDEE);
 * 		KSMaster - perform K* redesign;
 * 		doSinglePartFn - generate (bound or unbound) structures for the K* ensemble of a given protein-ligand complex;
 * 		doResEntropy - use SCMF to compute the residue entropy for each (non-Pro) residue in a protein.
 *
 */
public class KSParser
{
		
		boolean printSegID = false;
		
		boolean hElect = true; // should hydrogens be used in electrostatic energy calculations
		boolean hVDW = true; // should hydrogens be used in vdw energy calculations
		boolean hSteric = false; // should hydrogens be used in steric checks
		
		final double constRT = 1.9891/1000.0 * 298.15;   // in kCal / kelvin-mole (the T value here should be consistent with T in EEF1)
 	   
	// Config file name
		String cfgName = "KStar.cfg";
		
		ParamSet rParams = null; //the main KStar parameters
	
    // For soft vdw potential
		double softvdwMultiplier = 1.0;
	// For electrostatics
		boolean distDepDielect = true;
		double dielectConst = 1.0;
		
		//Determine if dihedral and solvation energies should be computed
		boolean doDihedE = false;
		boolean doSolvationE = false;
		double solvScale = 1.0;
		float stericThresh = -10000.0f; // allowed overlap between the vdW radii of two atoms, steric clash if larger overlap
		float softStericThresh = -10000.0f; // soft steric overlap threshold
		
		//the rotamer libraries
		RotamerLibrary rl = null; //the standard rotamer library for the protein (the Aa library)
		RotamerLibrary grl = null; //the rotamer library for the ligand (could be the AA or non-AA library)
	
		int numAAallowed = -1; //number of allowed AA types
		String resAllowed[] = null; //the type of allowed AA
		int rotamerIndexOffset[] = null; //the rotamer index offset for each allowed AA
		int totalNumRotamers = -1; //the total number of rotamers for the Lovell rotamer library
				
		
		final int regTag = 1; //regular tag for MPI messages
		final int updateTag = 2; //used in DACS for updating the best energy found for the different partitions
		int numProc = 1; //number of processors for MPI
		
		boolean mpiRun = false; //determines if this is an MPI run
	    
	    //the algorithm options that define what pruning criteria will be applied
		//	NOTE!!! These must be the same as in MutationManager.java
	    final int optSimple = 1;
	    final int optBounds = 2;
	    final int optSplit = 3;
	    final int optPairs = 4;
	    
	    //the assigned protein, ligand, and cofactor strand numbers
	    final int sysStrNum = 0; //the protein strand number is always 0
	    int ligStrNum = -1;
	    int cofStrNum = -1;
	 
	
	/** 
	 * Checks if this is an MPI run and calls the respective functions
	 * 
	 */
	public void checkMPI(String[] args) {
		
		if ((args.length>0)&&(args[0].equalsIgnoreCase("mpi"))) { //MPI run
			
			mpiRun = true;
			
			String tmp[] = new String[args.length-1]; //remove the mpi argument
			System.arraycopy(args, 1, tmp, 0, tmp.length);
			args = tmp;
			
			if ((args.length>0)&&(args[0].equalsIgnoreCase("-c"))){
				cfgName = args[1];
				String temp []= new String[args.length-2];
				System.arraycopy(args,2,temp,0,args.length-2);
				args = temp;
			}
			
			try{ handleDoMPI(args);} catch (Exception e){};
		}
		else {
			mpiRun = false;
			
			if ((args.length>0)&&(args[0].equalsIgnoreCase("-c"))){
				cfgName = args[1];
				String temp []= new String[args.length-2];
				System.arraycopy(args,2,temp,0,args.length-2);
				args = temp;
			}
			
			outputProgInfo(); //output program information
			setConfigPars(); //set the parameters from the configuration file
			
			parse(args); //parse the arguments
		}
	}
	
	/**
	* The main function which handles the OSPREY commands
	*/
	public void parse(String[] args) {
		
		boolean commandLineScript = false;
		boolean firstCommandLine = false;
		byte bytebuff[];
		String s = new String("");  // line being parsed
		
		if (args.length > 0) {
			commandLineScript = true;
			firstCommandLine = true;
		}
			
		bytebuff = new byte[150];
		if (!commandLineScript) {
			System.out.print("> ");
			try {
				System.in.read(bytebuff);
			}
			catch ( Exception e ){
				System.out.println("ERROR: An error occurred while reading input");
				System.exit(0);
			}
			s = new String(bytebuff);  // create a string from bytebuff
		}
		else if (commandLineScript && !firstCommandLine) {
			// If you were running a command line script and the file is over then quit
			s = new String("quit");
		}
		else if (firstCommandLine) {
			s = new String("");
			for(int i=0;i<args.length;i++)
				s = s.concat(args[i] + " ");
			firstCommandLine = false;
		}			
		 
		s = s.trim();  // remove whitespace from beginning and end of line
			
		StringTokenizer st = new StringTokenizer(s," ;\t\n\r\f");
		String firstToken = new String("");
		if (st.hasMoreTokens())
			firstToken = st.nextToken();  // snag a copy of the first token

		if (firstToken.equalsIgnoreCase("doSinglePairE"))
			doSinglePairE(s,null);
		else if (firstToken.equalsIgnoreCase("doResEntropy"))
			handleDoResEntropy(s,null);
		else if (firstToken.equalsIgnoreCase("selectResidues"))
			selectResidues(s);
		else if (firstToken.equalsIgnoreCase("compStericOverlap"))
			handleCompStericOverlap(s);
		else if (firstToken.equalsIgnoreCase("precomputeBackrubs"))
			handlePrecomputeBackrubs(s);
		
		else if (firstToken.equalsIgnoreCase("doDEE"))
			handleDoDEE(s);
		else if (firstToken.equalsIgnoreCase("doBranchDGMEC")) //added for BWM to work - Swati
			handleCompBranchDGMEC(s);
		else if (firstToken.equalsIgnoreCase("genStructDEE"))
			handleMinDEEApplyRot(s);
		else if (firstToken.equalsIgnoreCase("generateRandConfs"))
			generateRandConfs(s);
		else if (firstToken.equalsIgnoreCase("fitEparams"))
			fitEparams(s);
		
		else if (firstToken.equalsIgnoreCase("doSinglePartFn"))
			handleKSTest(s);
		else if (firstToken.equalsIgnoreCase("computeEnergyMol"))
			handleComputeEnergyMol(s);
		else if (firstToken.equalsIgnoreCase("KSMaster"))
			handleKSMaster(s);
		
		else if (firstToken.equalsIgnoreCase("genBackbones"))
			generateBackbones(s);
		
		if (mpiRun){ //exit from all slave nodes
			CommucObj cObj[] = new CommucObj[1];
			cObj[0] = null;
			for (int curProc=1; curProc<numProc; curProc++){
				try {MPI.COMM_WORLD.Send(cObj, 0, 1, MPI.OBJECT, curProc, regTag);} catch (Exception e){}
			}
		}

	} // End parse function	
	
	/**
	* Displays the program version and citations
	*/
	public void outputProgInfo() {
		
		System.out.println();
		System.out.println("OSPREY Protein Redesign Software Version 1.0");
		System.out.println("Copyright (C) 2001-2009 Bruce Donald Lab, Duke University");
		System.out.println("");
		System.out.println("This program is free software: you can redistribute it and/or modify");
		System.out.println("it under the terms of the GNU Lesser General Public License as");
		System.out.println("published by the Free Software Foundation, either version 3 of the"); 
		System.out.println("License, or (at your option) any later version.");
		System.out.println("");
		System.out.println("This program is distributed in the hope that it will be useful,");
		System.out.println("but WITHOUT ANY WARRANTY; without even the implied warranty of");
		System.out.println("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the");
		System.out.println("GNU Lesser General Public License for more details.");
		System.out.println("");
		System.out.println("There are additional restrictions imposed on the use and distribution");
		System.out.println("of this open-source code, including: (A) this header must be included");
		System.out.println("in any modification or extension of the code; (B) you are required to");
		System.out.println("cite our papers in any publications that use this code.  The citation");
		System.out.println("for the various different modules of our software, together with a");
		System.out.println("complete list of requirements and restrictions are found in the");
		System.out.println("document license.pdf enclosed with this distribution.");
		System.out.println("");
	
		System.out.println("OSPREY running on "+numProc+" processor(s)");
		System.out.println();
		
	}
	
	//Sets the parameters from the configuration file
	public void setConfigPars() {
		
		rParams = new ParamSet();
		rParams.addParamsFromFile(cfgName);
		
		hElect = (new Boolean((String)rParams.getValue("HELECT"))).booleanValue();
		hVDW = (new Boolean((String)rParams.getValue("HVDW"))).booleanValue();
		hSteric = (new Boolean((String)rParams.getValue("HSTERIC"))).booleanValue();
		distDepDielect = (new Boolean((String)rParams.getValue("DISTDEPDIELECT"))).booleanValue();
		dielectConst = (new Double((String)rParams.getValue("DIELECTCONST"))).doubleValue();
		doDihedE = (new Boolean((String)rParams.getValue("DODIHEDE"))).booleanValue();
		doSolvationE = (new Boolean((String)rParams.getValue("DOSOLVATIONE"))).booleanValue();
		solvScale = (new Double((String)rParams.getValue("SOLVSCALE"))).doubleValue();
		softvdwMultiplier = (new Double((String)rParams.getValue("VDWMULT"))).doubleValue();
		stericThresh = (new Float((String)rParams.getValue("STERICTHRESH"))).floatValue();
		softStericThresh = (new Float((String)rParams.getValue("SOFTSTERICTHRESH"))).floatValue();
		
		rl = new RotamerLibrary((String)rParams.getValue("ROTFILE"));
		numAAallowed = rl.getNumAAallowed();
		resAllowed = rl.getAAtypesAllowed();
		rotamerIndexOffset = rl.getRotamerIndexOffset();
		totalNumRotamers = rl.getTotalNumRotamers();
	}
	
	/******************************/
	// This function returns the number of tokens in string s
	private int numTokens(String s) {
		
		int curNum = 0;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");
		
		while (st.hasMoreTokens()) {
			curNum++;
		  st.nextToken();
		}
		return(curNum);
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
				// System.out.println("ERROR: Unable to access parameter " + x);
				return(new String(""));
			}
		}
		
		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken

	// Computes the factorial of input n
	public BigInteger factorial(int n){
	
		if (n==0)
			return BigInteger.valueOf(1);
	
		return (factorial(n-1).multiply(BigInteger.valueOf(n)));
	}

	public void saveMolecule(Molecule m, String fname, float energy){
	
		try{
			FileOutputStream fileOutputStream = new FileOutputStream(fname);
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( 
	fileOutputStream);
			PrintStream printStream = new PrintStream(bufferedOutputStream);
			Hashtable params = new Hashtable(7);
			params.put("printSegID",new Boolean(printSegID));
			params.put("comment","");
			params.put("energy", energy);
			params.put("showConnect",new Boolean(false));					
			new SaveMolecule(m, printStream, params); 
			printStream.close();
		}
		catch (IOException e) {
			System.out.println("ERROR: An io exception occurred while writing file");
			System.exit(0);
		}
		catch ( Exception e ){
			System.out.println(e.toString());
			System.out.println("ERROR: An exception occurred while writing file");
			System.exit(0);
		}
	}

	// This function generates all possible combinations of n choose m
	public void generateCombinations(int residueMutatable[][], int n, int m) {
	
		int curIndex[] = new int[1];
		int curComb[] = new int[n];
		curIndex[0] = 0;
		generateCombHelper(0,n,curIndex,residueMutatable,curComb,0,m);
	}
	private void generateCombHelper(int depth, int maxDepth, int curIndex[], int
		residueMutatable[][], int curComb[], int numUsed, int maxToUse){
		
		if (depth >= maxDepth){
			if (numUsed == maxToUse) {
				for (int i=0; i<maxDepth; i++) {
					residueMutatable[curIndex[0]][i] = curComb[i];
				}
				curIndex[0]++;
			}
			return;
		}
		
		curComb[depth] = 0;
		generateCombHelper(depth+1,maxDepth,curIndex,residueMutatable,curComb,numUsed,maxToUse);

		if (numUsed < maxToUse) {
			curComb[depth] = 1;
			generateCombHelper(depth+1,maxDepth,curIndex,residueMutatable,curComb,numUsed+1,maxToUse);
		}
	}
	// end combination code
	
	//Sets the allowables for AS residue position curPos
	private void setAllowablesHelper(RotamerSearch rs, ParamSet sParams, boolean addWT, int curPos, int residueMap[], String resDefault[]){
		String tempResAllow = (String)sParams.getValue("RESALLOWED"+curPos);
		for(int q=0;q<numTokens(tempResAllow);q++)
			rs.setAllowable(residueMap[curPos],getToken(tempResAllow,q+1));
		if (addWT)
			rs.setAllowable(residueMap[curPos],resDefault[curPos]); //the default type is set last
	}
	
	
	//Sets up the molecule system and returns the number of ligand rotamers
	private int setupMolSystem(Molecule m, ParamSet sParams, boolean useLig, String ligType){
		
		ligStrNum = -1;
		cofStrNum = -1;

		try{
			FileInputStream is = new FileInputStream((String)sParams.getValue("PDBNAME"));
			new PDBChemModel(m, is);
		}
		catch (Exception e){
			System.out.println("WARNING: An error occurred while reading file");
			System.out.println(e);
			System.exit(1);
		}
		
		m.strand[sysStrNum].isProtein = true;		// main protein
		
		int strNum = sysStrNum + 1; //the current strand number; 0 is reserved for the protein strand, the ligand strand is 1 (if present)

		//Get the ligand (if present and if it will be used)
		int pdbLigNum = (new Integer((String)sParams.getValue("PDBLIGNUM")));
		if (pdbLigNum>=0){ //with ligand in PDB
			
			int molLigNum = m.mapPDBresNumToMolResNum(pdbLigNum);
			Residue lig = m.residue[molLigNum]; // pull out the ligand
			m.deleteResidue(molLigNum);
			
			if (useLig) { //ligand will be used in design
				
				lig.renumberResidue();
				m.addStrand("L");
				ligStrNum = strNum;
				m.addResidue(ligStrNum,lig,true);
				strNum++;
				
				m.strand[ligStrNum].isProtein = (new Boolean((String)sParams.getValue("LIGAA"))).booleanValue();
				if (m.strand[ligStrNum].isProtein) //use the AA rotamer library for the ligand
					grl = rl;
				else //use the non-AA rotamer library for the ligand
					grl = new RotamerLibrary((String)rParams.getValue("GROTFILE"));
				
				// change ligand to the specified residue type
				if ( m.strand[ligStrNum].isProtein && !m.strand[ligStrNum].residue[0].name.equalsIgnoreCase(ligType) ) //not the same ligand type
					(new StrandRotamers(grl,m.strand[ligStrNum])).changeResidueType(m,0,ligType,true);
			}
		}
		else if (useLig){
			System.out.println("ERROR: Attempting to use a ligand, but ligand not found in system config file");
			System.exit(1);
		}
		
		//Get the cofactor (if present)
		int numCofactorRes = (new Integer((String)sParams.getValue("NUMCOFRES"))).intValue();
		if (numCofactorRes>0) { //there is a cofactor
			
			String cofMapString = (String)sParams.getValue("COFMAP");
			Residue cof[] = new Residue[numCofactorRes];
			for(int i=0;i<numCofactorRes;i++){
				int cofactorRes = m.mapPDBresNumToMolResNum((new Integer(getToken(cofMapString,i+1))).intValue());
				cof[i] = m.residue[cofactorRes];
				m.deleteResidue(cofactorRes);
				cof[i].renumberResidue();
			}
			
			m.addStrand("M");
			cofStrNum = strNum;
			for (int i=0; i<numCofactorRes; i++)
				m.addResidue(cofStrNum,cof[i],true);
			strNum++;
			
			m.strand[cofStrNum].isProtein = false;  // cofactor
		}
		
		//Determine the number of rotamers for the ligand (if used)
		int numLigRotamers = 0;
		if (useLig) {
			numLigRotamers = grl.getNumRotamers(m.strand[ligStrNum].residue[0].name);
			if (numLigRotamers == 0)
				numLigRotamers = 1;
		}
		return numLigRotamers;
	}

	/** 
	 * Computes the bound or unbound partition function for a given single protein-ligand complex and
	 * can compute the (energy-minimized) structures for all conformations in the partition function conformational ensemble
	*/
	public void handleKSTest(String s) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutataion search parameter filename (string)
		
		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters		

		// Pull search parameters
		int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
		String eMatrixNameMin = (String)sParams.getValue("MINENERGYMATRIXNAME");
		String eMatrixNameMax = (String)sParams.getValue("MAXENERGYMATRIXNAME");
		boolean doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE"))).booleanValue();
		boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB"))).booleanValue();
		boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS"))).booleanValue();
		boolean repeatSearch = (new Boolean((String)sParams.getValue("REPEATSEARCH"))).booleanValue();
		String backrubFile = (String)sParams.getValue("BACKRUBFILE");
		boolean scaleInt = (new Boolean((String)sParams.getValue("SCALEINT"))).booleanValue();
		float maxIntScale = (new Float((String)sParams.getValue("MAXINTSCALE"))).floatValue();
		float initEw = (new Float((String)sParams.getValue("INITEW"))).floatValue();
		double pruningE = (new Double((String)sParams.getValue("PRUNINGE"))).doubleValue();
		double stericE = (new Double((String)sParams.getValue("STERICE"))).doubleValue();
		boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
		String ligType = (String)sParams.getValue("LIGTYPE");		
		boolean saveConfs = (new Boolean((String)sParams.getValue("OUTPUTPDBS"))).booleanValue();
		String fName = (String)sParams.getValue("PDBPREFIX");
		String resMut = (String)sParams.getValue("RESMUT");
		float epsilon = (new Float((String)sParams.getValue("EPSILON"))).floatValue();
		float gamma = (new Float((String)sParams.getValue("GAMMA"))).floatValue();
		BigDecimal bestScore = new BigDecimal((String)sParams.getValue("BESTSCORE"));
		BigDecimal q_E = new BigDecimal((String)sParams.getValue("PROTPARTFN"));
		
		if (!doMinimize)
			minimizeBB = false;
		if (!minimizeBB)
			doBackrubs = false;
		
		if ( (!ligPresent) && ((new Boolean((String)sParams.getValue("USEUNBOUNDSTRUCT"))).booleanValue()) ) { //ligPresent, or a different input structure is used for the unbound partition function computation
			sParams.setValue("PDBNAME",sParams.getValue("UNBOUNDPDBNAME"));
			sParams.setValue("PDBLIGNUM","-1");
			eMatrixNameMin = sParams.getValue("MINENERGYMATRIXNAMEUNBOUND");
			eMatrixNameMax = sParams.getValue("MAXENERGYMATRIXNAMEUNBOUND");
		}
		
		//Setup the molecule system
		Molecule m = new Molecule();
		int numLigRotamers = setupMolSystem(m,sParams,ligPresent,ligType);

		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.getValue("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
			residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
		}
		System.out.println();
	
		RotamerSearch rs = new RotamerSearch(m,sysStrNum,ligStrNum,hElect,hVDW,hSteric,true,true,epsilon,stericThresh,softStericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,softvdwMultiplier,rl,grl);
		
		// Define the mutation amino-acid sequence
		System.out.print("Mutation Sequence:");
		String curSeq[] = new String[numInAS];
		for(int i=0;i<numInAS;i++){
			curSeq[i] = getToken(resMut,i+1);
			System.out.print(" "+curSeq[i]);
		}
		System.out.println();

		System.out.println("Beginning setAllowables");
		for(int i=0;i<numInAS;i++){
			rs.setAllowable(residueMap[i],curSeq[i]);
		}

		System.out.print("Loading precomputed min energy matrix...");
		loadPairwiseEnergyMatrices(sParams,rs,eMatrixNameMin+".dat",doMinimize,eMatrixNameMax+".dat");
		System.out.println("done");
		
		BigDecimal q_L = BigDecimal.ZERO;
		if (ligPresent)
			q_L = getLigPartFn(m,numInAS,ligType,eMatrixNameMin+".dat"); //compute the ligand partition function
		
		System.out.println("Before start");		
			
		boolean prunedRotAtRes[] = new boolean[numInAS*totalNumRotamers+numLigRotamers];
		for (int i=0; i<prunedRotAtRes.length; i++)
			prunedRotAtRes[i] = false;
		
		//Prune all rotamers that are incompatible with the template (intra E + re-to-template E >= stericE)
		prunedRotAtRes = rs.DoPruneStericTemplate(numInAS, totalNumRotamers, numLigRotamers, 
				residueMap, rotamerIndexOffset, prunedRotAtRes, stericE);
		
		if (doMinimize) //precompute the interval terms in the MinDEE criterion
			rs.doCompMinDEEIntervals(numInAS, totalNumRotamers, numLigRotamers, residueMap, 
					rotamerIndexOffset, prunedRotAtRes, scaleInt, maxIntScale);
		
		prunedRotAtRes = rs.DoDEEGoldstein(numInAS, totalNumRotamers, numLigRotamers, residueMap,
				rotamerIndexOffset, initEw, prunedRotAtRes, doMinimize, false, minimizeBB);
		
		//Prune with MinBounds
		prunedRotAtRes = rs.DoMinBounds(numInAS,totalNumRotamers,numLigRotamers,
				residueMap,rotamerIndexOffset,pruningE,prunedRotAtRes,initEw, false, false);
		
		//Compute the Ec value and prunedIsSteric[]
		rs.DoMinBounds(numInAS,totalNumRotamers,numLigRotamers,
				residueMap,rotamerIndexOffset,pruningE,prunedRotAtRes,initEw, false, true);
		
		BigDecimal initialBest = BigDecimal.ZERO;
		if (ligPresent)
			initialBest =  q_E.multiply(bestScore.multiply(q_L)).multiply(new BigDecimal(gamma * epsilon));
		
		//Do the rotamer search
		rs.slaveDoRotamerSearch(true,doMinimize,numInAS,numAAallowed,totalNumRotamers,rotamerIndexOffset,resAllowed,
				residueMap,ligPresent,initialBest,null,minimizeBB,saveConfs,fName,doBackrubs,backrubFile);
		
		if ((repeatSearch)&&(rs.repeatSearch)){ //the desired accuracy was not achieved, so repeat the search: the setup is already done
			
			System.out.println();
			System.out.println("Repeating search..");
			rs.repeatSearch = false; //reset the flag
			rs.slaveDoRotamerSearch(true,doMinimize,numInAS,numAAallowed,totalNumRotamers,rotamerIndexOffset,resAllowed,
					residueMap,ligPresent,initialBest,null,minimizeBB,saveConfs,fName,doBackrubs,backrubFile);
		}
	}
	
	// Finds the energy for a given input system (a molecule with specified flexible residues)
	public void handleComputeEnergyMol(String s) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Ligand (boolean), is true if present
		// 3: Amino acid type for ligand (if ligand is absent, write none or anything)

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters

		// Pull search parameters
		boolean ligPresent = (new Boolean(getToken(s,3))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = getToken(s,4);
		
		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,ligPresent,ligType);
		
		System.out.println("Starting energy computation");
		Amber96ext a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, softvdwMultiplier);
		a96ff.calculateTypesWithTemplates();
		a96ff.initializeCalculation();
		a96ff.setNBEval(hElect,hVDW);
		if (ligPresent)
			a96ff.setLigandNum((new Integer((String)sParams.getValue("PDBLIGNUM"))).intValue());
		
		double energy[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1);
		System.out.println("System energy: " + energy[0]+" (elect: "+energy[1]+" vdW: "+energy[2]+" solvation: "+energy[3]+")");
	}

	/**
	 * Performs K* redesign; sets up the K* computation from the input model and configuration files and distributes the
	 * candidate mutants for evaluation by the set of available processors.
	*/
	public void handleKSMaster(String s) {
	
		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutation search parameter filename (string)
		
		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters
		
		
		int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
		int numMutations = (new Integer((String)sParams.getValue("NUMMUTATIONS"))).intValue();
		String runName = (String)sParams.getValue("RUNNAME");
		String mutFileName = (String)sParams.getValue("MUTFILENAME");
		String eMatrixNameMin = (String)sParams.getValue("MINENERGYMATRIXNAME");
		String eMatrixNameMax = (String)sParams.getValue("MAXENERGYMATRIXNAME");
		boolean doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE"))).booleanValue();
		boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB"))).booleanValue();
		boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS"))).booleanValue();
		String backrubFile = (String)sParams.getValue("BACKRUBFILE");
		boolean repeatSearch = (new Boolean((String)sParams.getValue("REPEATSEARCH"))).booleanValue();
		boolean scaleInt = (new Boolean((String)sParams.getValue("SCALEINT"))).booleanValue();
		float maxIntScale = (new Float((String)sParams.getValue("MAXINTSCALE"))).floatValue();
		float initEw = (new Float((String)sParams.getValue("INITEW"))).floatValue();
		float pruningE = (new Float((String)sParams.getValue("PRUNINGE"))).floatValue();
		double stericE = (new Double((String)sParams.getValue("STERICE"))).doubleValue();
		boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
		String ligType = (String)sParams.getValue("LIGTYPE");
		float targetVol = (new Float((String)sParams.getValue("TARGETVOLUME"))).floatValue();
		float volWindow = (new Float((String)sParams.getValue("VOLUMEWINDOW"))).floatValue();
		boolean resumeSearch = (new Boolean((String)sParams.getValue("RESUMESEARCH"))).booleanValue();
		String resumeFilename = (String)sParams.getValue("RESUMEFILENAME");
		double gamma = (new Double((String)sParams.getValue("GAMMA"))).doubleValue();
		float epsilon = (new Float((String)sParams.getValue("EPSILON"))).floatValue();
		
		if (!mpiRun){
			System.out.println("ERROR: Distributed computation requires MPI");
			System.exit(1);
		}
		
		System.out.println("Run Name: "+runName);
		System.out.println("Precomputed Min Energy Matrix: "+eMatrixNameMin);
		System.out.println("Precomputed Max Energy Matrix: "+eMatrixNameMax);
		System.out.println("Ligand Type: "+ligType);
		System.out.println("Volume Center: "+targetVol);
		System.out.println("Volume Window Size: "+volWindow);
		System.out.println("Num Residues Allowed to Mutate: "+numMutations);
		
		if(resumeSearch) {
			System.out.println("** Resuming Search **");
			System.out.println("     resuming from file: "+resumeFilename);
		}
		
		// Create the mutation list with estimated energies
		OneMutation mutArray[] = new OneMutation[200000];
		
		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,ligPresent,ligType);

		// Generate all combinations (include (n choose m), (n choose m-1), ... , (n choose 1), and (n choose 0) )
		int numCombAll = 0;
		for (int i=numMutations; i>=0; i--)
			numCombAll += factorial(numInAS).divide(factorial(numInAS-i).multiply(factorial(i))).intValue();
		int residueMutatableAll[][] = new int[numCombAll][numInAS];
		int curInd = 0;
		for (int i=numMutations; i>=0; i--){
			int numCombCur = factorial(numInAS).divide(factorial(numInAS-i).multiply(factorial(i))).intValue();
			int residueMutatableCur[][] = new int[numCombCur][numInAS];
			generateCombinations(residueMutatableCur,numInAS,i);
			for (int j=0; j<numCombCur; j++){
				residueMutatableAll[curInd] = residueMutatableCur[j];
				curInd++;
			}
		}
		
		// At this point each row of residueMutatble is a 0/1 array, 1 indicates
		//  that that residues can mutate

		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.getValue("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
			residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
		}
		System.out.println();
		
		
		System.out.print("Checking if precomputed energy matrix is already computed...");
		RotamerSearch rs = new RotamerSearch(m,sysStrNum,ligStrNum,hElect,hVDW,hSteric,true,true,epsilon,stericThresh,softStericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,softvdwMultiplier,rl,grl);
		loadPairwiseEnergyMatrices(sParams,rs,eMatrixNameMin+".dat",doMinimize,eMatrixNameMax+".dat");
		rs = null;
		System.out.println("done");
		
		if ((new Boolean((String)sParams.getValue("USEUNBOUNDSTRUCT"))).booleanValue()){ //a different input structure is used for the unbound partition function computation
			ParamSet ubParams = new ParamSet(); //create a new parameter set, just for the unbound-case matrix computation; sParams must not be changed here
			ubParams.setParamsValues(sParams.getParams(), sParams.getValues(), sParams.getCurNum());
			ubParams.setValue("PDBNAME",ubParams.getValue("UNBOUNDPDBNAME"));
			ubParams.setValue("PDBLIGNUM","-1");
			ubParams.setValue("LIGPRESENT", "false");
			ubParams.setValue("MINENERGYMATRIXNAME", ubParams.getValue("MINENERGYMATRIXNAMEUNBOUND"));
			ubParams.setValue("MAXENERGYMATRIXNAME", ubParams.getValue("MAXENERGYMATRIXNAMEUNBOUND"));
			System.out.print("Checking if precomputed energy matrix (unbound) is already computed...");
			rs = new RotamerSearch(m,sysStrNum,ligStrNum,hElect,hVDW,hSteric,true,true,epsilon,stericThresh,softStericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,softvdwMultiplier,rl,grl);
			loadPairwiseEnergyMatrices(ubParams,rs,ubParams.getValue("MINENERGYMATRIXNAME")+".dat",doMinimize,ubParams.getValue("MAXENERGYMATRIXNAME")+".dat");
			rs = null;
			ubParams = null;
			System.out.println("done");
			m = new Molecule();
			setupMolSystem(m,sParams,ligPresent,ligType); //re-initialize, since some molecule-relative variables have changed (e.g., ligStrNum)
		}
		

		//Load mutation list for distribution
		mutArray = handleHybridKSLoadMutList(mutArray, mutFileName, numInAS, m, numCombAll, residueMutatableAll,
				sParams, residueMap, resDefault, resAllowed, numAAallowed, targetVol, volWindow);

		BigDecimal bestScore = new BigDecimal("0.0"); //for the resume results		
		// If doing a resume, read the initial results into a bunch of OneMutations
		if (resumeSearch) {
		
			bestScore = new BigDecimal("0.0"); //high scores are better
			
			OneMutation resumeResults[] = new OneMutation[mutArray.length];
			for(int q=0;q<mutArray.length;q++)
				resumeResults[q] = new OneMutation();
			resumeResults = readResumeFile(resumeResults,resumeFilename,numInAS,false,false,-1);
			System.out.println("Read "+resumeResults.length+" completed mutations");
			
			// Now filter removed mutations (already computed results
			//  are NOT written to file since you already have them)
			// We do need to maintain the best score
			int newIndex = 0;
			OneMutation newArray2[] = new OneMutation[mutArray.length];
			for(int q=0;q<mutArray.length;q++) {
				int w = findMutationIndex(resumeResults,mutArray[q].resTypes);
				if (w>=0)
					bestScore = bestScore.max(resumeResults[w].score); //higher scores are better for Hybrid MinDEE-K*
				else 
					newArray2[newIndex++] = mutArray[q];
			}
			mutArray = new OneMutation[newIndex];
			System.arraycopy(newArray2,0,mutArray,0,newIndex);
			System.out.println("Length of mutArray after removing already computed mutations: "+mutArray.length);
		}
		
		BigDecimal q_L = getLigPartFn(m,numInAS,ligType,eMatrixNameMin+".dat");
	
		MutationManager mutMan = new MutationManager(runName,mutArray,false);
		mutMan.setResidueMap(residueMap);
		mutMan.setResDefault(resDefault);
		mutMan.setRotamerIndexOffset(rotamerIndexOffset);
		mutMan.setLigType(ligType);
		mutMan.setLigPresent(ligPresent);
		mutMan.setNumMutations(numMutations);
		mutMan.setarpFilenameMin(eMatrixNameMin+".dat");
		mutMan.setarpFilenameMax(eMatrixNameMax+".dat");
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setDoBackrubs(doBackrubs);
		mutMan.setBackrubFile(backrubFile);
		mutMan.setRepeatSearch(repeatSearch);
		mutMan.setInitEw(initEw);
		mutMan.setGamma(gamma);
		mutMan.setEpsilon(epsilon);
		mutMan.setStericE(stericE);
		mutMan.setPruningE(pruningE);
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
		mutMan.setNumInAS(numInAS);
		mutMan.numTotalRotamers(totalNumRotamers);
		mutMan.numResAllowed(numAAallowed);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setCalculateVolumes(false);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setScaleInt(scaleInt);
		mutMan.setMaxIntScale(maxIntScale);
		mutMan.setLigPartFn(q_L);
		
		if (resumeSearch)
			mutMan.setBestScore(bestScore);	// Set the current best score from the partial results
		else
			mutMan.setBestScore(new BigDecimal("0.0")); //the initial best score is 0.0
		
		try{
			handleDoMPIMaster(mutMan,mutArray.length);
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}
		
		System.out.println("DONE: K* computation");
	}
	
	/**
	 * Computes the partition function for the ligand using the rotamers from the (ligand) rotamer library
	 */
	private BigDecimal getLigPartFn(Molecule m, int numInAS, String ligType, String eMatrixNameMin){
		
		float minMatrix[][][][][][] = (float [][][][][][])readObject(eMatrixNameMin);
		
		int numRot = grl.getNumRotamers(m.strand[ligStrNum].residue[0].name);
		if (numRot==0) //ALA or GLY
			numRot = 1;
		
		ExpFunction ef = new ExpFunction();
		
		BigDecimal q_L = new BigDecimal("0.0");
		
		for (int i=0; i<numRot; i++)
			q_L = q_L.add(ef.exp(-minMatrix[numInAS][grl.getAARotamerIndex(ligType)][i][numInAS][0][0]/constRT));
		
		System.out.println("Ligand partition function (double): "+q_L.doubleValue());
		
		return q_L;
	}
	
	/** 
	 * Loads the mutation sequence list for Hybrid MinDEE-K*; computes a list if one cannot be loaded
	 */
	private OneMutation[] handleHybridKSLoadMutList (OneMutation mutArray[], String mutFileName, int numInAS,
			Molecule m, int numComb, int residueMutatable[][], ParamSet sParams,int residueMap[], String resDefault[],
			String resAllowed[], int numResAllowed, float targetVol,float volWindow){
		
		// Look for previous mutation file
		System.out.println();
		System.out.print("Looking for mutation list file ");
		mutArray = loadMutationList(mutFileName,numInAS,false);

		if (mutArray == null) {
			
			rl.loadVolFile((String)rParams.getValue("VOLFILE")); //load the rotamer volume file
			
			// Create the mutation list with estimated energies
			mutArray = new OneMutation[200000];
			RotamerSearch rs = new RotamerSearch(m, sysStrNum, ligStrNum, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst,doDihedE,doSolvationE,solvScale,softvdwMultiplier,rl,grl);

			int curNumSeq = 0;
			for(int i=0; i<numComb; i++) {
				// Reset each amino acid type
				System.out.print("Starting mutation combination " + i + " ... ");
				rs.refreshSystemStrand(); // clears allowables and does some other stuff

				boolean addWT = (new Boolean((String)sParams.getValue("ADDWT"))).booleanValue();
				for(int j=0; j<numInAS; j++) {
					if (residueMutatable[i][j] == 1) {
						setAllowablesHelper(rs, sParams, addWT, j, residueMap, resDefault);
					}
					else {
						rs.setAllowable(residueMap[j],resDefault[j]);
					}
				}
				
				// Perform simple mutation search for this set of mutatable residues
				curNumSeq = rs.simpleMasterMutationSearch(residueMap,numInAS,
					resAllowed,numResAllowed,curNumSeq,mutArray,targetVol-volWindow,
					targetVol+volWindow);
				System.out.println("finished");
			}
			
			System.out.println("Sequences remaining after volume filter "+curNumSeq);
			
			// We now have all the mutations in mutArray, collapse the mutArray
			//  to the actual number of mutations we have.
			OneMutation newArray[] = new OneMutation[curNumSeq];
			System.out.println("Allocated newArray");
			System.out.println("Initial Length of mutArray: "+mutArray.length);
			System.arraycopy(mutArray,0,newArray,0,curNumSeq);
			mutArray = newArray;
			System.out.println("Trimmed Length of mutArray: "+mutArray.length);
			
			System.out.print("Removing duplicates...");
			mutArray = removeDuplicates(mutArray);
			System.out.println("done");
			
			System.out.println(mutArray.length+" unique mutation sequences found in volume range "+(targetVol-volWindow)+" to "+(targetVol+volWindow));
			BigInteger numConfs = BigInteger.ZERO;
			for (int i=0; i<mutArray.length; i++)
				numConfs = numConfs.add(mutArray[i].numConfUB.add(mutArray[i].numConfB));
			System.out.println("Total number of conformations (bound and unbound) for all sequences: "+numConfs);
			// Save mutation list
			saveMutationList(mutArray,mutFileName,false);
		}

	
		// Sort the mutation list
		// System.out.print("Sorting mutation list ... ");
		RyanQuickSort rqs = new RyanQuickSort();
		rqs.Sort(mutArray);
		rqs = null;
		// System.out.println("done");
		
		return mutArray;
	}


	/**
	 * Reads the results of a partially completed run into an array of CommucObj. The MutationManager then queries 
	 * this array before sending out a task.
	*/
	public OneMutation[] readResumeFile(OneMutation resumeResults[], String resumeFilename, int numInAS, boolean distrDACS, boolean PEMcomp, int initDepth) {
	
		BufferedReader bufread = null;
		try {
			File file = new File(resumeFilename);		
			FileReader fr = new FileReader(file);  
			bufread = new BufferedReader(fr); 
		}
		catch (FileNotFoundException e) {
			System.out.println("ERROR: Resume File Not Found");
			return(null);
		}

		boolean done = false;
		String str = null;
		int resultNum = 0;

		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
			 System.out.println("ERROR: An error occurred while reading input");
			 System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else if (getToken(str,1).equalsIgnoreCase("Completed")) {
				if (PEMcomp) {//PEM computation
					resumeResults[resultNum].resMut = new int[numInAS];
					for(int q=0;q<numInAS;q++) 
						resumeResults[resultNum].resMut[q] = new Integer(getToken(str,8+q)).intValue();
					resumeResults[resultNum].flagMutType = getToken(str,8+numInAS);
				}
				else { //mutation search
					if (!distrDACS){ //Hybrid-K* or MinDEE/A* resume
						resumeResults[resultNum].score = new BigDecimal(getToken(str,5));
						resumeResults[resultNum].resTypes = new String[numInAS];
						for(int q=0;q<numInAS;q++) {
							resumeResults[resultNum].resTypes[q] = getToken(str,17+q);	
						}
					}
					else {//distributed DACS resume
						resumeResults[resultNum].mutNum = new Integer(getToken(str,3)).intValue();
						resumeResults[resultNum].score = new BigDecimal(getToken(str,7));
						resumeResults[resultNum].resMut = new int[initDepth];
						for(int q=0;q<initDepth;q++)
							resumeResults[resultNum].resMut[q] = new Integer(getToken(str,9+q));
					}
				}
				resultNum++;
			}
		}
					
		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}

		// Resize completed mutation array
		OneMutation temp[] = new OneMutation[resultNum];
		System.arraycopy(resumeResults,0,temp,0,resultNum);
		resumeResults = temp;
		return(resumeResults);
	}
	
	
	// Finds the index of the mutation in resumeResults with the same
	//  mutation sequence as the targetMutation. If none are found, -1
	//  is returned.
	public int findMutationIndex(OneMutation resumeResults[],
								 String targetMutation[]) {
		
		for(int q=0;q<resumeResults.length;q++) {
			if (resumeResults[q].isSame(targetMutation))
				return(q);
		}
		return(-1);
	}

	
	// Attempts to read a list of mutations from file
	public OneMutation[] loadMutationList(String fName, int numInAS, boolean PEMcomp) {
		
		BufferedReader bufread = null;
		try {
			File file = new File(fName);
			FileReader fr = new FileReader(file);
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			System.out.println(" ... no mutation list file found. Computing one.");
			return(null);
		}

		boolean done = false;
		String str = null;
		int resultNum = 0;
		OneMutation mutList[] = new OneMutation[1];
		
		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
			 System.out.println("ERROR: An error occurred while reading input");
			 System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else {
				if (PEMcomp) {//PEM computation
					mutList[resultNum] = new OneMutation();
					mutList[resultNum].resMut = new int[numInAS];
					
					for(int q=0;q<numInAS;q++) 
						mutList[resultNum].resMut[q] = new Integer(getToken(str,1+q)).intValue();
					
					mutList[resultNum].flagMutType = getToken(str,1+numInAS);
				}
				else {//mutation search
					mutList[resultNum] = new OneMutation();
					mutList[resultNum].score = new BigDecimal(getToken(str,1));
					mutList[resultNum].vol = new Float(getToken(str,2)).floatValue();
					mutList[resultNum].resTypes = new String[numInAS];
					for(int q=0;q<numInAS;q++) {
						mutList[resultNum].resTypes[q] = getToken(str,3+q);	
					}					
				}
				
				resultNum++;
				if (resultNum >= mutList.length){
					OneMutation newArray[] = new OneMutation[mutList.length+1000];
					System.arraycopy(mutList,0,newArray,0,resultNum);
					mutList = newArray;
				}
			}
		}
					
		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}

		// Resize completed mutation array
		OneMutation temp[] = new OneMutation[resultNum];
		System.arraycopy(mutList,0,temp,0,resultNum);
		mutList = temp;
		System.out.println(" ... read "+mutList.length+" mutations from mutation list "+fName);
		return(mutList);
	}	

	// Saves the list of mutations so that a PEM computation/mutation search
	//  doesn't need to recompute these during a resume. Thus
	//  the resume can go more quickly.
	public void saveMutationList(OneMutation mutList[], String fName, boolean PEMcomp) {

		if (mutList.length == 0)
			return;
		
		int numInAS = 0;
		if (PEMcomp)
			numInAS = mutList[0].resMut.length;
		else
			numInAS = mutList[0].resTypes.length;

		PrintStream printStream = setupOutputFile(fName);
		for(int q=0;q<mutList.length;q++) {
			if (PEMcomp) {//PEM computation
				for(int w=0;w<numInAS;w++) {
					printStream.print(" "+mutList[q].resMut[w]);
				}
				printStream.print(" "+mutList[q].flagMutType);
				printStream.println();
			}
			else { //mutation search
				printStream.print(mutList[q].score + " " + mutList[q].vol);
				for(int w=0;w<numInAS;w++) {
					printStream.print(" "+mutList[q].resTypes[w]);
				}
				printStream.println();
			}
		}
		printStream.close();
	}

	//Removes duplicate mutations (for which the mutation sequence is the same) from a given list
	public OneMutation [] removeDuplicates(OneMutation mutArray[]){
		
		//First, sort the list alphabetically, according to the mutation sequence
		RyanQuickSort rqs = new RyanQuickSort();
		rqs.Sort(mutArray);
		rqs = null;
		
		//Copy mutArray into nArray, excluding duplicate entries
		OneMutation nArray[] = new OneMutation[mutArray.length];
		
		//Copy the first element
		nArray[0] = mutArray[0];
		int nAIndex = 1;
		
		//Compare each mutation with the previous one in the list
		for (int i=1; i<mutArray.length; i++){ //for each mutation
			if (!(mutArray[i].isSame(mutArray[i-1].resTypes))){ //different sequence
				nArray[nAIndex] = mutArray[i];
				nAIndex++;
			}
		}
		
		mutArray = new OneMutation[nAIndex];
		System.arraycopy(nArray,0,mutArray,0,mutArray.length);
		
		return mutArray;//return the reduced list
	}

	// Mutation search Slave function
	public CommucObj handleKSSlave(CommucObj cObj) {

		if (cObj.PEMcomp){ //PEM computation
			if (!cObj.entropyComp) //PEM computation
				cObj = handleComputeAllPairwiseRotamerEnergiesSlave(cObj);
			else //entropy E matrix computation
				cObj = handleDoResEntropySlave(cObj);
		}		
		else { //distributed mutation search
			if (cObj.distrDACS){ //running distributed DACS
				cObj = doDistrDACSSlave(cObj);
			}
			else if (cObj.distrDEE){ //running distributed DEE
				cObj = doDistrDEESlave(cObj);
			}
			else { //running Hybrid MinDEE-K*
				cObj = hybridKScompute(cObj);
			}
		}
		return cObj;
	}
	
	/**
	 * Handles the computation of the K* score for a single mutation sequence with the target ligand.
	 * The 'cObj' parameter contains the mutation search input distributed by the main processor.
	 * Returns the results of the computation to the main processor.
	 */
	private CommucObj hybridKScompute(CommucObj cObj){
		
		for(int runNum = 0; runNum<2; runNum++) {
			long startTime = System.currentTimeMillis();

			// First compute q_E, then q_EL
			boolean ligPresent = (runNum == 1);
			
			ParamSet params = null;
			String minEmatrixFile = null;
			String maxEmatrixFile = null;
			if ( (!ligPresent) && ((new Boolean((String)cObj.params.getValue("USEUNBOUNDSTRUCT"))).booleanValue()) ) { //use a different input PDB structure for the unbound case
				params = new ParamSet();
				params.setParamsValues(cObj.params.getParams(), cObj.params.getValues(), cObj.params.getCurNum());
				params.setValue("PDBNAME",params.getValue("UNBOUNDPDBNAME"));
				params.setValue("PDBLIGNUM","-1");
				minEmatrixFile = params.getValue("MINENERGYMATRIXNAMEUNBOUND")+".dat";
				maxEmatrixFile = params.getValue("MAXENERGYMATRIXNAMEUNBOUND")+".dat";
			}
			else { //a single input PDB structure is used for the bound and unbound computations
				params = cObj.params;
				minEmatrixFile = cObj.arpFilenameMin;
				maxEmatrixFile = cObj.arpFilenameMax;
			}

			//Setup the molecule system
			Molecule m = new Molecule();
			int numLigRotamers = setupMolSystem(m,params,ligPresent,cObj.ligType);

			RotamerSearch rs = new RotamerSearch(m,sysStrNum,ligStrNum, hElect, hVDW, hSteric, true,
						true, cObj.epsilon, cObj.stericThresh, cObj.softStericThresh, cObj.distDepDielect, cObj.dielectConst,cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult,rl,grl);		
			
			System.out.print("Loading precomputed min energy matrix...");
			rs.loadPairwiseEnergyMatrices(minEmatrixFile,true);
			System.out.println("done");
			
			if (cObj.doMinimization){ //MinDEE, so load the max matrix
				System.out.print("MinDEE: Loading precomputed max energy matrix...");
				rs.loadPairwiseEnergyMatrices(maxEmatrixFile,false);
				System.out.println("done");
			}

			System.out.println("Beginning setAllowables");
			// Ligand allowable set in the RotamerSearch() constructor
			for(int q=0; q<cObj.numInAS; q++)
				rs.setAllowable(cObj.residueMap[q],cObj.currentMutation[q]);
			
			
			//Initially, no rotamers have been pruned
			boolean prunedRotAtRes[] = new boolean[cObj.numInAS*cObj.numTotalRotamers+numLigRotamers];
			for (int i=0; i<prunedRotAtRes.length; i++)
				prunedRotAtRes[i] = false;
			
			
			//Prune all rotamers that are incompatible with the template (intra E + re-to-template E >= stericE)
			prunedRotAtRes = rs.DoPruneStericTemplate(cObj.numInAS, cObj.numTotalRotamers, numLigRotamers, 
					cObj.residueMap, rotamerIndexOffset, prunedRotAtRes, cObj.stericE);			
		
			//Perform the DEE pruning
			if (cObj.doMinimization) //compute the MinDEE interval terms
				rs.doCompMinDEEIntervals(cObj.numInAS, cObj.numTotalRotamers, numLigRotamers, cObj.residueMap, 
						rotamerIndexOffset, prunedRotAtRes, cObj.scaleInt, cObj.maxIntScale);
			
			prunedRotAtRes = rs.DoDEEGoldstein(cObj.numInAS, cObj.numTotalRotamers, numLigRotamers, 
					cObj.residueMap,cObj.rotamerIndexOffset, cObj.initEw, prunedRotAtRes, cObj.doMinimization, false, cObj.minimizeBB);
			
			//Prune with MinBounds (last parameter is false)
			prunedRotAtRes = rs.DoMinBounds(cObj.numInAS,cObj.numTotalRotamers,numLigRotamers,
					cObj.residueMap,cObj.rotamerIndexOffset,cObj.pruningE,prunedRotAtRes,cObj.initEw, false, false);
			
			//Compute the Ec value and prunedIsSteric[] (last parameter is true)
			rs.DoMinBounds(cObj.numInAS,cObj.numTotalRotamers,numLigRotamers,
					cObj.residueMap,cObj.rotamerIndexOffset,cObj.pruningE,prunedRotAtRes,cObj.initEw, false, true);
		
			
			boolean usingInitialBest = ligPresent;
			BigDecimal initialBest = (new BigDecimal("0.0"));
			
			if (usingInitialBest)
				initialBest =  cObj.q_E.multiply(cObj.bestScore.multiply(cObj.q_L)).multiply(new BigDecimal(cObj.gamma * cObj.epsilon));

			rs.slaveDoRotamerSearch(cObj.computeEVEnergy,cObj.doMinimization,cObj.numInAS,
				cObj.numResAllowed,cObj.numTotalRotamers,cObj.rotamerIndexOffset,resAllowed,
				cObj.residueMap,usingInitialBest,initialBest,cObj,cObj.minimizeBB,false,null,cObj.doBackrubs,cObj.backrubFile);
			
			if ((cObj.repeatSearch)&&(rs.repeatSearch)){ //the desired accuracy was not achieved, so repeat the search: the setup is already done
				
				rs.repeatSearch = false; //reset the flag
				if (ligPresent)
					cObj.EL_repeatEw = true; //set the flag
				else
					cObj.E_repeatEw = true;
				
				rs.slaveDoRotamerSearch(cObj.computeEVEnergy,cObj.doMinimization,cObj.numInAS,
						cObj.numResAllowed,cObj.numTotalRotamers,cObj.rotamerIndexOffset,resAllowed,
						cObj.residueMap,usingInitialBest,initialBest,cObj,cObj.minimizeBB,false,null,cObj.doBackrubs,cObj.backrubFile);
			}
				
			long stopTime = System.currentTimeMillis();
			if(runNum == 0)
				cObj.q_E_Time = Math.round((stopTime - startTime) / 1000.0f);
			else
				cObj.q_EL_Time = Math.round((stopTime - startTime) / 1000.0f);
		} // end for(runNum)
		
		return cObj;
	}
	
	//Load the pairwise energy matrices; if not computed, compute, and the load
	private void loadPairwiseEnergyMatrices(ParamSet sParams, RotamerSearch rs, String minMatrixFile, boolean doMinimize, String maxMatrixFile){
		
		rs.loadPairwiseEnergyMatrices(minMatrixFile,true);
		if (doMinimize)
			rs.loadPairwiseEnergyMatrices(maxMatrixFile,false);
		
		if ( (rs.getMinMatrix()==null) || ( doMinimize && rs.getMaxMatrix()==null ) ) { //at least one of the matrices not computed, so compute
			
			System.out.println("Precomputed energy matrices not available..");
			
			long startTime = System.currentTimeMillis();
			
			handleComputeAllPairwiseRotamerEnergiesMaster(sParams);
			
			long stopTime = System.currentTimeMillis();
			System.out.println("PEM execution time: "+((stopTime-startTime)/(60.0*1000.0)));
			System.out.println("DONE: Pairwise energy matrix precomputation");
			
			rs.loadPairwiseEnergyMatrices(minMatrixFile,true);
			if (doMinimize)
				rs.loadPairwiseEnergyMatrices(maxMatrixFile,false);
		}
	}
	
/////////////////////////////////////////////////////////////////////////
// MIN and MAX pairwise energy matrices computation
/////////////////////////////////////////////////////////////////////////
	
	/**
	 * This function sets up the computation for all min and max pairwise rotamer interaction energies and stores 
	 * these energies into user-specified precomputed energy matrices. The computation is distributed to the available processors.
	 * For each rotamer, for each allowed amino acid type, for each flexible residue position, the following energies
	 * are computed: intra-rotamer, rotamer-to-template, rotamer-rotamer, and template-only energies. If a ligand is
	 * present, all energies involving the ligand are also computed.
	 */
	public void handleComputeAllPairwiseRotamerEnergiesMaster(ParamSet sParams) {	
		
		int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
		int numMutations = 2; //pairwise energies are computed
		boolean doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE"))).booleanValue();
		boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB"))).booleanValue();
		boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS"))).booleanValue();
		String backrubFile = (String)sParams.getValue("BACKRUBFILE");
		String runName = (String)sParams.getValue("RUNNAME");
		String minEMatrixName = (String)sParams.getValue("MINENERGYMATRIXNAME");
		String maxEMatrixName = (String)sParams.getValue("MAXENERGYMATRIXNAME");
		
		boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = (String)sParams.getValue("LIGTYPE");
		
		if (!doMinimize) //no minimization
			minimizeBB = false;
		if (!minimizeBB) //not backbone minimization
			doBackrubs = false;
		
		System.out.println("Run Name: "+runName);
		System.out.println("Precomputed Minimum Energy Matrix: "+minEMatrixName);
		System.out.println("Precomputed Maximum Energy Matrix: "+maxEMatrixName);
		System.out.println("Ligand Type: "+ligType);
		System.out.println("Num Residues Allowed to Mutate: "+numMutations);

		
		System.out.println("Computing _All_ Rotamer-Rotamer Energies");
		
		System.out.println("Starting minimum and maximum bound energy computation");
		
		//Setup the molecule system
		Molecule m = new Molecule();
		int numLigRotamers = setupMolSystem(m,sParams,ligPresent,ligType);
					
		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.getValue("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
			residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
		}
		System.out.println();
		
		RotamerSearch rs = new RotamerSearch(m,sysStrNum,ligStrNum, hElect, hVDW, hSteric, true,
				true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, rl, grl);
		
		int resMut[] = new int[numInAS];
		for (int i=0; i<resMut.length; i++)
			resMut[i] = 1;
		
		System.out.println("Beginning setAllowables");
		//Ligand allowable set in the RotamerSearch() constructor	
		//Set the allowable AAs for each AS residue
		boolean addWT = (new Boolean((String)sParams.getValue("ADDWT"))).booleanValue();
		for (int j=0; j<numInAS; j++){
			setAllowablesHelper(rs, sParams, addWT, j, residueMap, resDefault);
		}
		
		//initialize the pairwise energy matrices (full initialization - for all residues in residueMap[], the ligand, and the template)
		PEMHandler pemH = new PEMHandler();
		float mutationEnergiesMin[][][][][][] = pemH.initializePairEMatrix(numInAS,ligPresent,resMut,residueMap,rs,ligPresent,ligType,true,true,rl,grl,numAAallowed,true);
		float mutationEnergiesMax[][][][][][] = pemH.copyMultiDimArray(mutationEnergiesMin);
		
		
		OneMutation mutArray[] = getMutArrayPairEcomp(numInAS,ligPresent,minimizeBB);

		/*//Sort the mutation list
		System.out.print("Sorting mutation list ... ");
		RyanQuickSort rqs = new RyanQuickSort();
		rqs.Sort(mutArray);
		rqs = null;
		System.out.println("done");*/
		
		MutationManager mutMan = new MutationManager(null,mutArray,true);
		mutMan.setResidueMap(residueMap);
		mutMan.setResDefault(resDefault);
		mutMan.setRotamerIndexOffset(rotamerIndexOffset);
		mutMan.setLigType(ligType);
		mutMan.setarpFilenameMin(minEMatrixName);
		mutMan.setPairEMatrixMin(mutationEnergiesMin);
		if (doMinimize){
			mutMan.setarpFilenameMax(maxEMatrixName);
			mutMan.setPairEMatrixMax(mutationEnergiesMax);
		}
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
		mutMan.setNumInAS(numInAS);
		mutMan.setnumLigRotamers(numLigRotamers);
		mutMan.numTotalRotamers(totalNumRotamers);
		mutMan.numResAllowed(numAAallowed);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setDoBackrubs(doBackrubs);
		mutMan.setBackrubFile(backrubFile);
		mutMan.setCalculateVolumes(false);
		mutMan.setLigPresent(ligPresent);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);

		try{
			handleDoMPIMaster(mutMan,mutArray.length);
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}
		
		outputObject(mutMan.getMinEmatrix(),minEMatrixName+".dat");
		if (doMinimize)
			outputObject(mutMan.getMaxEmatrix(),maxEMatrixName+".dat");
		
		System.out.println("DONE: Pairwise energy matrix precomputation..");
	}
	
	/**
	 * Computes a specific part of the pairwise energy matrices, as specified by the parameters in the 'cObj' parameter,
	 * distributed by the main processor. Returns the results of the computation to the main processor.
	 */
	public CommucObj handleComputeAllPairwiseRotamerEnergiesSlave(CommucObj cObj) {
				
		long startTime = System.currentTimeMillis();
						
		boolean ligPresent = cObj.ligPresent;

		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,cObj.params,ligPresent,cObj.ligType);
		
		RotamerSearch rs = new RotamerSearch(m,sysStrNum,ligStrNum, hElect, hVDW, hSteric, true,
					true, cObj.epsilon, cObj.stericThresh, cObj.softStericThresh, cObj.distDepDielect, cObj.dielectConst, cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult,rl,grl);
	
		boolean useLig = ligPresent;
		if (((cObj.flagMutType.compareTo("AS-AS")==0)||(cObj.flagMutType.compareTo("SHL-AS")==0)||(cObj.flagMutType.compareTo("TEMPL")==0))&&(ligPresent)){
			useLig = false;
			m.deleteStrand(ligStrNum);//we do not need the ligand for these runs
			rs = new RotamerSearch(m,sysStrNum,-1, hElect, hVDW, hSteric, true,
					true, cObj.epsilon, cObj.stericThresh, cObj.softStericThresh, cObj.distDepDielect, cObj.dielectConst, cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult,rl,grl);					
		}
		
		System.out.println("Beginning setAllowables");
		// Ligand allowable set in the RotamerSearch() constructor
		boolean addWT = (new Boolean((String)cObj.params.getValue("ADDWT"))).booleanValue();
		for(int j=0; j<cObj.numInAS; j++) {
			if (cObj.resMut[j] == 1) {
				setAllowablesHelper(rs, cObj.params, addWT, j, cObj.residueMap, cObj.resDefault);
			}
		}
		
		
		boolean shellRun = false;
		boolean intraRun = false;
		boolean templateOnly = false;		
		
		if (cObj.flagMutType.compareTo("TEMPL")==0){
			
			// **** Normal Active Site residue runs ****
			// Computes active site residue to active site residue pair energies					
			shellRun = true;ligPresent = false;intraRun = false;templateOnly = true;
		}
		else if (cObj.flagMutType.compareTo("AS-AS")==0){
			
			// **** Normal Active Site residue runs ****
			// Computes active site residue to active site residue pair energies					
			shellRun = false;ligPresent = false;intraRun = false;templateOnly = false;
		}	
		else if (cObj.flagMutType.compareTo("SHL-AS")==0){
			
			// Then shell runs for the active site residues
			// Computes the active site residue rotamers to shell energies					
			shellRun = true;ligPresent = false;intraRun = false;templateOnly = false;				
		}
		else if (cObj.flagMutType.compareTo("INTRA")==0){
			
			// Compute all intra-residue energies					
			shellRun = false;intraRun = true;templateOnly = false;
		}				
		else if (cObj.flagMutType.compareTo("LIG-AS")==0){
			
			// **** Ligand present runs ****
			// This section computes the inter-residue energies between
			//  active site residues and the ligand
			shellRun = false;intraRun = false;templateOnly = false;			
		}
		else { //(cObj.flagMutType.compareTo("LIG-SHL")==0)
			
			// Computes ligand rotamer to shell energies
			shellRun = true; intraRun = false;templateOnly = false;
		}
		
		// The goal is that the total energy of a system can be bounded by the sum of 
		//  all pairwise active site residue entries plus the entry for each active site
		//  residue's shell run plus each active site residue's self intra-residue energy.
		//  If a ligand is present then one should add the ligand to shell energy, the
		//  ligand to each active site residue pairwise energy, and the ligand self intra-
		//  residue energy.
		
		//initialize the pairwise energy matrices (partial initialization - only for the residues involved in this computation, e.g., AS-AS)
		PEMHandler pemH = new PEMHandler();
		float minEmatrix[][][][][][] = pemH.initializePairEMatrix(cObj.numInAS,cObj.ligPresent,cObj.resMut,cObj.residueMap,
				rs,useLig,cObj.ligType,shellRun,intraRun,rl,grl,numAAallowed,false);
		float maxEmatrix[][][][][][] = pemH.copyMultiDimArray(minEmatrix);
		
		
		//Compute the corresponding matrix entries
		rs.simplePairwiseMutationAllRotamerSearch(cObj.residueMap,cObj.numInAS,cObj.rotamerIndexOffset,
				cObj.numTotalRotamers,cObj.doMinimization,ligPresent,shellRun,intraRun,
				cObj.resMut,minEmatrix,maxEmatrix,cObj.minimizeBB,cObj.doBackrubs,
				templateOnly,cObj.backrubFile);
		
		
		long stopTime = System.currentTimeMillis();
		cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);
		
		
		//Store the information in less space to allow the master node to buffer several cObj at once
		cObj.compEE = pemH.generateCompEE(minEmatrix, maxEmatrix);
		
		return cObj;
	}
	
	//Generates and saves to file the mutation list for the pairwise energy matrix computation
	private OneMutation[] getMutArrayPairEcomp(int numInAS, boolean ligPresent, boolean minimizeBB){
		
		final int numMutations = 2; //pairwise energy computation
		
		// Generate all combinations
		int numComb = factorial(numInAS).divide(factorial(numInAS-numMutations).multiply(factorial(numMutations))).intValue();
		int residueMutatable[][] = new int[numComb][numInAS];
		generateCombinations(residueMutatable,numInAS,numMutations);
		// At this point each row of residueMutatble is a 0/1 array which specifies a mutation 
		//  pair combination, 1 indicates that that residue can mutate in the specified combination
		
		System.out.println("Number of possible mutation combinations: "+numComb);
		
		// Create the mutation list with estimated energies
		OneMutation mutArray[] = new OneMutation[numComb];			
		
		//Set the AS-AS mutations
		int curMutNum = 0;
		for(int i=0; i<numComb; i++) {
			
			mutArray[i] = new OneMutation();
			mutArray[i].flagMutType = "AS-AS";
			mutArray[i].resMut = new int[numInAS];
			for(int j=0; j<numInAS; j++) {
				mutArray[i].resMut[j] = residueMutatable[i][j];
			}
			
			// Perform simple mutation search for this set of mutatable residues
			curMutNum++;
		}
		
		//Add the runs for template only, AS-shell, AS-ligand, intra-residue energies, and ligand-shell
		int t=0;
		if (minimizeBB)
			t = 1;
		int numOtherMut;
		if (ligPresent)
			numOtherMut = 2+t+2*numInAS;
		else
			numOtherMut = 1+t+numInAS;
		OneMutation otherMutArray[] = new OneMutation[numOtherMut];
		
		for (int i=0; i<numOtherMut; i++){
			otherMutArray[i] = new OneMutation();
			otherMutArray[i].resMut = new int[numInAS];
		}
		
		//Set the AS-shell mutations
		for (int i=0; i<numInAS; i++){
			otherMutArray[i].flagMutType = "SHL-AS";
			for (int j=0; j<numInAS; j++){
				if (i==j)
					otherMutArray[i].resMut[j] = 1;
				else
					otherMutArray[i].resMut[j] = 0;
			}
		}
		
		//Set the intra-residue energies run
		otherMutArray[numInAS].flagMutType = "INTRA";
		for (int j=0; j<numInAS; j++)
			otherMutArray[numInAS].resMut[j] = 1;
		
		//Set the template energy run
		if (minimizeBB){
			otherMutArray[numInAS+1].flagMutType = "TEMPL";
			for (int j=0; j<numInAS; j++)
				otherMutArray[numInAS+1].resMut[j] = 0;
		}
		
		if (ligPresent){//if the ligand is present, set the corresponding runs
			
			//Set the AS-ligand mutations
			for (int i=1+t+numInAS; i<=2*numInAS+t; i++){
				otherMutArray[i].flagMutType = "LIG-AS";
				for (int j=0; j<numInAS; j++){
					if ((i-1-t-numInAS)==j)
						otherMutArray[i].resMut[j] = 1;
					else
						otherMutArray[i].resMut[j] = 0;
				}
			}
			
			//Set the ligand-shell run
			otherMutArray[1+t+2*numInAS].flagMutType = "LIG-SHL";
			for (int j=0; j<numInAS; j++)
				otherMutArray[1+t+2*numInAS].resMut[j] = 0;
		}
		
		// We now have all the mutations in mutArray, collapse the mutArray
		//  to the actual number of mutations we have.
		OneMutation newArray[] = new OneMutation[curMutNum+numOtherMut];
		System.arraycopy(otherMutArray,0,newArray,0,numOtherMut);//add the other mutations first
		System.arraycopy(mutArray,0,newArray,numOtherMut,curMutNum);//then add the AS-AS mutations
		
		mutArray = newArray;
		System.out.println("Length of mutArray: "+mutArray.length);
		
		return mutArray;
	}
	
	// Finds the index of the mutation in resumeResults with the same
	//  mutation sequence as resMut. If none are found, -1 is returned.
	public int sampFindMutationIndex(OneMutation resumeResults[], String flMutType, int mutResidues[]) {
		
		for(int q=0;q<resumeResults.length;q++) 
			if ((resumeResults[q].flagMutType.compareTo(flMutType)==0) && (sameSeq(resumeResults[q].resMut,mutResidues)))
				return(q);
			
		return(-1);
	}

	//Determines if the residues that can mutate are the same for two mutation sequences
	private boolean sameSeq (int computedResMut[], int allResMut[]){
		
		boolean found = true;
		for (int i=0; i<computedResMut.length; i++){
			if (computedResMut[i]!=allResMut[i])
				found = false;
		}
		return found;
	}
///////////////////////////////////////////////////////////////////////////
//	End of MIN and MAX Pairwise Energy Precomputation
///////////////////////////////////////////////////////////////////////////
	
////////////////////////////////////////////////////////////////
//	 Compute minimized-GMEC section
////////////////////////////////////////////////////////////////	
	/** 
	 * Computes the (energy-minimized) structure for the rotameric conformations specified by the input file.
	 * This function is used to generate structures for a set of output conformations from a DEE/A* search.
	 */
	public void handleMinDEEApplyRot(String s) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutation search parameter filename (string)

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

		// Pull search parameters
		int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
		String confResFile = (String)sParams.getValue("CONFRESFILE");
		int numResults = (new Integer((String)sParams.getValue("NUMRESULTS"))).intValue();
		boolean doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE"))).booleanValue();
		boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB"))).booleanValue();
		boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS"))).booleanValue();
		String backrubFile = (String)sParams.getValue("BACKRUBFILE");
		boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
		String ligType = (String)sParams.getValue("LIGTYPE");
		boolean outputPDB = (new Boolean((String)sParams.getValue("OUTPUTPDBS"))).booleanValue();
		
		int numResidues;
		if (ligPresent)
			numResidues = numInAS+1;
		else
			numResidues = numInAS;
		
		//Read the results file into the AA and rot matrices
		String AAtypes[] = new String[numResults*numResidues];
		int rotNums[] = new int[numResults*numResidues];
		readRotResFile(confResFile,AAtypes,rotNums,numResults,numResidues);
		
		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,ligPresent,ligType);

		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.getValue("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
			residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
		}
		System.out.println();
		
		String curSeq[] = new String[numInAS];
		
		int numSaved = 0;
		for (int curResult=0; curResult<numResults; curResult++){
			System.out.print("Starting minimization of result "+(curResult+1)+"..");
			
			Amber96ext a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, softvdwMultiplier);
			
			StrandRotamers sysLR = null;
			StrandRotamers ligLR2 = null;
			
			sysLR = new StrandRotamers(rl,m.strand[sysStrNum]);
			if(ligPresent)
				ligLR2 = new StrandRotamers(grl,m.strand[ligStrNum]);
			
			//the starting index for the current result within AAtypes[] and rotNums[]
			int startInd = numResidues*curResult;
		
			//Set the allowables for the AS residues
			for(int j=0;j<numInAS;j++){
				curSeq[j] = AAtypes[startInd+j];
				sysLR.setAllowable(residueMap[j],curSeq[j]);
				sysLR.changeResidueType(m,residueMap[j],curSeq[j],true,true);
				m.strand[sysStrNum].residue[residueMap[j]].flexible = true;
			}
			if(ligPresent){
				ligLR2.setAllowable(0,ligType); //the only allowable AA type for the ligand
				if (m.strand[ligStrNum].isProtein) //apply mutation
					ligLR2.changeResidueType(m,0,ligType,true,true);
				m.strand[ligStrNum].residue[0].flexible = true;
			}
			
			a96ff.calculateTypesWithTemplates();
			a96ff.initializeCalculation();
			a96ff.setNBEval(hElect,hVDW);
			a96ff.setLigandNum((new Integer((String)sParams.getValue("PDBLIGNUM"))).intValue());
			
			int curAA[] = new int[m.numberOfResidues];
			for(int j=0;j<m.strand[sysStrNum].numberOfResidues;j++){
				curAA[j] = sysLR.getIndexOfNthAllowable(j,0);
			}
			
			SimpleMinimizer simpMin = null;
			BBMinimizer bbMin = null;
			BackrubMinimizer brMin = null;
			if ( doMinimize && (!minimizeBB) ){
				simpMin = new SimpleMinimizer();
				if(ligPresent)
					simpMin.initialize(m,0,1,a96ff,sysLR,ligLR2,curAA,ligLR2.getIndexOfNthAllowable(0,0),doDihedE,rl,grl);
				else
					simpMin.initialize(m,0,a96ff,sysLR,curAA,doDihedE,rl);
			}
			else if (minimizeBB) {
				if (!doBackrubs){
					bbMin = new BBMinimizer();
					if (ligPresent)
						bbMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum);
					else
						bbMin.initialize(m, a96ff, residueMap, sysStrNum);
				}
				else {
					brMin = new BackrubMinimizer();
					if (ligPresent)
						brMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum, backrubFile, hSteric, stericThresh);
					else
						brMin.initialize(m, a96ff, residueMap, sysStrNum, backrubFile, hSteric, stericThresh);
				}
			}
			
			m.backupAtomCoord();
			//Apply the corresponding rotamers
			for (int j=0; j<numInAS; j++){										
				if (rl.getNumRotForAAtype(curAA[residueMap[j]])!=0){//not GLY or ALA
					int curRot = rotNums[startInd+j];
					sysLR.applyRotamer(m, residueMap[j], curRot);
				}
			}				
			if (ligPresent){ //apply the ligand rotamer
				if (grl.getNumRotForAAtype(ligLR2.getIndexOfNthAllowable(0,0))!=0){//not GLY or ALA
					int curRot = rotNums[startInd+numInAS]; //the ligand rotamer
					ligLR2.applyRotamer(m, 0, curRot);//the ligand level
				}
			}	
			
			double unMinE[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //the energy before minimization
			
			if ( doMinimize && (!minimizeBB) )
				simpMin.minimize(35,false);
			else if (minimizeBB) {
				if (!doBackrubs)
					bbMin.minimizeFull(false);
				else
					brMin.minimizeFull();
			}
			
			double minE[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //the energy after minimization
			if ((doMinimize)&&(!minimizeBB)&&(doDihedE)) //add dihedral energies
				minE[0] += simpMin.computeDihedEnergy();
			
			if (outputPDB){ //save molecule
				saveMolecule(m,"pdbs/savedMol"+(numSaved+1)+".pdb",(float)minE[0]);
				numSaved++;
			}
			
			m.restoreAtomCoord();
			m.updateCoordinates();	
			
			if (minE[0]>unMinE[0])
				minE = unMinE;
			
			System.out.println("done");
		}
		System.out.println("done");
	}
	
	private void readRotResFile (String resFile, String AAtypes[], int rotNums[], int numResults, int numResidues){
		
		BufferedReader bufread = null;
		try {
			File file = new File(resFile);		
			FileReader fr = new FileReader(file);  
			bufread = new BufferedReader(fr); 
		}
		catch (FileNotFoundException e) {
			System.out.println("ERROR: Results File Not Found");
		}

		boolean done = false;
		String str = null;
		int resultNum = 0;

		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
			 System.out.println("ERROR: An error occurred while reading input");
			 System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else {
				for (int i=0; i<numResidues; i++){
					AAtypes[resultNum*numResidues+i] = getToken(str,2+i);
					rotNums[resultNum*numResidues+i] = new Integer((String)getToken(str,2+numResidues+i)).intValue();
				}
				resultNum++;
			}
		}
		
		if (numResults!=resultNum){
			System.out.println("Error: Not all results available for reading");
			System.exit(0);
		}
					
		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}
	}
////////////////////////////////////////////////////////////////
//	 End of Compute minimized-GMEC section
////////////////////////////////////////////////////////////////

	
	
	
///////////////////////////////////////////////////////////////////////////
//	DEE section
///////////////////////////////////////////////////////////////////////////
	/**
	 * Performs a DEE (Traditional, MinDEE, BD, or BRDEE) pruning with A* search, with or without DACS;
	 * the only parameter 's' (String) includes the command-line arguments specifying the filenames of the two input configuration files.
	 * If distributed DACS is performed, the computation is distributed to the available processors for evaluation.
	*/
	public void handleDoDEE(String s) {
		
		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: DEE config filename (string)
		
		//long startTimeAll = System.currentTimeMillis(); //line commented by Swati
		long startTimeAllCPU = CPUTime(); // next three line added to calculate the exact time - Swati
		long startTimeAllUser = UserTime();
		long startTimeAllSystem = SystemTime();

		System.out.println("Performing DEE");

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters
		
		// Pull search parameters
		int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
		String runName = ((String)sParams.getValue("RUNNAME"));
		int numMaxMut = (new Integer((String)sParams.getValue("NUMMAXMUT"))).intValue();
		int algOption = (new Integer((String)sParams.getValue("ALGOPTION"))).intValue();
		boolean doDACS = (new Boolean((String)sParams.getValue("DODACS"))).booleanValue();
		boolean useFlags = (new Boolean((String)sParams.getValue("SPLITFLAGS"))).booleanValue();
		boolean distrDACS = (new Boolean((String)sParams.getValue("DISTRDACS"))).booleanValue();
		//boolean distrDEE = (new Boolean((String)sParams.getValue("DISTRDEE"))).booleanValue();
		boolean distrDEE = false; //the distributed DEE section is outdated and must be carefully checked before called
		boolean doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE"))).booleanValue();
		boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB"))).booleanValue();
		boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS"))).booleanValue();
		String backrubFile = ((String)sParams.getValue("BACKRUBFILE"));
		boolean approxMinGMEC = (new Boolean((String)sParams.getValue("APPROXMINGMEC"))).booleanValue();
		boolean preprocPairs = (new Boolean((String)sParams.getValue("PREPROCPAIRS"))).booleanValue();
		boolean scaleInt = (new Boolean((String)sParams.getValue("SCALEINT"))).booleanValue();
		String runNameEMatrixMin = (String)(sParams.getValue("MINENERGYMATRIXNAME"));
		String runNameEMatrixMax = (String)(sParams.getValue("MAXENERGYMATRIXNAME"));
		float initEw = (new Float((String)sParams.getValue("INITEW"))).floatValue();
		float lambda = (new Float((String)sParams.getValue("LAMBDA"))).floatValue();
		float pruningE = (new Float((String)sParams.getValue("PRUNINGE"))).floatValue();
		double stericE = (new Double((String)sParams.getValue("STERICE"))).doubleValue();
		float pairSt = (new Float((String)sParams.getValue("PAIRST"))).floatValue();
		float maxIntScale = (new Float((String)sParams.getValue("MAXINTSCALE"))).floatValue();
		double minRatioDiff = (new Double((String)sParams.getValue("MINRATIODIFF"))).doubleValue();
		int initDepth = (new Integer((String)sParams.getValue("INITDEPTH"))).intValue();
		int subDepth = (new Integer((String)sParams.getValue("SUBDEPTH"))).intValue();
		int diffFact = (new Integer((String)sParams.getValue("DIFFFACT"))).intValue();	
		boolean genInteractionGraph = (new Boolean((String)sParams.getValue("GENINTERACTIONGRAPH"))).booleanValue();
		float distCutoff = (new Float((String)sParams.getValue("DISTCUTOFF"))).floatValue();
		float eInteractionCutoff = (new Float((String)sParams.getValue("EINTERACTIONCUTOFF"))).floatValue();
		String outputConfInfo = (String)(sParams.getValue("OUTPUTCONFINFO"));
		String outputPruneInfo = (String)(sParams.getValue("OUTPUTPRUNEINFO"));
		boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = (String)(sParams.getValue("LIGTYPE"));
		boolean useEref = (new Boolean((String)sParams.getValue("USEEREF"))).booleanValue();
		boolean resumeSearch = (new Boolean((String)sParams.getValue("RESUMESEARCH"))).booleanValue();
		String resumeFilename = ((String)sParams.getValue("RESUMEFILENAME"));
		
		
		if ((!mpiRun)&&((distrDACS)||distrDEE)){
			System.out.println("ERROR: Distributed computation requires MPI");
			System.exit(1);
		}
		
		if (!doMinimize) //no minimization
			minimizeBB = false;
		if (!minimizeBB) //not backbone minimization
			doBackrubs = false;
		
		if (genInteractionGraph) //DACS is not performed when generating the interaction graph
			doDACS = false;
			
		//Setup the molecule system
		Molecule m = new Molecule();
		int numLigRotamers = setupMolSystem(m,sParams,ligPresent,ligType);
					
		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.getValue("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
			residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
		}
		System.out.println();
		
		RotamerSearch rs = new RotamerSearch(m,sysStrNum,ligStrNum, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, rl, grl);


		System.out.print("Loading precomputed energy matrix...");
		loadPairwiseEnergyMatrices(sParams,rs,runNameEMatrixMin+".dat",doMinimize,runNameEMatrixMax+".dat");
		System.out.println("done");
		
		
		/////////////////////////////////////////////////////////////
		// DEE section
		
		//long startTime = System.currentTimeMillis(); //line commented by Swati
		long startTimeCPU = CPUTime(); //three lines added for calculating exact time - Swati
		long startTimeUser = UserTime();
		long startTimeSystem = SystemTime();
		
		//Set the allowable AAs for each AS residue
		boolean addWT = (new Boolean((String)sParams.getValue("ADDWT"))).booleanValue();
		for (int j=0; j<numInAS; j++){
			setAllowablesHelper(rs, sParams, addWT, j, residueMap, resDefault);
		}
		
		float eRef[] = null;
		if (useEref) { //add the reference energies to the min (and max) intra-energies
			eRef = getResEntropyEmatricesEref(useEref,rs.getMinMatrix(),rs.sysLR,residueMap,null,numInAS);
			rs.addEref(eRef, doMinimize, ligPresent, numInAS, residueMap);
		}
		
		boolean prunedRotAtRes[] = new boolean [numInAS*totalNumRotamers + numLigRotamers];
		
		final String rotFile = ("rot_out"+System.currentTimeMillis()); //output the pruned rotamers (for distributed DEE)
		final String sfFile = ("sf_matrix"+System.currentTimeMillis()); //output the split flags (for distributed DEE)
		
		
		//first prune all rotamers that are incompatible with the template (intra E + re-to-template E >= stericE)
		System.out.println();
		System.out.println("Pruning all rotamers incompatible with the template..");			
		prunedRotAtRes = rs.DoPruneStericTemplate(numInAS, totalNumRotamers, numLigRotamers, residueMap,
			rotamerIndexOffset, prunedRotAtRes, stericE);
		System.out.println();
		
		//preprocess pairs of rotamers (mark pairwise energies greater than the cutoff as steric clashes)
		if (preprocPairs){
			System.out.println("Preprocessing pairs of rotamers, cutoff of "+pairSt);
			rs.preprocessPairs(pairSt, numInAS, ligPresent, residueMap);
			System.out.println();
		}
		
		
		//Setup and do the DEE pruning
		if ((useFlags)||(algOption>=3))
			rs.setSplitFlags(prunedRotAtRes.length);//initialize the split flags
		
		int numPrunedRot = 0;
		int numPrunedPairs = 0;
		int numPrunedRotThisRun = 0;
		int numPrunedPairsThisRun = 0;
		boolean done = false;
		int numRuns = 1;
		
		while (!done){ //repeat the pruning cycle until no more rotamers are pruned	
			
			numPrunedRotThisRun = 0; 
			numPrunedPairsThisRun = 0;
			
			System.out.println("Starting DEE cycle run: "+numRuns);
			
			if (doMinimize) //precompute the interval terms in the MinDEE criterion
				rs.doCompMinDEEIntervals(numInAS, totalNumRotamers, numLigRotamers, residueMap, rotamerIndexOffset, 
						prunedRotAtRes, scaleInt, maxIntScale);
		
			//Depending on the chosen algorithm option, apply the corresponding pruning criteria;			
			System.out.println("Starting pruning with DEE (simple Goldstein)");		
			prunedRotAtRes = rs.DoDEEGoldstein(numInAS, totalNumRotamers, numLigRotamers, residueMap,
					rotamerIndexOffset, initEw, prunedRotAtRes, doMinimize, useFlags, minimizeBB);			
			System.out.println();
			
			if ((algOption>=3)){ //simple Goldstein pairs
				System.out.println("Starting pruning with DEE (mb pairs)");			
				rs.DoDEEPairs(numInAS, totalNumRotamers, numLigRotamers, residueMap,
						rotamerIndexOffset, initEw, prunedRotAtRes, null, doMinimize, useFlags, true, false, minimizeBB, scaleInt, maxIntScale);
				System.out.println();
			}
			
			if ((useFlags)||(algOption>=3)){
				System.out.println("Starting pruning with Bounding Flags");			
				rs.DoBoundFlags(numInAS, totalNumRotamers, numLigRotamers, residueMap,
					rotamerIndexOffset, pruningE, prunedRotAtRes, initEw, useFlags);
				System.out.println();
			}
			
			System.out.println("Starting pruning with DEE (1-sp split-DEE)");			
			prunedRotAtRes = rs.DoDEEConfSplitting(numInAS, totalNumRotamers, numLigRotamers, residueMap,
					rotamerIndexOffset, initEw, prunedRotAtRes, null, doMinimize, useFlags, 1, false, minimizeBB);
			System.out.println();
			
			System.out.println("Starting pruning with Bounds");			
			prunedRotAtRes = rs.DoMinBounds(numInAS, totalNumRotamers, numLigRotamers, residueMap,
				rotamerIndexOffset, pruningE, prunedRotAtRes, initEw, useFlags, false);
			System.out.println();
			
			//check how many rotamers/pairs are pruned this run
			int numTotalPrunedRot = countPrunedRot(prunedRotAtRes);
			int numTotalPrunedPairs = 0;
			if ((useFlags)||(algOption>=3)) //pairs pruning is performed
				numTotalPrunedPairs = countPrunedPairs(rs.getSplitFlags());
			
			if ((numTotalPrunedRot!=numPrunedRot)||(numTotalPrunedPairs!=numPrunedPairs)) { //new rotamers/pairs pruned this run
				
				numPrunedRotThisRun += numTotalPrunedRot - numPrunedRot;
				numPrunedPairsThisRun += numTotalPrunedPairs - numPrunedPairs;
				
				numPrunedRot = numTotalPrunedRot; //update the number of pruned rotamers
				numPrunedPairs = numTotalPrunedPairs;
				numRuns++;
			}
			else { //no more rotamers pruned, so perform the computationally-expensive 2-sp split-DEE and pairs
				
				if (!doDACS){ //DACS will not be performed
					
					if ((algOption>=3)){ //simple Goldstein pairs
						System.out.println("Starting pruning with DEE (full pairs)");	
						
						if (distrDEE){ //distributed full Goldstein pairs DEE
							doDistrDEEMaster(numInAS, numLigRotamers, residueMap, 
									ligType, sParams, ligPresent, resDefault,
									runNameEMatrixMin, runNameEMatrixMax, initEw, prunedRotAtRes, doMinimize, 
									rs, sfFile, rotFile, useFlags, -1, optPairs, minimizeBB, scaleInt, maxIntScale, useEref, eRef);
							
							rs.setSplitFlags((boolean [][])readObject(sfFile)); //get the DE pairs from the distributed run
						}
						else { //perform on a single processor
							rs.DoDEEPairs(numInAS, totalNumRotamers, numLigRotamers, residueMap,
									rotamerIndexOffset, initEw, prunedRotAtRes, null, doMinimize, useFlags, 
									false, false, minimizeBB, scaleInt, maxIntScale);
						}
						System.out.println();
					}
					
					if ((algOption>=2)){ //2-sp conf splitting
						System.out.println("Starting pruning with DEE (2-sp split-DEE)");
						
						if (distrDEE){ //distributed 2-sp split-DEE
							doDistrDEEMaster(numInAS, numLigRotamers, residueMap, 
									ligType, sParams, ligPresent, resDefault,
									runNameEMatrixMin, runNameEMatrixMax, initEw, prunedRotAtRes, doMinimize, 
									rs, sfFile, rotFile, useFlags, 2, optSplit, minimizeBB, scaleInt, maxIntScale, useEref, eRef);
							
							prunedRotAtRes = ((boolean [])readObject(rotFile));
						}
						else { //perform on a single processor			
							prunedRotAtRes = rs.DoDEEConfSplitting(numInAS, totalNumRotamers, numLigRotamers, residueMap,
									rotamerIndexOffset, initEw, prunedRotAtRes, null, doMinimize, useFlags, 2, false, minimizeBB);
						}
						System.out.println();
					}
					
					//check if 2-sp split-DEE and pairs pruned new rotamers
					numTotalPrunedRot = countPrunedRot(prunedRotAtRes);
					numTotalPrunedPairs = 0;
					if ((useFlags)||(algOption>=3)) //pairs pruning is performed
						numTotalPrunedPairs = countPrunedPairs(rs.getSplitFlags());
					
					if ((numTotalPrunedRot==numPrunedRot)&&(numTotalPrunedPairs==numPrunedPairs)) //no more rotamers/pairs pruned
						done = true;
					else { //new rotamers pruned this run
						numPrunedRotThisRun += numTotalPrunedRot - numPrunedRot;
						numPrunedPairsThisRun += numTotalPrunedPairs - numPrunedPairs;
						
						numPrunedRot = numTotalPrunedRot; //update the number of pruned rotamers
						numPrunedPairs = numTotalPrunedPairs;
						numRuns++;
					}
				}
				else //DACS will be performed
					done = true;
			}
			System.out.println("Num pruned rot this run: "+numPrunedRotThisRun);
			System.out.println("Num pruned pairs this run: "+numPrunedPairsThisRun);
			System.out.println();
		}
		
		//long pruneTime = System.currentTimeMillis(); //this line is commented by Swati
		long pruneTimeCPU = CPUTime(); // these lines added to caculate the exact time - Swati
		long pruneTimeUser = UserTime();
		long pruneTimeSystem = SystemTime();
			
		if (!doDACS){ //DACS will not be performed
			
			if (genInteractionGraph) //generate interaction graph
				genInteractionGraph(numInAS, rs, prunedRotAtRes, ligPresent, ligType, numLigRotamers, runName, residueMap, eInteractionCutoff, distCutoff, m, preprocPairs, pairSt);
			
			else { //perform A* search to enumerate conformations
				double bestScore = Math.pow(10,38); //for DEE/A*, the initial best score is the highest possible
				rs.doAStarGMEC(outputConfInfo,true,doMinimize,
						numInAS,totalNumRotamers,rotamerIndexOffset,residueMap,resDefault,numMaxMut,initEw,
						bestScore,null,approxMinGMEC,lambda,minimizeBB,useEref,eRef,doBackrubs,backrubFile);
			}
		}
		else { //DACS			
			int numRotForRes[] = compNumRotForRes(numInAS, rs, numLigRotamers, residueMap);
			BigInteger numInitUnprunedConfs = compNumUnprunedConfs(numInAS, totalNumRotamers, prunedRotAtRes, numLigRotamers,
					numRotForRes, ligPresent);			
			
			int msp[] = new int[numInAS]; //split positions (for DACS)
			for (int i=0; i<msp.length; i++)
				msp[i] = -1;
			
			if (distrDACS) { //distributed DACS (only for level 0)
				
				if (initDepth<=0){
					System.out.println("ERROR: distributed DACS called with 'initDepth="+initDepth+"' partitioning positions; use a positive value");
					System.exit(1);
				}
				if (subDepth<0)
					subDepth = 0;
				
				distrDEE = false; //do not perform both distributed DACS and distributed DEE
				
				//choose the major splitting positions
				for (int i=0; i<initDepth; i++)
					msp[i] = chooseSplitPos(numInAS,prunedRotAtRes,totalNumRotamers,numRotForRes, msp, i, minRatioDiff);
				
				int maxNumPartitions = 1;
				for (int i=0; i<initDepth; i++)
					maxNumPartitions *= numRotForRes[msp[i]];
				
				OneMutation resumeResults[] = null;			
				if (resumeSearch){ //read resume results
					System.out.println("Reading resume results..");
					resumeResults = new OneMutation[maxNumPartitions];
					for(int q=0;q<resumeResults.length;q++)
						resumeResults[q] = new OneMutation();
					resumeResults = readResumeFile(resumeResults,resumeFilename,numInAS,true,false,initDepth);
					System.out.println("Read "+resumeResults.length+" completed partitions.");
					System.out.println();
				}
				
				doDistrDACSMaster(runName, numInAS, rs, numLigRotamers, residueMap, ligPresent, resDefault,
						rotFile, prunedRotAtRes, algOption, sfFile, useFlags, initEw, pruningE, initDepth, msp,
						numInitUnprunedConfs, diffFact, outputPruneInfo, outputConfInfo, minRatioDiff, doMinimize,
						runNameEMatrixMin, runNameEMatrixMax, ligType, sParams, approxMinGMEC, 
						lambda, numRotForRes, resumeResults, minimizeBB, numMaxMut, scaleInt, maxIntScale, 
						useEref, eRef, doBackrubs, backrubFile, subDepth);
			}
			else { //single-processor DACS
				
				initDepth = 0; //only used for distributed DACS
				if (subDepth<=0){
					System.out.println("ERROR: single-processor DACS called with 'subDepth="+subDepth+"' partitioning positions; use a positive value");
					System.exit(1);
				}
				
				PrintStream logPS = setupOutputFile(outputPruneInfo);
			
				doDACS(numInAS, rs, numLigRotamers, residueMap, ligPresent, resDefault,
						rotFile, prunedRotAtRes, algOption, sfFile, useFlags, initEw, pruningE, initDepth, 0, logPS, msp,
						numInitUnprunedConfs, diffFact, outputConfInfo, minRatioDiff, doMinimize, runNameEMatrixMin,
						runNameEMatrixMax, distrDEE, ligType, sParams, approxMinGMEC, lambda, null, 
						null, minimizeBB, numMaxMut, scaleInt, maxIntScale, useEref, eRef, doBackrubs, backrubFile, subDepth);
			}
		}
		
		/*long stopTime = System.currentTimeMillis(); This section commented by Swati
		
		System.out.println("Pruning time: "+((pruneTime-startTime)/(60.0*1000.0)));
		if (genInteractionGraph)
			System.out.println("Graph generation time: "+((stopTime-pruneTime)/(60.0*1000.0)));
		else
			System.out.println("Enumeration/DACS time: "+((stopTime-pruneTime)/(60.0*1000.0)));
		System.out.println("DEE execution time: "+((stopTime-startTime)/(60.0*1000.0)));
		System.out.println("DEE done");
		System.out.println("Total execution time: "+((stopTime-startTimeAll)/(60.0*1000.0)));*/

		long stopTimeCPU = CPUTime(); // these all lines are added to calculate the exact time and print it - Swati 
		long stopTimeUser = UserTime();
		long stopTimeSystem = SystemTime();
		
		System.out.println("Pruning time: CPU "+((pruneTimeCPU-startTimeCPU)/(60.0*1000000000.0)));
		System.out.println("Pruning time: User "+((pruneTimeUser-startTimeUser)/(60.0*1000000000.0)));
		System.out.println("Pruning time: System "+((pruneTimeSystem-startTimeSystem)/(60.0*1000000000.0)));
		if (genInteractionGraph)
		{
			System.out.println("Graph generation time: CPU "+((stopTimeCPU-pruneTimeCPU)/(60.0*1000000000.0)));
			System.out.println("Graph generation time: User "+((stopTimeUser-pruneTimeUser)/(60.0*1000000000.0)));
			System.out.println("Graph generation time: System "+((stopTimeSystem-pruneTimeSystem)/(60.0*1000000000.0)));
		}
		else
		{
			System.out.println("Enumeration/DACS time: CPU "+((stopTimeCPU-pruneTimeCPU)/(60.0*1000000000.0)));
			System.out.println("Enumeration/DACS time: User "+((stopTimeUser-pruneTimeUser)/(60.0*1000000000.0)));
			System.out.println("Enumeration/DACS time: System "+((stopTimeSystem-pruneTimeSystem)/(60.0*1000000000.0)));
		}
		System.out.println("DEE execution time: CPU "+((stopTimeCPU-startTimeCPU)/(60.0*1000000000.0)));
		System.out.println("DEE execution time: User "+((stopTimeUser-startTimeUser)/(60.0*1000000000.0)));
		System.out.println("DEE execution time: System "+((stopTimeSystem-startTimeSystem)/(60.0*1000000000.0)));
		System.out.println("DEE done");
		System.out.println("Total execution time: CPU "+((stopTimeCPU-startTimeAllCPU)/(60.0*1000000000.0)));
		System.out.println("Total execution time: User "+((stopTimeUser-startTimeAllUser)/(60.0*1000000000.0)));
		System.out.println("Total execution time: System "+((stopTimeSystem-startTimeAllSystem)/(60.0*1000000000.0)));
		//end of DEE section
		/////////////////////////////////////////////////////////////
	}
	
	/**
	 * Performs the DACS partition-specific computation. The parameters 'rs' and 'rs.sysLR' must be valid.
	 */
	private void doDACS(int numInAS, RotamerSearch rs, int numLigRotamers, int residueMap[],
			boolean ligPresent, String resDefault[], String rotFile,
			boolean prunedRotAtRes[], int algOption, String sfFile, boolean useFlags, float initEw, float pruningE,
			int initDepth, int curDepth, PrintStream logPS, int majorSplitPos[], BigInteger numInitUnprunedConfs,
			int diffFact, String outputConfInfo, double minRatioDiff, boolean doMinimize, String minPEM, 
			String maxPEM, boolean distrDEE, String ligType, ParamSet sParams, 
			boolean approxMinGMEC, float lambda, int partIndex[], CommucObj cObj, 
			boolean minimizeBB, int numMaxMut, boolean scaleInt, float maxIntScale, boolean useEref, float eRef[], 
			boolean doBackrubs, String backrubFile, int subDepth){
		
		if (curDepth>=(initDepth+subDepth))
			return;
		
		System.out.println("Starting pruning with DEE (DACS).");
		
		//prunedRotAtRes[] should not be modified here, in order to be able to distinguish
		//	newly pruned rotamers and rotamers pruned by Bounds or DEE
		
		//the num rotamers for each AS residue and the ligand (if present);
		//	sysLR in rs must be valid (with all the possible AA's for each residue position)
		int numRotForRes[] = compNumRotForRes(numInAS, rs, numLigRotamers, residueMap);
		
		if (curDepth>=initDepth) {//sub-partition, so majorSplitPos[curDepth] is unknown; compute it
			majorSplitPos[curDepth] = chooseSplitPos(numInAS,prunedRotAtRes,totalNumRotamers,numRotForRes,
					majorSplitPos, curDepth, minRatioDiff);//the splitting position
		}
		
		int numPartitions = numRotForRes[majorSplitPos[curDepth]]; //the number of partitions
		int numPrunedPartitions = 0; //the number of partitions with no unpruned confs
		
		System.out.println("Current depth: "+curDepth);
		System.out.println("Major splitting residue number: "+majorSplitPos[curDepth]);
		System.out.println();
		
		//map rotamer index to rot num for the splitting residue
		int indexMap[] = getIndexMap(numPartitions,rs,majorSplitPos[curDepth],residueMap,totalNumRotamers,rotamerIndexOffset);
		
		//count the total num confs and the num unpruned confs
		BigInteger numConfsTotalForPartition[] = new BigInteger[numPartitions];
		BigInteger numConfsUnprunedForPartition[] = new BigInteger[numPartitions];
		//BigInteger numInitConfsUnprunedForPartition[] = new BigInteger[numPartitions];
		BigInteger numTotalConfs = new BigInteger("0");
		BigInteger numUnprunedConfs = new BigInteger("0");
		BigInteger numEvaluatedConfs = new BigInteger("0");
		
		//update the best energy found
		pruningE = (float)Math.min(pruningE, rs.getBestE());
		double bestScore = pruningE;
		
		boolean savedSpFlags[][] = null;
		if ((useFlags)||(algOption>=3))
			savedSpFlags = rs.getSplitFlags(savedSpFlags);
		
		//determine the prunings for each of the sub-solutions (the partitions)
		boolean prunedForPartition[][] = new boolean[numPartitions][prunedRotAtRes.length];
		for (int i=0; i<prunedForPartition.length; i++){
			
			if ((curDepth>=initDepth)||(indexMap[i]==partIndex[curDepth])) { //sub-partitions or current partition is the partition distributed for computation
			
				//copy the prunings from before conf splitting (from Bounds and simple traditional DEE)
				System.arraycopy(prunedRotAtRes,0,prunedForPartition[i],0,prunedRotAtRes.length);
				
				int curPartIndex = indexMap[i]; //the rotamer index of the partitioning rotamer
				
				//artificially set all rotamers at the splitting position, other than the rotamer 
				//	for the current partition, to pruned, so that there will be only one rotamer at
				//	that residue position when the conf splitting criterion is applied;
				//	when done, subtract the artifcially pruned rotamers from the total number of pruned
				boolean indToUnprune[] = new boolean[numPartitions];
				for (int j=0; j<indToUnprune.length; j++)
					indToUnprune[j] = false;
				
				//check the partition only if the current partitioning rotamer is not already pruned
				if (!prunedRotAtRes[curPartIndex]){
					
					for (int j=0; j<numPartitions; j++){
						if (j!=i) {//not the rotamer for the current partition
							int curInd = indexMap[j]; //the index of the current rotamer
							if (!prunedForPartition[i][curInd]){ //not pruned by the other DEE methods
								prunedForPartition[i][curInd] = true;
								indToUnprune[j] = true;
							}
						}
					}
					
					if ((useFlags)||(algOption>=3))
						rs.setSplitFlags(savedSpFlags); //for each partition, reset the flags to the globally valid ones
					
					int numPrunedRot = countPrunedRot(prunedForPartition[i]);
					int numPrunedPairs = 0;
					if ((useFlags)||(algOption>=3))
						numPrunedPairs = countPrunedPairs(rs.getSplitFlags());
					int numPrunedRotThisRun = 0;
					int numPrunedPairsThisRun = 0;
					boolean done = false;
					int numRuns = 1;
					
					while (!done){ //repeat the pruning cycle until no more rotamers are pruned	
						
						numPrunedRotThisRun = 0; 
						numPrunedPairsThisRun = 0;
						
						System.out.println("Starting DEE cycle run: "+numRuns);
						
						if (doMinimize) //precompute the interval terms in the MinDEE criterion
							rs.doCompMinDEEIntervals(numInAS, totalNumRotamers, numLigRotamers, residueMap, rotamerIndexOffset, 
									prunedForPartition[i], scaleInt, maxIntScale);
						
						System.out.println("Starting pruning with DEE (simple Goldstein)");		
						prunedForPartition[i] = rs.DoDEEGoldstein(numInAS, totalNumRotamers, numLigRotamers, residueMap,
								rotamerIndexOffset, initEw, prunedForPartition[i], doMinimize, useFlags, minimizeBB);			
						System.out.println();
						
						if ((algOption>=3)){ //simple Goldstein pairs
							System.out.println("Starting pruning with DEE (mb pairs)");			
							rs.DoDEEPairs(numInAS, totalNumRotamers, numLigRotamers, residueMap,
									rotamerIndexOffset, initEw, prunedForPartition[i], null, doMinimize, 
									useFlags, true, false, minimizeBB, scaleInt, maxIntScale);
							System.out.println();
						}
						
						if ((useFlags)||(algOption>=3)){
							System.out.println("Starting pruning with Bounding Flags");			
							rs.DoBoundFlags(numInAS, totalNumRotamers, numLigRotamers, residueMap,
								rotamerIndexOffset, pruningE, prunedForPartition[i], initEw, useFlags);
							System.out.println();
						}
							
						System.out.println("Starting pruning with DEE (1-sp split-DEE)");
						prunedForPartition[i] = rs.DoDEEConfSplitting(numInAS, totalNumRotamers, numLigRotamers, residueMap,
								rotamerIndexOffset, initEw, prunedForPartition[i], null, doMinimize, 
								useFlags, 1, false, minimizeBB);
						System.out.println();
						
						System.out.println("Starting pruning with Bounds");			
						prunedForPartition[i]= rs.DoMinBounds(numInAS, totalNumRotamers, numLigRotamers, residueMap,
							rotamerIndexOffset, pruningE, prunedForPartition[i], initEw, useFlags, false);
						System.out.println();
						
						//check how many rotamers/pairs are pruned this run
						int numTotalPrunedRot = countPrunedRot(prunedForPartition[i]);
						int numTotalPrunedPairs = 0;
						if ((useFlags)||(algOption>=3))
							numTotalPrunedPairs = countPrunedPairs(rs.getSplitFlags());
						
						if ((numTotalPrunedRot!=numPrunedRot)||(numTotalPrunedPairs!=numPrunedPairs)) { //new rotamers/pairs pruned this run
							
							numPrunedRotThisRun += numTotalPrunedRot - numPrunedRot;
							numPrunedPairsThisRun += numTotalPrunedPairs - numPrunedPairs;
							
							numPrunedRot = numTotalPrunedRot; //update the number of pruned rotamers
							numPrunedPairs = numTotalPrunedPairs;
							numRuns++;
						}
						else { //no more rotamers pruned, so perform the computationally-expensive 2-sp split-DEE and pairs
							
							if ((algOption>=3)&&(curDepth>=(initDepth+subDepth-1))){ //simple Goldstein pairs					
								System.out.println("Starting pruning with DEE (full pairs)");	
								
								if (distrDEE){ //distributed full Goldstein pairs DEE
									doDistrDEEMaster(numInAS, numLigRotamers, residueMap, 
											ligType, sParams, ligPresent, resDefault, 
											minPEM, maxPEM, initEw, prunedForPartition[i], doMinimize, rs, sfFile, rotFile, 
											useFlags, -1, optPairs, minimizeBB, scaleInt, maxIntScale, useEref, eRef);
									
									rs.setSplitFlags((boolean [][])readObject(sfFile)); //get the DE pairs from the distributed run
								}
								else { //perform on a single processor
									rs.DoDEEPairs(numInAS, totalNumRotamers, numLigRotamers, residueMap,
											rotamerIndexOffset, initEw, prunedForPartition[i], null, doMinimize, 
											useFlags, false, false, minimizeBB, scaleInt, maxIntScale);
								}
								System.out.println();
							}
							
							if ((algOption>=2)){ //2-sp conf splitting					
								System.out.println("Starting pruning with DEE (2-sp split-DEE)");
								
								if (distrDEE) { //distributed 2-sp split-DEE
									doDistrDEEMaster(numInAS, numLigRotamers, residueMap, 
											ligType, sParams, ligPresent, resDefault,
											minPEM, maxPEM, initEw, prunedForPartition[i], doMinimize, rs, sfFile, rotFile, 
											useFlags, 2, optSplit, minimizeBB, scaleInt, maxIntScale, useEref, eRef);
									
									prunedForPartition[i] = ((boolean [])readObject(rotFile));
								}
								else {
									prunedForPartition[i] = rs.DoDEEConfSplitting(numInAS, totalNumRotamers, numLigRotamers, residueMap,
											rotamerIndexOffset, initEw, prunedForPartition[i], null, doMinimize, 
											useFlags, 2, false, minimizeBB);
								}
								System.out.println();
							}
							
							//check if 2-sp split-DEE and pairs pruned new rotamers/pairs
							numTotalPrunedRot = countPrunedRot(prunedForPartition[i]);
							if ((useFlags)||(algOption>=3))
								numTotalPrunedPairs = countPrunedPairs(rs.getSplitFlags());
							
							if ((numTotalPrunedRot==numPrunedRot)&&(numTotalPrunedPairs==numPrunedPairs)) //no more rotamers/pairs pruned
								done = true;
							else { //new rotamers pruned this run
								numPrunedRotThisRun += numTotalPrunedRot - numPrunedRot;
								numPrunedPairsThisRun += numTotalPrunedPairs - numPrunedPairs;
								
								numPrunedRot = numTotalPrunedRot; //update the number of pruned rotamers
								numPrunedPairs = numTotalPrunedPairs;
								numRuns++;
							}
						}
						System.out.println("Num pruned rot this run: "+numPrunedRotThisRun);
						System.out.println("Num pruned pairs this run: "+numPrunedPairsThisRun);
						System.out.println();
					}
				}
					
				//count the number of pruned rotamers for each residue (except for the splitting residue,
				//	which always has only 1 available rotamer for the current partition)
				int numPrunedRotForRes[] = new int[numInAS]; //after pruning for this partition
				//int numInitPrunedRotForRes[] = new int[numInAS]; //initial prunings, before pruning for this partition
				for (int j=0; j<numInAS; j++){
					numPrunedRotForRes[j] = 0;
					if (j!=majorSplitPos[curDepth]){
						for (int k=0; k<totalNumRotamers; k++){ //pruned rot are true and must be in the current set of allowed AA
							if (prunedForPartition[i][j*totalNumRotamers + k])
								numPrunedRotForRes[j]++;
							//if (prunedRotAtRes[j*totalNumRotamers + k])
								//numInitPrunedRotForRes[j]++;
						}
					}
					else {// j==majorSplitPos
						numPrunedRotForRes[j] = 0;
						//numInitPrunedRotForRes[j] = 0;
					}
				}
				int numPrunedLigRot = 0;
				//int numInitPrunedLigRot = 0;
				for (int k=0; k<numLigRotamers; k++){
					if (prunedForPartition[i][numInAS*totalNumRotamers + k])
						numPrunedLigRot++;
					//if (prunedRotAtRes[numInAS*totalNumRotamers + k])
					//	numInitPrunedLigRot++;
				}
				
				//count the total num confs and the num unpruned confs
				numConfsTotalForPartition[i] = new BigInteger("1");
				numConfsUnprunedForPartition[i] = new BigInteger("1");
				//numInitConfsUnprunedForPartition[i] = new BigInteger("1");
				if (prunedRotAtRes[curPartIndex]){ //current partitioning rotamer already pruned, so no unpruned confs for this partition
					numConfsUnprunedForPartition[i] = new BigInteger("0");
					numPrunedPartitions++;
				}
				
				for (int j=0; j<numInAS; j++){
					if (!(isSplitRes(j,majorSplitPos,curDepth))){ //the split residues contribute only 1 rotamer
						numConfsTotalForPartition[i] = numConfsTotalForPartition[i].multiply(BigInteger.valueOf(numRotForRes[j]));
						numConfsUnprunedForPartition[i] = numConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numRotForRes[j]-numPrunedRotForRes[j]));
						//numInitConfsUnprunedForPartition[i] = numInitConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numRotForRes[j]-numInitPrunedRotForRes[j]));
					}
				}
				if(ligPresent){
					numConfsTotalForPartition[i] = numConfsTotalForPartition[i].multiply(BigInteger.valueOf(numLigRotamers));
					numConfsUnprunedForPartition[i] = numConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numLigRotamers-numPrunedLigRot));
					//numInitConfsUnprunedForPartition[i] = numInitConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numLigRotamers-numInitPrunedLigRot));
				}			
	
				numTotalConfs = numTotalConfs.add(numConfsTotalForPartition[i]);
				numUnprunedConfs = numUnprunedConfs.add(numConfsUnprunedForPartition[i]);
				
				BigInteger pruneDiffFact = BigInteger.valueOf(10).pow(diffFact);
				
				System.out.println("Num unpruned confs: "+numConfsUnprunedForPartition[i]+" diffFact: "+pruneDiffFact);
				System.out.println();
				System.out.println();
				
				//output pruning info to file
				logPS.print("curDepth: "+curDepth+" curPartition: "+i+" majorSplitPos: "+majorSplitPos[curDepth]+" ");
				logPS.print(numConfsTotalForPartition[i]+" "+numInitUnprunedConfs+" "+numConfsUnprunedForPartition[i]);
				logPS.println();
				logPS.println();logPS.flush();
				
				//if ((curDepth+1<maxDepth)&&(numConfsUnprunedForPartition[i].compareTo(numInitUnprunedConfs.divide(pruneDiffFact))==1)){ //not enough pruned, so partition at new depth
				if ( (curDepth+1<initDepth) || 
						( (curDepth+1<(initDepth+subDepth)) && (numConfsUnprunedForPartition[i].compareTo(pruneDiffFact)==1) ) ){ //not enough pruned, so partition at new depth
					doDACS(numInAS, rs, numLigRotamers, residueMap, ligPresent, resDefault,
							rotFile, prunedForPartition[i], algOption, sfFile, useFlags, initEw, pruningE, initDepth, curDepth+1, 
							logPS, majorSplitPos, numConfsUnprunedForPartition[i], diffFact, outputConfInfo, minRatioDiff,
							doMinimize, minPEM, maxPEM, distrDEE, ligType, sParams, approxMinGMEC, 
							lambda, partIndex, cObj, minimizeBB, numMaxMut, scaleInt, maxIntScale, useEref, eRef, doBackrubs, backrubFile, subDepth);
				}
				else if (!prunedRotAtRes[curPartIndex]){ //if enough pruned or maxDepth partitioning reached, do the rotamer search
					
					bestScore = Math.min(bestScore,rs.getBestE());//best E for the partitions so far
					
					//Do the rotamer search
					rs.doAStarGMEC(outputConfInfo,true,doMinimize,numInAS,totalNumRotamers,rotamerIndexOffset,residueMap,
							resDefault,numMaxMut,initEw,bestScore,null,approxMinGMEC,lambda,minimizeBB,useEref,eRef,doBackrubs,backrubFile);
					
					numEvaluatedConfs = numEvaluatedConfs.add(rs.numConfsEvaluated); //add the evaluated confs for this partition
					pruningE = (float)Math.min(pruningE,rs.getBestE());//update cutoff energy for MinBounds
				}
				
				//unprune the artificially pruned indices
				for (int j=0; j<indToUnprune.length; j++){
					if (indToUnprune[j]){
						prunedForPartition[i][indexMap[j]] = false;
					}
				}
			}
		}
		
		if (cObj==null){ //not distributed DACS
		
			System.out.println("numTotalConfs: "+numTotalConfs+"; numUnprunedConfs: "+numUnprunedConfs+"; numEvaluatedConfs: "+numEvaluatedConfs);
			for (int i=0; i<numPartitions; i++)System.out.print(numConfsTotalForPartition[i]+" ");System.out.println();
			for (int i=0; i<numPartitions; i++)System.out.print(numConfsUnprunedForPartition[i]+" ");
			System.out.println();
			System.out.println("Major splitting residue number: "+majorSplitPos[curDepth]);
			System.out.println("Number of partitions: "+numPartitions);
			System.out.println("Number of non-zero partitions: "+(numPartitions-numPrunedPartitions));
			System.out.println("Additional pruned rotamers: ");
			
			//count the number of partitions for which a rotamer is pruned (counting only the rotamers
			//	not pruned by Bounds or simple traditional DEE)
			int countNumPartitionsPrunedRot[] = new int[prunedRotAtRes.length];
			for (int i=0; i<prunedRotAtRes.length; i++){//for each rotamer
				countNumPartitionsPrunedRot[i] = 0;
				if (!prunedRotAtRes[i]){ //only if not pruned by the other two methods
					for (int j=0; j<numPartitions; j++){ //check for each partition
						if (prunedForPartition[j][i])
							countNumPartitionsPrunedRot[i]++;
					}
				}
				
				//output information
				if (countNumPartitionsPrunedRot[i]>0)
					System.out.println("index: "+i+"; num partitions in which pruned: "+countNumPartitionsPrunedRot[i]);
			}
		}
		else { //distributed DACS
			cObj.bestScore = cObj.bestScore.min(BigDecimal.valueOf(rs.getBestE()));
		}
	}
	
	/**
	 * Implements threads for DACS: allows the current best score among all partitions to be distributed to every partition;
	 * This thread performs the DACS computation, while the main thread monitors for updates of the best energy from the other partitions;
	 * The communication is performed via the common RotamerSearch object (synchronized access to the bestEMin variable)
	*/
	private class DACSthread implements Runnable {
		private RotamerSearch rs = null;
		private CommucObj cObj = null;
		private PrintStream logPS = null;
		String outputConfInfo = null;
		DACSthread(RotamerSearch rsP, CommucObj cObjP, PrintStream lP, String ociP){
			rs = rsP;
			cObj = cObjP;
			logPS = lP;
			outputConfInfo = ociP;
		}
		public void run(){
			//Perform DACS
			doDACS(cObj.numInAS, rs, cObj.numLigRotamers, cObj.residueMap,
					cObj.ligPresent, cObj.resDefault, null,
					cObj.prunedRot, cObj.algOption, cObj.sfFileIn, cObj.useSF, cObj.initEw, cObj.pruningE,
					cObj.initDepth, 0, logPS, cObj.msp, cObj.numInitUnprunedConf,
					cObj.diffFact, outputConfInfo, cObj.minRatioDiff, cObj.doMinimization, null, 
					null, false, cObj.ligType, cObj.params,  
					cObj.approxMinGMEC, cObj.lambda, cObj.partIndex, cObj, cObj.minimizeBB, cObj.numMutations,
					cObj.scaleInt, cObj.maxIntScale, cObj.useEref, cObj.eRef, cObj.doBackrubs, cObj.backrubFile, cObj.subDepth);
		}
	}
	
	//Compute the number of unpruned conformations 
	private BigInteger compNumUnprunedConfs(int numInAS, int totalNumRotamers, boolean prunedRotAtRes[], int numLigRotamers,
			int numRotForRes[], boolean ligPresent) {
		
		int numPrunedRotForRes[] = new int[numInAS]; //after pruning for this partition
		for (int j=0; j<numInAS; j++){
			numPrunedRotForRes[j] = 0;
			for (int k=0; k<totalNumRotamers; k++){ //pruned rot are true and must be in the current set of allowed AA
				if (prunedRotAtRes[j*totalNumRotamers + k])
					numPrunedRotForRes[j]++;
			}
		}
		int numPrunedLigRot = 0;
		for (int k=0; k<numLigRotamers; k++){
			if (prunedRotAtRes[numInAS*totalNumRotamers + k])
				numPrunedLigRot++;
		}
		
		//count the total num confs and the num unpruned confs
		BigInteger numConfsTotal = new BigInteger("1");
		BigInteger numConfsUnpruned = new BigInteger("1");
		
		for (int j=0; j<numInAS; j++){
			numConfsTotal = numConfsTotal.multiply(BigInteger.valueOf(numRotForRes[j]));
			numConfsUnpruned = numConfsUnpruned.multiply(BigInteger.valueOf(numRotForRes[j]-numPrunedRotForRes[j]));
		}
		if(ligPresent){
			numConfsTotal = numConfsTotal.multiply(BigInteger.valueOf(numLigRotamers));
			numConfsUnpruned = numConfsUnpruned.multiply(BigInteger.valueOf(numLigRotamers-numPrunedLigRot));
		}	
		
		return numConfsUnpruned;
	}
	
	//Choose the splitting position (use only AS positions; the ligand is chosen not to be a splitting position)
	private int chooseSplitPosRandom(int numInAS){		
		Random randNum = new Random();
		return randNum.nextInt(numInAS);
	}
	
	//Choose the splitting position (use only AS positions; the ligand is chosen not to be a splitting position);
	//	Choose the AS residue position with the smallest fraction of pruned rotamers (from MinBounds and simple MinDEE)
	private int chooseSplitPos(int numInAS, boolean prunedRotAtRes[], int numTotalRotamers, int numRotForRes[],
			int majorSplitPos[], int curDepth, double minRatioDiff){		
		
		final int minPartitions = 5; //the min number of rotamers that a splitting residue can have
		
		double pruneRatio[] = new double[numInAS];
		int minPos = -1;
		double minRatio = (double)Math.pow(10,38);
		for (int curRes=0; curRes<numInAS; curRes++){
			
			if (!(isSplitRes(curRes,majorSplitPos,curDepth))){
				if (numRotForRes[curRes]>=minPartitions){ //do not split at residues with very small number of rotamers
					int curPruned = 0;
					for (int curRot=0; curRot<numTotalRotamers; curRot++){
						if (prunedRotAtRes[curRes*numTotalRotamers+curRot]){ //cur rot is pruned (pruned rotamers are necessarily in the cur set of allowed AAs)
							curPruned++;
						}
					}
					pruneRatio[curRes] = (double)curPruned/numRotForRes[curRes];
					if (minRatio>=pruneRatio[curRes]){
						if ((minPos==-1)||(curRes<minPos)||(minRatio>=pruneRatio[curRes]+minRatioDiff)) {//preference to split at lower-numbered residues
							minRatio = pruneRatio[curRes];
							minPos = curRes;
						}
					}
				}
			}
		}
		
		if (minPos!=-1){
			//System.out.println("minPos: "+minPos);
			//for (int i=0;i<numInAS;i++)System.out.print(pruneRatio[i]+" ");System.out.println();
			return minPos;
		}
		else //if split position not hosen, choose randomly
			return chooseSplitPosRandom(numInAS);
	}
	
	//Check if the residue curRes is one of the splitRes
	private boolean isSplitRes(int curRes, int majorSplitPos[], int curDepth){
		for (int i=0; i<=curDepth; i++){
			if (curRes==majorSplitPos[i])
				return true;
		}
		return false;
	}
	
	//Compute the number of rotamers for each residue position (assign to numRotForRes[])
	private int [] compNumRotForRes(int numInAS, RotamerSearch rs, int numLigRot, int residueMap[]){
		
		int numRotForRes[] = new int[numInAS];
		boolean ligPresent = (numLigRot!=0); //ligand present
		int treeLevels = numInAS;
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
				for (int i=0; i<rs.sysLR.getNumAllowable(residueMap[curLevel]); i++){ //add the rot for all allowable AA at this residue
					int newRot = rl.getNumRotForAAtype(rs.sysLR.getIndexOfNthAllowable(residueMap[curLevel],i));
					if (newRot==0) //GLY or ALA
						newRot = 1;
					curNumRot += newRot; 
				}
			}
			numRotForRes[curLevel] = curNumRot;
		}
		return numRotForRes;
	}
	
	//Get the mapping between rotamer indices (into the pruning matrix) and the number of the
	//	current rotamer for the giveen residue; assumes sysLR in rs is valid (all allowables for the AS residues)
	private int [] getIndexMap(int numPartitions, RotamerSearch rs, int curRes, int residueMap[], int numTotalRot,
			int rotamerIndexOffset[]){
		
		int indexMap[] = new int[numPartitions];
		int indNum = 0;
		for (int AA=0; AA<rs.sysLR.getNumAllowable(residueMap[curRes]); AA++){ //for each AA for the given AS residue
			int curAA = rs.sysLR.getIndexOfNthAllowable(residueMap[curRes],AA);
			int numRotForAA = rl.getNumRotForAAtype(curAA);
			if (numRotForAA==0) //GLY or ALA
				numRotForAA = 1;
			
			for (int curRot=0; curRot<numRotForAA; curRot++){ //for each rot for the given AA
				indexMap[indNum] = curRes*numTotalRot + rotamerIndexOffset[curAA] + curRot;
				indNum++;
			}
		}
		return indexMap;
	}
	
	//Setup the file with name filename for output
	private PrintStream setupOutputFile(String fileName){
		PrintStream logPS = null; //the output file for conf info
		try {			
			FileOutputStream fileOutputStream = new FileOutputStream(fileName);
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			logPS = new PrintStream( bufferedOutputStream );
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
		}
		return logPS;
	}	
///////////////////////////////////////////////////////////////////////////
//	End of DEE section
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
//	Distributed DACS section
///////////////////////////////////////////////////////////////////////////
	/**
	 * Handles the distribution of the DACS computation to the set of available processors
	*/
	private void doDistrDACSMaster(String runName, int numInAS, RotamerSearch rs, int numLigRotamers, int residueMap[],
			boolean ligPresent, String resDefault[], String rotFile,
			boolean prunedRotAtRes[], int algOption, String sfFile, boolean useFlags, float initEw, float pruningE,
			int initDepth, int majorSplitPos[], BigInteger numInitUnprunedConfs,
			int diffFact, String outputPruneInfo, String outputConfInfo, double minRatioDiff, boolean doMinimize, String minPEM, 
			String maxPEM, String ligType, ParamSet sParams, 
			boolean approxMinGMEC, float lambda, int numRotForRes[], OneMutation resumeResults[], boolean minimizeBB, int numMaxMut,
			boolean scaleInt, float maxIntScale, boolean useEref, float eRef[], boolean doBackrubs, String backrubFile, int subDepth){
		
		System.out.println("Starting DACS (distributed)");		
		System.out.println("Forming DACS partitions..");
		
		int indexMap[][] = new int[initDepth][];
		int numPartitions[] = new int[indexMap.length];
		int maxNumPartitions = 1;
		
		//map rotamer index to rot num for the splitting residues
		for (int i=0; i<indexMap.length; i++){
			numPartitions[i] = numRotForRes[majorSplitPos[i]]; //the number of partitions
			indexMap[i] = getIndexMap(numPartitions[i],rs,majorSplitPos[i],residueMap,totalNumRotamers,rotamerIndexOffset);
			maxNumPartitions *= numPartitions[i];
		}
	
		//get the partitions
		OneMutation mutArray[] = formDACSpartitions(maxNumPartitions, initDepth, indexMap, numPartitions, prunedRotAtRes, majorSplitPos);		
		
		if (resumeResults!=null){ //remove completed partitions
			int curMut = 0;
			OneMutation tmpArray2[] = new OneMutation[mutArray.length];
			for (int i=0; i<mutArray.length; i++){
				boolean partFound = true;
				for (int j=0; j<resumeResults.length; j++){
					partFound = true;
					for (int k=0; k<initDepth; k++){
						if (mutArray[i].resMut[k]!=resumeResults[j].resMut[k]){
							partFound = false;
							break;
						}
					}
					if (partFound) //partition already computed
						break;
				}
				if (!partFound){
					tmpArray2[curMut] = new OneMutation();
					tmpArray2[curMut].mutNum = mutArray[i].mutNum;
					tmpArray2[curMut].resMut = mutArray[i].resMut;
					curMut++;
				}
			}
			OneMutation tmpArray[] = new OneMutation[curMut]; //trim the size of the partition array
			System.arraycopy(tmpArray2, 0, tmpArray, 0, curMut);
			mutArray = tmpArray;
			
			System.out.println("Number non-zero partitions after removing completed results: "+mutArray.length);
			
			//update the best energy so far
			for (int j=0; j<resumeResults.length; j++){
				pruningE = (float)Math.min(pruningE, resumeResults[j].score.doubleValue());
			}
		}
		
		//output the rs object
		outputObject(rs,rotFile);
		
		//sort the partitions
		sortDACSpartitions(mutArray,initDepth,majorSplitPos,numPartitions,prunedRotAtRes,indexMap,numInAS,
				numLigRotamers,residueMap,pruningE,initEw,rs);
		
		System.out.println();
		
		MutationManager mutMan = new MutationManager(runName,mutArray,false);
		mutMan.setResidueMap(residueMap);
		mutMan.setResDefault(resDefault);
		mutMan.setRotamerIndexOffset(rotamerIndexOffset);
		mutMan.setLigType(ligType);
		mutMan.setarpFilenameMin(minPEM);
		mutMan.setarpFilenameMax(maxPEM);
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
		mutMan.setNumInAS(numInAS);
		mutMan.setnumLigRotamers(numLigRotamers);
		mutMan.numTotalRotamers(totalNumRotamers);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setDoBackrubs(doBackrubs);
		mutMan.setBackrubFile(backrubFile);
		mutMan.setCalculateVolumes(false);
		mutMan.setNumMutations(numMaxMut);
		mutMan.setInitEw(initEw);
		mutMan.setPruningE(pruningE);
		mutMan.setLigPresent(ligPresent);
		mutMan.setUseSF(useFlags);
		mutMan.setPrunedRot(prunedRotAtRes);
		mutMan.setSfFile(sfFile);
		mutMan.setRotFile(rotFile);
		mutMan.setDistrDACS(true);
		mutMan.setDistrDEE(false);
		mutMan.setBestScore(new BigDecimal(pruningE)); //the best E initially is the pruningE read from the parameter file
		mutMan.setAlgOption(algOption);
		mutMan.setInitDepth(initDepth);
		mutMan.setSubDepth(subDepth);
		mutMan.setDiffFact(diffFact);
		mutMan.setMinRatioDiff(minRatioDiff);
		mutMan.setNumInitUnprunedConf(numInitUnprunedConfs);
		mutMan.setOutputPruneInfo(outputPruneInfo);
		mutMan.setOutputConfInfo(outputConfInfo);
		mutMan.setMSP(majorSplitPos);
		mutMan.setApproxMinGMEC(approxMinGMEC);
		mutMan.setLambda(lambda);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setScaleInt(scaleInt);
		mutMan.setMaxIntScale(maxIntScale);
		mutMan.setUseEref(useEref);
		mutMan.setEref(eRef);

		try{
			handleDoMPIMaster(mutMan,mutArray.length);
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}
		
		//Delete the temporary rs file
		if (!(new File(rotFile)).delete()){
			System.out.println("ERROR: cannot delete file "+rotFile);
			System.exit(1);
		}
	}
	
	// Distributed DACS Slave function
	private CommucObj doDistrDACSSlave(CommucObj cObj) {
				
		long startTime = System.currentTimeMillis();
						
		boolean ligPresent = cObj.ligPresent;

		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,cObj.params,ligPresent,cObj.ligType);
		
		RotamerSearch rs = (RotamerSearch)readObject(cObj.rotFileIn); //load the saved rs from master
		
		String fn = "";
		for (int i=0; i<cObj.partIndex.length; i++)
			fn += ("_"+cObj.partIndex[i]);
		String outputConfInfo = ("./conf_info/"+cObj.outputConfInfo+fn);
		PrintStream logPS = setupOutputFile("./conf_info/"+cObj.outputPruneInfo+fn);
		
		float otherBestE = cObj.pruningE; //the best energy from other partitions
		
		Thread t = new Thread(new DACSthread(rs,cObj,logPS,outputConfInfo)); //create a new thread for performing the DACS search
		t.start();
		long waitTime = 300000; //five minutes
		while (t.isAlive()){ //DACS search thread is still running
			
			try{ t.join(waitTime);} catch (Exception e){} //wait for waitTime before the next interruption of the DACS search		
			
			float rsBestE = rs.getBestE(); //the best energy from the current partition
			
			//System.out.println("partition "+cObj.partIndex+": curBestE "+bestE+" rsBestE "+rsBestE);
			
			if (rsBestE<otherBestE){ //new best energy for the current partition; update
				
				CommucObj c[] = new CommucObj[1];
				c[0] = new CommucObj();
				c[0].pruningE = rsBestE;
				
				//System.out.println("partition "+cObj.partIndex+": sending update to main node..");
				
				try { MPI.COMM_WORLD.Send(c, 0, 1, MPI.OBJECT, 0, updateTag);} catch (Exception e){}; //send back updated best energy
				
				otherBestE = rsBestE;
			}
			
			//check if there are updates for the best energy from the other partitions
			try {				
				//System.out.println("partition "+cObj.partIndex+": checking for update from main node..");
				
				float c[] = new float[1];
				while (MPI.COMM_WORLD.Iprobe(0, updateTag)!=null) { //new update message received
					
					//System.out.println("partition "+cObj.partIndex+": update from main node received..");
					
					MPI.COMM_WORLD.Recv(c, 0, 1, MPI.FLOAT, 0, updateTag);
					
					//System.out.println("partition "+cObj.partIndex+": updateE "+c[0]+" curBestE: "+bestE);
					
					rs.updateBestE(c[0]);
					otherBestE = Math.min(otherBestE,c[0]);				
				}
			}
			catch (Exception e){};
		}		
		
		logPS.flush();
		logPS.close();		
		
		rs = null;
		cObj.prunedRot = null;//smaller object, for sending back		
		
		long stopTime = System.currentTimeMillis();
		cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);
		
		//System.out.println("Partition "+cObj.partIndex+" done, time: "+cObj.elapsedTime/60.0f+" minutes; sending results to main node..");
		
		return cObj;
	}	
	
	//Forms the DACS partitions based on the splitting residue positions
	private OneMutation [] formDACSpartitions(int maxNumPartitions, int initDepth, int indexMap[][], int numPartitions[],
			boolean prunedRotAtRes[], int majorSplitPos[]){
		
		OneMutation mutArray[] = new OneMutation[maxNumPartitions];
		int curInd[] = new int[initDepth];
		int curMut[] = new int[1];
		curMut[0] = 0;
		
		formDACSpartitionsHelper(initDepth, indexMap, numPartitions, prunedRotAtRes, mutArray, curMut, 0, curInd);
		
		OneMutation tmpArray[] = new OneMutation[curMut[0]]; //trim the size of the partition array
		System.arraycopy(mutArray, 0, tmpArray, 0, curMut[0]);
		mutArray = tmpArray;
		
		System.out.print("Partitioning residues: ");
		for (int i=0; i<initDepth; i++)
			System.out.print(majorSplitPos[i]+" ");
		System.out.println();
		System.out.println("Number of non-zero partitions: "+curMut[0]);
		
		return mutArray;
	}
	
	//Determines all non-pruned partitions for DACS deistribution
	//Called by formDACSpartitions(.)
	private void formDACSpartitionsHelper(int initDepth, int indexMap[][], int numPartitions[],
			boolean prunedRotAtRes[], OneMutation mutArray[], int curMut[], int curDepth, int curInd[]){
		
		if (curDepth>=initDepth){ //new partition
			mutArray[curMut[0]] = new OneMutation();
			mutArray[curMut[0]].mutNum = curMut[0];
			mutArray[curMut[0]].resMut = new int[initDepth];
			for (int i=0; i<initDepth; i++)
				mutArray[curMut[0]].resMut[i] = curInd[i];
			curMut[0]++;
		}
		else {
			for (int i=0; i<numPartitions[curDepth]; i++){
				if (!prunedRotAtRes[indexMap[curDepth][i]]){ //only consider non-pruned partitions for distribution
					curInd[curDepth] = indexMap[curDepth][i];
					formDACSpartitionsHelper(initDepth, indexMap, numPartitions, prunedRotAtRes, mutArray, curMut, curDepth+1, curInd);
				}
			}
		}
	}
	
	//Sorts the DACS partitions by lower energy bounds;
	//The sorted array is returned in mutArray[]
	private void sortDACSpartitions(OneMutation mutArray[], int initDepth, int majorSplitPos[], int numPartitions[], 
			boolean prunedRotAtRes[], int indexMap[][], int numInAS, int numLigRotamers, int residueMap[],
			float pruningE, float initEw, RotamerSearch rsP) {
		
		RotamerSearch rs = rsP; //no changes should be made to the original RotamerSearch object
		rsP = null;
		
		System.out.print("Computing a lower bound on the conformational energy for each partition..");
		for (int m=0; m<mutArray.length; m++){ //for each partition
			
			boolean prunedForPartition[] = new boolean[prunedRotAtRes.length]; //no changes should be made to prunedRotAtRes[]
			System.arraycopy(prunedRotAtRes, 0, prunedForPartition, 0, prunedRotAtRes.length);
			
			mutArray[m].setSortScores(true); //sort by lower bounds
			
			//artificially set to pruned all rotamers at the splitting position, other than the current partitioning rotamer for the given partitioning position
			for (int i=0; i<initDepth; i++){ //first, set all other rotamers for the partitioning positions to pruned
				int curPart = mutArray[m].resMut[i];
				for (int j=0; j<numPartitions[i]; j++){
					int curInd = indexMap[i][j];					
					if (curInd!=curPart) { //not the rotamer for the current partition
						if (!prunedForPartition[curInd]) //rotamer not already pruned
							prunedForPartition[curInd] = true;
					}
				}
			}
			
			//compute a lower bound on the conformational energies for this partition
			rs.DoMinBounds(numInAS,totalNumRotamers,numLigRotamers,residueMap,rotamerIndexOffset,pruningE,prunedForPartition,initEw, false, false, true);
			mutArray[m].score = new BigDecimal(rs.getBoundForPartition());
		}
		System.out.println("done");
		
		//sort the partitions
		System.out.print("Sorting partitions by their lower energy bounds..");
		RyanQuickSort rqs = new RyanQuickSort();
		rqs.Sort(mutArray);
		rqs = null;
		System.out.println("done");
	}
///////////////////////////////////////////////////////////////////////////
//	End of Distributed DACS section
///////////////////////////////////////////////////////////////////////////
	
	
///////////////////////////////////////////////////////////////////////////////////////////////
//	Distributed DEE section
// WARNING: the Distributed DEE section is outdated and must be carefully checked before called
///////////////////////////////////////////////////////////////////////////////////////////////
	//Handles the distribution of the DEE computation to the slave nodes
	private void doDistrDEEMaster(int numInAS, int numLigRotamers, int residueMap[],
			String ligType, ParamSet sParams, boolean ligPresent, String resDefault[],
			String minPEM, String maxPEM, float initEw, boolean prunedRotAtRes[],
			boolean doMinimize, RotamerSearch rs, String sfFile, String rotFile, boolean useSF, 
			int numSpPos, int typeDEE, boolean minimizeBB, boolean scaleInt, float maxIntScale, boolean useEref, float eRef[]){
		
		//the total number of residues (active site + ligand, if present)
		int totalNumRes = numInAS;
		if (ligPresent) //ligand is present
			totalNumRes++;
		
		// Generate all combinations
		int numMutRes = 1; //singles
		if (typeDEE==optPairs)
			numMutRes = 2; //pairs criterion
		int numComb = factorial(totalNumRes).divide(factorial(totalNumRes-numMutRes).multiply(factorial(numMutRes))).intValue();
		int residueMutatable[][] = new int[numComb][totalNumRes];
		generateCombinations(residueMutatable,totalNumRes,numMutRes);
		
		OneMutation mutArray[] = new OneMutation[numComb];
		for (int curMut=0; curMut<mutArray.length; curMut++){
			mutArray[curMut] = new OneMutation();
			mutArray[curMut].mutNum = curMut;
			mutArray[curMut].resMut = new int[totalNumRes];
			for (int curRes=0; curRes<totalNumRes; curRes++)
				mutArray[curMut].resMut[curRes] = residueMutatable[curMut][curRes];
		}
		
		boolean splitFlags[][] = null;
		splitFlags = rs.getSplitFlags(splitFlags);//get the current DE pairs
		outputObject(splitFlags,sfFile);
		
		String AAallowed[] = new String[numInAS]; //the AA's to which each residue can mutate
		for (int i=0; i<numInAS; i++)
			AAallowed[i] = (String)sParams.getValue("RESALLOWED"+i);
		
		MutationManager mutMan = new MutationManager(null,mutArray,false);
		mutMan.setResidueMap(residueMap);
		mutMan.setResDefault(resDefault);
		mutMan.setRotamerIndexOffset(rotamerIndexOffset);
		mutMan.setLigType(ligType);
		mutMan.setarpFilenameMin(minPEM);
		mutMan.setarpFilenameMax(maxPEM);
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
		mutMan.setNumInAS(numInAS);
		mutMan.setnumLigRotamers(numLigRotamers);
		mutMan.numTotalRotamers(totalNumRotamers);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setCalculateVolumes(false);
		mutMan.setInitEw(initEw);
		mutMan.setLigPresent(ligPresent);
		mutMan.setUseSF(useSF);
		mutMan.setSpFlags(splitFlags);
		mutMan.setPrunedRot(prunedRotAtRes);
		mutMan.setSfFile(sfFile);
		mutMan.setRotFile(rotFile);
		mutMan.setDistrDACS(false);
		mutMan.setDistrDEE(true);
		mutMan.setAAallowed(AAallowed);
		mutMan.setNumSpPos(numSpPos);
		mutMan.setTypeDEE(typeDEE);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setScaleInt(scaleInt);
		mutMan.setMaxIntScale(maxIntScale);
		mutMan.setUseEref(useEref);
		mutMan.setEref(eRef);

		try{
			handleDoMPIMaster(mutMan,mutArray.length);
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}
	}
	
	// Distributed DEE Slave function
	private CommucObj doDistrDEESlave(CommucObj cObj) {
				
		long startTime = System.currentTimeMillis();
						
		boolean ligPresent = cObj.ligPresent;

		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,cObj.params,ligPresent,cObj.ligType);
		
		RotamerSearch rs = new RotamerSearch(m,sysStrNum,ligStrNum, hElect, hVDW, hSteric, true,
					true, cObj.epsilon, cObj.stericThresh, cObj.softStericThresh, cObj.distDepDielect, cObj.dielectConst, cObj.doDihedE, cObj.doSolvationE, cObj.solvScale, cObj.vdwMult, rl, grl);
		
		System.out.print("Loading precomputed energy matrix...");

		rs.loadPairwiseEnergyMatrices(new String(cObj.arpFilenameMin+".dat"),true);
		if (cObj.doMinimization)
			rs.loadPairwiseEnergyMatrices(new String(cObj.arpFilenameMax+".dat"),false);
		System.out.println("done");
		
		if (cObj.useEref) { //add the reference energies to the min (and max) intra-energies
			rs.addEref(cObj.eRef, cObj.doMinimization, ligPresent, cObj.numInAS, cObj.residueMap);
		}
		
		
		int totalNumRes = cObj.numInAS;
		if (ligPresent)
			totalNumRes++;				
		
		//determine the two residues in the pair (for pairs DE) or the one residue (split-DEE)
		boolean resInMut[] = new boolean[totalNumRes];
		for (int i=0; i<totalNumRes; i++){
			if (cObj.resMut[i]==1)
				resInMut[i] = true;
			else
				resInMut[i] = false;
		}
		
		
		//Set the allowable AA for each residue;
		// 		the ligand allowable set in the RotamerSearch() constructor	
		boolean addWT = (new Boolean((String)cObj.params.getValue("ADDWT"))).booleanValue();
		for (int j=0; j<cObj.numInAS; j++){
			setAllowablesHelper(rs, cObj.params, addWT, j, cObj.residueMap, cObj.resDefault);
		}
	
		
		//Perform DEE pairs
		boolean splitFlags[][] = (boolean [][])readObject(cObj.sfFileIn);
		rs.setSplitFlags(splitFlags.length);
		rs.setSplitFlags(splitFlags);//initialize the split flags
		
		if (cObj.typeDEE==optPairs){ //simple Goldstein pairs	
			rs.DoDEEPairs(cObj.numInAS, cObj.numTotalRotamers, cObj.numLigRotamers, cObj.residueMap,
					rotamerIndexOffset, cObj.initEw, cObj.prunedRot, resInMut, cObj.doMinimization,
					cObj.useSF, false, true, cObj.minimizeBB, cObj.scaleInt, cObj.maxIntScale);
		}
		else if (cObj.typeDEE==optSplit){ //1- or 2-sp split-DEE
			//Precompute the MinDEE intervals
			rs.doCompMinDEEIntervals(cObj.numInAS, cObj.numTotalRotamers, cObj.numLigRotamers, cObj.residueMap, 
					rotamerIndexOffset, cObj.prunedRot, cObj.scaleInt, cObj.maxIntScale);
			
			cObj.prunedRot = rs.DoDEEConfSplitting(cObj.numInAS, cObj.numTotalRotamers, cObj.numLigRotamers, 
					cObj.residueMap, rotamerIndexOffset, cObj.initEw, cObj.prunedRot, resInMut, 
					cObj.doMinimization, true, cObj.numSpPos, true, cObj.minimizeBB);
		}
		
		long stopTime = System.currentTimeMillis();
		cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);
		
		
		//We cannot send back the whole splitFlags matrix, since it is too big; even just sending the 
		//		newly-identified DE pairs can be infeasible for large systems, so we output the 
		//		new splitFlags matrix to a temp file, which is then read by the master node
		if (cObj.typeDEE==optPairs){
			cObj.sfFileOut = ("tmp_"+cObj.mutationNumber+"_"+stopTime);
			splitFlags = rs.getSplitFlags(splitFlags);
			outputObject(splitFlags,cObj.sfFileOut);
		}
		
		return cObj;
	}	
///////////////////////////////////////////////////////////////////////////
//	End of Distributed DEE section
///////////////////////////////////////////////////////////////////////////
	
	private Object readObject(String inFile){
		return readObject(inFile,true);
	}
	
	private Object readObject(String inFile, boolean repeat){
		Object inObj = null;
		boolean done = false;
		while (!done){
			try{
				ObjectInputStream in = new ObjectInputStream(new FileInputStream(inFile));
				inObj = in.readObject();
				in.close();
				done = true;
			}
			catch (Exception e){
				//System.out.println(e.toString());
				//System.out.println("ERROR: An exception occurred while reading from object file");
				if (repeat)
					done = false;
				else
					done = true;
			}
		}
		return inObj;
	}
	
	private void outputObject(Object outObj, String outFile){
		try{
			FileOutputStream fout = new FileOutputStream(outFile);
			ObjectOutputStream out = new ObjectOutputStream(fout);
			out.writeObject(outObj);
			out.close();
		}
		catch (Exception e){
			System.out.println(e.toString());
			System.out.println("ERROR: An exception occurred while writing object file");
			System.exit(0);
		}
	}
	
	private int countPrunedRot(boolean prunedRot[]){
		int countPruned = 0;
		for (int i=0; i<prunedRot.length; i++){
			if (prunedRot[i])
				countPruned++;
		}
		return countPruned;
	}
	
	private int countPrunedPairs(boolean prunedPairs[][]){
		int countPruned = 0;
		for (int i=0; i<prunedPairs.length; i++){
			for (int j=i+1; j<prunedPairs.length; j++){
				if (prunedPairs[i][j])
					countPruned++;
			}
		}
		return countPruned;
	}
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MPI section
//		The tutorial "One-step Tutorial: MPI: It's easy to get started" 
//			(http://www.lam-mpi.org/tutorials/one-step/ezstart.php ; accessed Oct 23, 2006) was used as MPI code reference
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Setup and do MPI
	//	The parameter is only used by the master node
	public void handleDoMPI(String args[]) throws MPIException {
		
		MPI.Init(args);
		
		int procRank = MPI.COMM_WORLD.Rank();
		numProc = MPI.COMM_WORLD.Size();	
		
		//System.out.println("Node rank: "+procRank+" of "+numProc);
		MPI.COMM_WORLD.Barrier();
		
		setConfigPars();
		MPI.COMM_WORLD.Barrier();
		
		if (procRank==0){ //master node
			
			outputProgInfo();
			parse(args);
		}
		else {//slave node
			handleDoMPISlave();
		}
		
		MPI.Finalize();
	}
	
	//Do MPI for the master node
	public void handleDoMPIMaster(MutationManager mutMan, int size) throws MPIException {
		
		CommucObj cObjArray[] = new CommucObj[size];
		int numFinished = 0;
		
		int curMut = 0;
		for (int curProc=1; curProc<numProc; curProc++){ //distribute a single mutation per processor, for all processors
			
			if (curMut<cObjArray.length){ //more mutations to distribute
				
				System.out.println("Retrieving "+curMut+" of "+(cObjArray.length));
				cObjArray[curMut] = mutMan.getNextComObj(curMut);
				
				MPI.COMM_WORLD.Send(cObjArray, curMut, 1, MPI.OBJECT, curProc, regTag);
				curMut++;
				
				System.out.println("Sent to proc "+curProc);
				System.out.println();
			}
			else
				break;
		}
		
		boolean distrDACS = mutMan.getDistrDACS(); //distributed DACS computation
		
		while (numFinished<cObjArray.length){ //distribute and receive all remaining mutations
			
			CommucObj cObj[] = new CommucObj[1];
			cObj[0] = new CommucObj();
			
			//System.out.println("Receiving message on main node..");
			
			Status s = MPI.COMM_WORLD.Recv(cObj, 0, 1, MPI.OBJECT, MPI.ANY_SOURCE, MPI.ANY_TAG);
			
			//System.out.println("Received message on main node: tag "+s.tag+" source "+s.source);
			
			if (distrDACS){ //DACS computation, so check if the new energy is better than the best energy so far
				
				float curBestE = mutMan.getPruningE();
				
				if (cObj[0].pruningE<curBestE){ //the new energy is better
					
					//System.out.println("Updating best energy on main node from "+curBestE+" to "+cObj[0].pruningE+", source (partition): "+cObj[0].partIndex);
					
					mutMan.setBestScore(new BigDecimal(cObj[0].pruningE));
					mutMan.setPruningE(cObj[0].pruningE);
					
					float c[] = new float[1];
					c[0] = cObj[0].pruningE;
					for (int curProc=1; curProc<numProc; curProc++){ //update the best energy in each partition
						MPI.COMM_WORLD.Isend(c, 0, 1, MPI.FLOAT, curProc, updateTag);
					}
				}
			}
			
			if (s.tag==regTag){ //completed job
				mutMan.processFinishedMutation(cObj[0]);
				numFinished++;
				
				System.out.println("Finished: "+cObj[0].mutationNumber+", Time: "+(cObj[0].elapsedTime/60.0));
				
				if (curMut<cObjArray.length){
					
					System.out.print("Retrieving "+curMut+" of "+(cObjArray.length));
					cObjArray[curMut] = mutMan.getNextComObj(curMut);
					
					MPI.COMM_WORLD.Send(cObjArray, curMut, 1, MPI.OBJECT, s.source, regTag);
					curMut++;
					
					System.out.println(", Sent to proc "+s.source);
					System.out.println();
				}
			}
		}
	}
	
	//Do MPI for a slave node
	public void handleDoMPISlave() throws MPIException {
		
		int rank = MPI.COMM_WORLD.Rank();
		
		while (true){
			
			Status s = MPI.COMM_WORLD.Probe(0, MPI.ANY_TAG);
			if (s!=null) {
				if (s.tag==regTag){ //new computation			
					CommucObj cObj[] = new CommucObj[1];
					cObj[0] = new CommucObj();
					
					//System.out.println("node "+rank+" receiving message from main node..");
					
					MPI.COMM_WORLD.Recv(cObj, 0, 1, MPI.OBJECT, 0, MPI.ANY_TAG);
					
					//System.out.println("node "+rank+" received message from main node..");
					
					if (cObj[0]==null) //computation is done
						return;
					
					cObj[0] = handleKSSlave(cObj[0]); //perform computation
					MPI.COMM_WORLD.Send(cObj, 0, 1, MPI.OBJECT, 0, regTag); //send back result
				}
				else { //(s.tag==updateTag), so discard
					float c[] = new float[1];
					MPI.COMM_WORLD.Recv(c, 0, 1, MPI.FLOAT, 0, MPI.ANY_TAG);
				}
			}
		}
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	 End of MPI section
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
///////////////////////////////////////////////////////////////////////////
//	Backbone Flexibility Section
///////////////////////////////////////////////////////////////////////////
	public void generateBackbones(String s){
		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Backbone config filename (string)

		System.out.println("Performing Backbone Generation");

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters
		
		// Pull search parameters
		int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();	
		String runName = ((String)sParams.getValue("RUNNAME"));
		boolean sysSampling = (new Boolean((String)sParams.getValue("SYSSAMPLING"))).booleanValue();
		double theta = (new Float((String)sParams.getValue("THETA"))).floatValue();
		double alpha = (new Float((String)sParams.getValue("ALPHA"))).floatValue();
		int numSamples = (new Integer((String)sParams.getValue("NUMSAMPLES"))).intValue();
		boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = (String)(sParams.getValue("LIGTYPE"));
		
		if (theta%alpha!=0){
			System.out.println("ERROR: Choose theta = k*alpha, for k - an integer.");
			System.exit(1);
		}
			
		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,ligPresent,ligType);
					
		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.getValue("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
			residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
		}
		System.out.println();
		
		RotamerSearch rs = new RotamerSearch(m,sysStrNum,ligStrNum, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, rl, grl);
		
		rs.doGenBackbones(runName, numInAS, residueMap, theta, alpha, numSamples, sysSampling);
	}
///////////////////////////////////////////////////////////////////////////
//	End of Backbone Flexibility Section
///////////////////////////////////////////////////////////////////////////
	
///////////////////////////////////////////////////////////////////////////
//	Generate Random Conformations Section
///////////////////////////////////////////////////////////////////////////
	//Generates a random set of mutations/conformations for a given system
	public void generateRandConfs(String s) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Run name (for output files)
		// 3: Ligand is present (boolean)
		// 4: Ligand type (if present)
		// 5: Number of conformations to be generated

		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters

		// Pull search parameters
		int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
		String runName = getToken(s,3);
		boolean ligPresent = (new Boolean(getToken(s,4))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = getToken(s,5);
		int num = (new Integer(getToken(s,6))).intValue();
		
		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,ligPresent,ligType);

		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.getValue("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
			residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
		}
		System.out.println();		
		
		PrintStream logPS = setupOutputFile(runName);
		
		Random r = new Random();
		for (int i=0; i<num; i++){
			String AAs[] = new String[numInAS];
			int rot[] = new int[numInAS+1];
			for (int a=0; a<numInAS; a++){
				AAs[a] = rl.getAAName(r.nextInt(numAAallowed));
				int n = rl.getNumRotamers(AAs[a]);
				if (n<=0)
					n = 1;
				rot[a] = r.nextInt(n);
			}
			if (ligPresent){
				int n = grl.getNumRotamers(ligType);
				if (n<=0)
					n = 1;
				rot[numInAS] = r.nextInt(n);
			}
			
			logPS.print(i+" ");
			for (int a=0; a<numInAS; a++)
				logPS.print(AAs[a]+" ");
			if (ligPresent)
				logPS.print(ligType+" ");
			
			for (int a=0; a<numInAS; a++)
				logPS.print(rot[a]+" ");
			if (ligPresent)
				logPS.print(rot[numInAS]+" ");
			
			logPS.println();
			logPS.flush();
		}
		logPS.close();
	}
///////////////////////////////////////////////////////////////////////////
//	End of Generate Random Conformations Section
///////////////////////////////////////////////////////////////////////////
	
	public void doSinglePairE(String s, ParamSet sParams) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutation search parameter filename (string)		
		
		//Only read system and mutation files if sParams is null
		
		if (sParams==null){ //parameter files not read yet, so read them in
			sParams = new ParamSet();
			sParams.addParamsFromFile(getToken(s,2)); //read system parameters
			sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters
		}
		
		int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
		int numMutations = 2; //pairwise energies are computed
		boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB"))).booleanValue();
		boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS"))).booleanValue();
		String backrubFile = (String)sParams.getValue("BACKRUBFILE");
		String runName = (String)sParams.getValue("RUNNAME");
		String minEMatrixName = (String)sParams.getValue("MINENERGYMATRIXNAME");
		String maxEMatrixName = (String)sParams.getValue("MAXENERGYMATRIXNAME");
		
		boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
		boolean inputSysWithLig = ligPresent;
		String ligType = null;
		if (ligPresent)
			ligType = (String)sParams.getValue("LIGTYPE");

		boolean resumeSearch = (new Boolean((String)sParams.getValue("RESUMESEARCH"))).booleanValue();
		String resumeFilename = (String)sParams.getValue("RESUMEFILENAME");
		
		System.out.println("Run Name: "+runName);
		System.out.println("Precomputed Minimum Energy Matrix: "+minEMatrixName);
		System.out.println("Precomputed Maximum Energy Matrix: "+maxEMatrixName);
		System.out.println("Ligand Type: "+ligType);
		System.out.println("Num Residues Allowed to Mutate: "+numMutations);

		
		System.out.println("Computing _All_ Rotamer-Rotamer Energies");
		
		System.out.println("Starting minimum and maximum bound energy computation");
		
		if(resumeSearch) {
			System.out.println("** Resuming Search **");
			System.out.println("     resuming from file: "+resumeFilename);
		}
		
		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,ligPresent,ligType);
					
		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.getValue("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
			residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
		}
		System.out.println();
		
		RotamerSearch rs = new RotamerSearch(m,sysStrNum,ligStrNum, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, rl, grl);
		
		
		
		int resMut[] = new int[numInAS];
		for (int i=0; i<numInAS; i++)
			resMut[i] = 0;
		
		String flagMutType = "AS-AS";
		resMut[0] = 1; resMut[1] = 1;
		
		
		
		boolean useLig = ligPresent;
		
		if (((flagMutType.compareTo("AS-AS")==0)||(flagMutType.compareTo("SHL-AS")==0)||(flagMutType.compareTo("TEMPL")==0))&&(ligPresent)){
			useLig = false;
			m.deleteStrand(1);//we do not need the ligand for these runs
			rs = new RotamerSearch(m, sysStrNum, -1, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, rl, grl);					
		}
		
		System.out.println("Beginning setAllowables");
		// Ligand allowable set in the RotamerSearch() constructor
		boolean addWT = (new Boolean((String)sParams.getValue("ADDWT"))).booleanValue();
		for(int j=0; j<numInAS; j++) {
			if (resMut[j] == 1) {
				setAllowablesHelper(rs, sParams, addWT, j, residueMap, resDefault);
			}
		}
		
		
		boolean shellRun = false;
		boolean intraRun = false;
		boolean templateOnly = false;		
		
		if (flagMutType.compareTo("TEMPL")==0){
			
			// **** Normal Active Site residue runs ****
			// Computes active site residue to active site residue pair energies					
			shellRun = true;ligPresent = false;intraRun = false;templateOnly = true;
		}
		else if (flagMutType.compareTo("AS-AS")==0){
			
			// **** Normal Active Site residue runs ****
			// Computes active site residue to active site residue pair energies					
			shellRun = false;ligPresent = false;intraRun = false;templateOnly = false;
		}	
		else if (flagMutType.compareTo("SHL-AS")==0){
			
			// Then shell runs for the active site residues
			// Computes the active site residue rotamers to shell energies					
			shellRun = true;ligPresent = false;intraRun = false;templateOnly = false;				
		}
		else if (flagMutType.compareTo("INTRA")==0){
			
			// Compute all intra-residue energies					
			shellRun = false;intraRun = true;templateOnly = false;
		}				
		else if (flagMutType.compareTo("LIG-AS")==0){
			
			// **** Ligand present runs ****
			// This section computes the inter-residue energies between
			//  active site residues and the ligand
			shellRun = false;intraRun = false;templateOnly = false;			
		}
		else { //(cObj.flagMutType.compareTo("LIG-SHL")==0)
			
			// Computes ligand rotamer to shell energies
			shellRun = true; intraRun = false;templateOnly = false;
		}
		
		// The goal is that the total energy of a system can be bounded by the sum of 
		//  all pairwise active site residue entries plus the entry for each active site
		//  residue's shell run plus each active site residue's self intra-residue energy.
		//  If a ligand is present then one should add the ligand to shell energy, the
		//  ligand to each active site residue pairwise energy, and the ligand self intra-
		//  residue energy.
		
		//initialize the pairwise energy matrices (partial initialization - only for the residues involved in this computation, e.g., AS-AS)
		PEMHandler pemH = new PEMHandler();
		float minEmatrix[][][][][][] = pemH.initializePairEMatrix(numInAS,inputSysWithLig,resMut,residueMap,
				rs,useLig,ligType,shellRun,intraRun,rl,grl,numAAallowed,false);
		float maxEmatrix[][][][][][] = pemH.copyMultiDimArray(minEmatrix);
		
		//Compute the corresponding matrix entries
		rs.simplePairwiseMutationAllRotamerSearch(residueMap,numInAS,rotamerIndexOffset,
				totalNumRotamers,true,ligPresent,shellRun,intraRun,
				resMut,minEmatrix,maxEmatrix,minimizeBB,doBackrubs,templateOnly,backrubFile);
		
		return;
	}
	
	//Computes conformation energies for different combinations of the energy function parameters
	private void fitEparams(String s){
		
		String firstParam = getToken(s,1);
		
		String sysFile = getToken(s,2);

		// Pull search parameters
		String confResFile = getToken(s,3);
		String runName = getToken(s,4);
		boolean ligPresent = (new Boolean(getToken(s,5))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = getToken(s,6);
		int numResults = (new Integer(getToken(s,7))).intValue();
		boolean minimizeBB = (new Boolean(getToken(s,8))).booleanValue();
		
		int numSteps = 10;
		double maxVdwMult = 1.05;
		double maxSolvScale = 0.3;
		
		solvScale = 0.0;
		double initVdwMult = 0.63;
				
		double vdwDelta = (maxVdwMult-initVdwMult)/numSteps;
		double solvDelta = (maxSolvScale-solvScale)/numSteps;
		
		for (int i1=0; i1<numSteps; i1++){
			solvScale += solvDelta;
			for (int i2=0; i2<2; i2++){
				if (i2==0)
					distDepDielect = true;
				else
					distDepDielect = false;
				
				for (int i3=0; i3<=numSteps; i3++){
					if (i3==0)
						dielectConst = 1.0;
					else
						dielectConst = 4*i3;
					
					softvdwMultiplier = initVdwMult-vdwDelta;
					for (int i4=0; i4<=numSteps; i4++){
						softvdwMultiplier += vdwDelta;
						
						String runNameParams = (runName+"_"+solvScale+"_"+distDepDielect+"_"+dielectConst+"_"+softvdwMultiplier);
						
						String s1 = (firstParam+" "+sysFile+" "+confResFile+" "+runNameParams+" false none "+numResults+" "+minimizeBB);
						handleMinDEEApplyRot(s1);
						
						if (ligPresent){
							runNameParams = (runNameParams+"_lig");
							s1 = (firstParam+" "+sysFile+" "+(confResFile+"_lig")+" "+runNameParams+" true "+ligType+" "+numResults+" "+minimizeBB);
							handleMinDEEApplyRot(s1);
						}
					}
				}
			}
		}
	}
	
	
////////////////////////////////////////////////////////////////
// Compute Residue Entropy Section
////////////////////////////////////////////////////////////////
	/**
	 * Handles the SCMF residue entropy computation.
	 */
	public void handleDoResEntropy(String s, ParamSet sParams) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutation search parameter filename (string)		
		
		//Only read system and mutation files if sParams is null
		
		if (sParams==null){ //parameter files not read yet, so read them in
			sParams = new ParamSet();
			sParams.addParamsFromFile(getToken(s,2)); //read system parameters
			sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters
		}
		
		String runName = (String)sParams.getValue("RUNNAME");
		
		String matrixName = (String)sParams.getValue("MATRIXNAME");
		boolean useEref = (new Boolean((String)sParams.getValue("USEEREF"))).booleanValue();
		float dist = (new Float((String)sParams.getValue("DIST"))).floatValue();
		String rotProbFile = (String)sParams.getValue("ROTPROBFILE");
		float stericE = (new Float((String)sParams.getValue("STERICE"))).floatValue();
		float maxPairE = (new Float((String)sParams.getValue("MAXPAIRE"))).floatValue();
		
		
		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,false,null);
		
		int numRes = m.strand[sysStrNum].numberOfResidues;
		
		String resDefault[] = new String[numRes];
		int defResNum[] = new int[numRes];
		for(int i=0;i<numRes;i++){
			resDefault[i] = m.strand[sysStrNum].residue[i].name;
			defResNum[i] = m.strand[sysStrNum].residue[i].getResNumber();
		}
		
		
		rotProbFile = (rotProbFile+".dat");
		float rotProb[][] = (float [][])readObject(rotProbFile,false);
		if (rotProb==null){ //Perform SCMF to compute the rotamer probabilities		
			
			//read in or compute all of the energy matrices;
			//NOTE: backbone-to-backbone and rotamer-to-backbone energies are included in the pairwise rotamer energies
			float asasE[][][][] = getResEntropyEmatricesPair(matrixName, sParams, numRes, resDefault, m, runName, dist);
			float intraEnergies[][] = getResEntropyEmatricesIntra(matrixName, sParams, numRes, m, runName);
			float eRef[] = getResEntropyEmatricesEref(useEref, null, null, null, intraEnergies, numRes);
			
			rotProb = compRotProbSCMF(numRes, intraEnergies, asasE, eRef, rotProbFile, resDefault, stericE, maxPairE);
		}
		
		int numProx[] = new int[numRes]; //get the number of proximate residues for each residue position
		for (int i=0; i<numProx.length; i++)
			numProx[i] = 0;
		String asDistFile = matrixName+"_dist.dat";
		boolean as[][] = (boolean [][])readObject(asDistFile,false);
		for (int i=0; i<numRes; i++){
			for (int j=i+1; j<numRes; j++){
				if (as[i][j]){
					numProx[i]++;
					numProx[j]++;
				}
			}
		}
		
		m = null;		
		
		PrintStream logPS = setupOutputFile(runName);
		
		logPS.print("resNum pdbResNum resDefault entropy"+" ");
		for (int j=0; j<resAllowed.length; j++)
			logPS.print(resAllowed[j]+" ");
		logPS.println("numProx");
		logPS.flush();
		
		
		//Compute the AA probabilities for each residue position (as a function of the rotamer probabilities);
		//Compute the entropy at each position as a function of the amino acid probabilities for that position
		final double kB = 1.0;
		for (int i=0; i<numRes; i++){
			
			if ( !resDefault[i].equalsIgnoreCase("PRO") ){
			
				float aaProbBBE[] = new float[resAllowed.length];
				int curInd = 0;
				for (int j=0; j<aaProbBBE.length; j++){
					
					int numCurRot = rl.getNumRotamers(resAllowed[j]);
					if (numCurRot==0)
						numCurRot = 1;
					
					aaProbBBE[j] = 0.0f;
					for (int r=0; r<numCurRot; r++){
						aaProbBBE[j] += rotProb[i][curInd];
						curInd++;
					}
				}
				
				//Compute the unnormalized AA probabilities as a weighted function of energies and PDB stats (if included)
				double aaProbUnNorm[] = new double[resAllowed.length];
				double aaNorm = 0.0;
				for (int j=0; j<resAllowed.length; j++){
					
					aaProbUnNorm[j] = 1.0;
					aaProbUnNorm[j] *= aaProbBBE[j];
					
					aaNorm += aaProbUnNorm[j];
				}
				
				//Normalize the probabilities
				double aaProb[] = new double[resAllowed.length];
				for (int j=0; j<resAllowed.length; j++){
					if (aaNorm!=0.0)
						aaProb[j] = aaProbUnNorm[j]/aaNorm;
					else
						aaProb[j] = 0.0;
				}
				
				//Compute the entropy for the current residue position
				double sumAA = 0.0;
				for (int j=0; j<aaProb.length; j++){
					if (aaProb[j]>0.0)
						sumAA += aaProb[j]*Math.log(aaProb[j]);
				}
				
				double entropy = -kB * sumAA;
				
				logPS.print(i+" "+defResNum[i]+" "+resDefault[i]+" "+entropy+" ");
				for (int j=0; j<aaProb.length; j++)
					logPS.print(aaProb[j]+" ");
				
				logPS.println(numProx[i]);
				logPS.flush();
			}
			else {
				logPS.println(i+" "+defResNum[i]+" "+resDefault[i]+" "+0.0); //only for residue positions with wildtype Pro
			}
		}
		logPS.close();
	}
	
	//Computes the rotamer probabilities for all rotamers at all residue positions using SCMF
	private float[][] compRotProbSCMF(int numRes, float intraEnergies[][],
			float asasE[][][][], float eRef[], String rotProbFile, String resDefault[], float stericE, float maxPairE){
		
		final float constR = (float)(1.9891/1000.0);//the gas constant
		float T = 50000; //initial temperature
		final float endT = 298.15f; //the minimum temperature for annealing
		float tStepSize = 100.0f; //the temperature step size for annealing
		final float eps = 0.0001f; //the convergence threshold
		final float lambda = 0.5f; //scaling factor for updating the rotamer probabilities
		
		
		for (int i=0; i<asasE.length; i++){//Set the max energy for any element in asasE[][][][] to maxPairE
			if (asasE[i]!=null){
				for (int j=0; j<asasE[i].length; j++){
					if (asasE[i][j]!=null){
						for (int k=0; k<asasE[i][j].length; k++){
							if (asasE[i][j][k]!=null){
								for (int l=0; l<asasE[i][j][k].length; l++){
									if (asasE[i][j][k][l]>maxPairE)
										asasE[i][j][k][l] = maxPairE;
								}
							}
						}
					}
				}
			}
		}
		
		int numPrunedRot = 0;
		boolean prunedRot[][] = new boolean[numRes][totalNumRotamers];
		for (int i=0; i<numRes; i++){
			for (int j=0; j<totalNumRotamers; j++){
				if ( (intraEnergies[1+i*totalNumRotamers+j][0]) > stericE){
					prunedRot[i][j] = true;
					numPrunedRot++;
				}
				else
					prunedRot[i][j] = false;
			}
		}
		System.out.println("Num rotamers pruned due to incompatibility with the template: "+numPrunedRot);
		
		
		//For each residue, compute the probability of each rotamer for that residue
		float Emf[][] = new float[numRes][totalNumRotamers];
		float rotProb[][] = new float[numRes][totalNumRotamers];
		float oldProb[][] = new float[numRes][totalNumRotamers];
		for (int i=0; i<numRes; i++){
			for (int j=0; j<totalNumRotamers; j++) {
				if (!prunedRot[i][j])
					rotProb[i][j] = 1.0f/totalNumRotamers;
				else
					rotProb[i][j] = 0.0f;
				
				oldProb[i][j] = rotProb[i][j];
			}
		}
		
		while (T>=endT){ //perform annealing
			
			System.out.println("Starting run at T = "+T);
			
			boolean done = false;
			while (!done){
				
				//Compute the new mean-field energy for each rotamer
				for (int i=0; i<numRes; i++){
					
					if ( !resDefault[i].equalsIgnoreCase("PRO") ){
					
						for (int j=0; j<totalNumRotamers; j++){
							
							if (!prunedRot[i][j]){
								
								Emf[i][j] = intraEnergies[1+i*totalNumRotamers+j][0] - eRef[getAAindFromRotNum(j)];
								
								if (asasE[i]!=null){
									
									for (int k=0; k<numRes; k++){
										if ( (k!=i) && (!resDefault[k].equalsIgnoreCase("PRO")) ){ //for all residues with which i_j has contact
											if ( (i<k) && (asasE[i][k]!=null) ){
												for (int l=0; l<totalNumRotamers; l++){
													if (!prunedRot[k][l])
														Emf[i][j] += asasE[i][k][j][l]*rotProb[k][l];
												}
											}
											else if ( (i>k) && (asasE[k][i]!=null) ){
												for (int l=0; l<totalNumRotamers; l++){
													if (!prunedRot[k][l])
														Emf[i][j] += asasE[k][i][l][j]*rotProb[k][l];
												}
											}
										}
									}
								}
							}
						}
					}
				}
				
				//Update the rotamer probabilities			
				for (int i=0; i<numRes; i++){
					
					if ( !resDefault[i].equalsIgnoreCase("PRO") ){
					
						float normFactor = 0.0f;
						for (int j=0; j<totalNumRotamers; j++){
							
							if (!prunedRot[i][j])
								normFactor += (float)Math.exp( -Emf[i][j] / (constR*T));	
						}
						
						for (int j=0; j<totalNumRotamers; j++){
							
							if (!prunedRot[i][j]){
								
								oldProb[i][j] = rotProb[i][j]; //the probability before the update
								
								if (normFactor!=0.0f)
									rotProb[i][j] = lambda*((float)Math.exp( -Emf[i][j] / (constR*T)) / normFactor) + (1-lambda)*oldProb[i][j];
								else
									rotProb[i][j] = 0.0f;
							}
						}
					}
				}
				
				float rms = checkRotProbConvEntropy(rotProb,oldProb,resDefault,prunedRot);
				
				if (rms>eps)
					done = false;
				else
					done = true;
			}
			
			T -= tStepSize;
		}
		
		outputObject(rotProb,rotProbFile);
		
		return rotProb;
	}
	
	//Checks if the rotamer probabilities for the entropy computation have converged
	private float checkRotProbConvEntropy(float rotProb[][], float oldProb[][], String resDefault[], boolean prunedRot[][]){
		
		float sum = 0.0f;
		for (int i=0; i<rotProb.length; i++){
			if ( !resDefault[i].equalsIgnoreCase("PRO") ){
				for (int j=0; j<rotProb[i].length; j++){
					if (!prunedRot[i][j])
						sum += (float)Math.pow( (rotProb[i][j]-oldProb[i][j]) , 2.0);
				}
			}
		}
		float rms = (float)Math.sqrt(sum);
		
		System.out.println("RMS: "+rms);
		
		return rms;
	}
	
	//Reads in (if computed) or computes the energy matrices for rot-to-rot pairwise energies;
	//This is not setup to handle energy minimization properly (mainly due to the approach used for computing only a small
	//		subset of the pairwise rot-to-rot energies), so doMinimize and minimizeBB should be false
	private float [][][][] getResEntropyEmatricesPair(String matrixName, ParamSet sParams, int numRes, String resDefault[], Molecule m,
			String runName, float dist){
		
		String origPDB = (String)sParams.getValue("PDBNAME");
		
		//Check for the pairwise rot-to-rot file;
		String pairName = matrixName+"_pair.dat";
		
		float asasE[][][][] = readPairMatrixEntropy(pairName,numRes);
		if (asasE==null){ //compute the rot-to-template energy matrix
			
			//For each residue position i, get all residue positions that are within dist
			String asDistFile = matrixName+"_dist.dat";
			boolean as[][] = (boolean [][])readObject(asDistFile,false);
			if (as==null){
				as = new boolean[numRes][numRes];
				computeEntropyEmatrixMaster(numRes,runName+".log",asDistFile,sParams,false,false,false,null,
						null,null,true,dist,as,false); //the distances are returned in as[][]
			}
			
			int numPairs = 0;
			asasE = new float[numRes][numRes][][];
			for (int i=0; i<numRes; i++){
				if (!resDefault[i].equalsIgnoreCase("PRO")){
					for (int j=i+1; j<numRes; j++){
						if (!resDefault[j].equalsIgnoreCase("PRO")){
							if (as[i][j]){
								asasE[i][j] = new float[1][];
								numPairs++;
							}
						}
					}
				}
			}
						
			computeEntropyEmatrixMaster(numRes,runName+".log",pairName,sParams,false,false,false,null,
					null,asasE,false,0.0f,null,false);
			
			sParams.setValue("PDBNAME", origPDB);
			m = new Molecule();
			setupMolSystem(m,sParams,false,null);
			
			if (ligStrNum>=0) //the ligand is not used here
				m.deleteStrand(ligStrNum);
			
			asasE = readPairMatrixEntropy(pairName,numRes);
		}		
		
		return asasE;
	}
	
	//Reads in (if computed) or computes the energy matrices for intra-rot energies;
	//This is not setup to handle energy minimization properly (mainly due to the approach used for computing only a small
	//		subset of the pairwise rot-to-rot energies), so doMinimize and minimizeBB should be false
	private float [][] getResEntropyEmatricesIntra(String matrixName, ParamSet sParams, int numRes, Molecule m, String runName){
		
		String origPDB = (String)sParams.getValue("PDBNAME");
		
		//Check for the intra energies file
		String intraName = matrixName+"_intra.dat";
		
		float intraEnergies[][] = (float [][])readObject(intraName,false); //check if already computed			
		
		if (intraEnergies==null) //compute the intra-rotamer energy matrix
			computeEntropyEmatrixMaster(numRes,runName+".log",intraName,sParams,false,false,false,null,
					intraEnergies,null,false,0.0f,null,true); 
		
		sParams.setValue("PDBNAME", origPDB);
		m = new Molecule();
		setupMolSystem(m,sParams,false,null);
		
		if (ligStrNum>=0) //the ligand is not used here
			m.deleteStrand(ligStrNum);
		
		intraEnergies = (float [][])readObject(intraName,false);
		
		return intraEnergies;
	}
	
	//Reads in the amino acid reference energies (if used);
	private float [] getResEntropyEmatricesEref(boolean useEref, float intraEnergies[][][][][][], StrandRotamers sysLR, int residueMap[], float intraEnergiesEntropy[][], int numRes){
		
		if ( (intraEnergies==null && intraEnergiesEntropy==null) || (intraEnergies!=null && intraEnergiesEntropy!=null) ){ //exactly one of the two matrices should be non-null
			System.out.println("ERROR: exactly one matrix can be used for the reference energy computation.");
			System.exit(1);
		}
		
		float eRef[] = new float[resAllowed.length];
		
		if (useEref){ //use AA reference energies		
			eRef = compEref(intraEnergies,sysLR,residueMap,intraEnergiesEntropy,numRes);
		}
		else {
			for (int i=0; i<eRef.length; i++)
				eRef[i] = 0.0f;
		}
		
		return eRef;
	}
	
	//Reads in the pairwise energy matrix for the entropy computation
	private float [][][][] readPairMatrixEntropy(String fName, int numRes){
		
		BufferedReader bufread = null;
		try {
			File file = new File(fName);
			FileReader fr = new FileReader(file);
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			return null;
		}
		
		float asasE[][][][] = new float[numRes][][][];

		boolean done = false;
		String str = null;
		
		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
			 System.out.println("ERROR: An error occurred while reading input");
			 System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else if (getToken(str,1).charAt(0)=='%') //skip comment lines
				continue;
			else {				
				int i = new Integer(getToken(str,1)).intValue();
				String name = getToken(str,2);
				asasE[i] = (float [][][])readObject(name,false);
				if (asasE[i]==null){
					System.out.println("ERROR: Could not read data from file "+name);
					System.exit(1);
				}
			}
		}
					
		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}
		
		return asasE;
	}
	
	//Determines how many residue positions in the system strand numbered from (pos+1...) (pos is strand-relative numbering)
	//		are within dist from residue position pos
	public boolean [] getProxAS(Molecule m, int pos, float dist, boolean as[]){
		
		int residueMap[] = new int[2];
		residueMap[0] = pos;
		
		for (int i=pos+1; i<m.strand[sysStrNum].numberOfResidues; i++){
			
			if (i!=pos){
				
				boolean done = false;
				
				residueMap[1] = i;
				
				Molecule m1 = getASASMolEntropy(m, residueMap);
				
				StrandRotamers sysLR = new StrandRotamers(rl,m1.strand[0]);
				
				for (int j=0; j<m1.numberOfResidues; j++){
					for (int r=0; r<resAllowed.length; r++){
						sysLR.setAllowable(j,resAllowed[r]);
						m1.residue[j].flexible = true;
					}
				}
				
				for(int q1=0;q1<sysLR.getNumAllowable(0);q1++) {
					
					int AAindex1 = sysLR.getIndexOfNthAllowable(0,q1);
					sysLR.changeResidueType(m1,0,rl.getAAName(AAindex1),true,true);
					
					for(int q2=0;q2<sysLR.getNumAllowable(1);q2++) {
						
						int AAindex2 = sysLR.getIndexOfNthAllowable(1,q2);
						sysLR.changeResidueType(m1,1,rl.getAAName(AAindex2),true,true);
						
						int numRot1 = rl.getNumRotForAAtype(AAindex1);
						
						int w1 = 0;
						if (numRot1<=0)
							w1 = -1;
						
						while ((w1<numRot1)&&(!done)){
							
							if (w1!=-1)
								sysLR.applyRotamer(m1, 0, w1);
							
							int numRot2 = rl.getNumRotForAAtype(AAindex2);
							
							int w2 = 0;
							if (numRot2<=0)
								w2= -1;
							
							while ((w2<numRot2)&&(!done)){
								
								if (w2!=-1)
									sysLR.applyRotamer(m1, 1, w2);
								
								Residue r1 = m1.strand[0].residue[0];
								Residue r2 = m1.strand[0].residue[1];
								
								if (r1.getDist(r2,true)<=dist){
									as[i] = true;
									done = true;
								}
								
								w2++;
							}
							
							w1++;
						}
						if (done)
							break;
					}
					if (done)
						break;
				}
				if (!done)
					as[i] = false;
			}
		}
		
		return as;
	}
	
	//Returns the AA index into rotamerIndexOffset to which rotNum belongs
	private int getAAindFromRotNum(int rotNum){
		for (int i=0; i<rotamerIndexOffset.length-1; i++){
			if ( (rotNum>=rotamerIndexOffset[i]) && (rotNum<rotamerIndexOffset[i+1]) )
				return i;
		}
		if (!(rotNum>=totalNumRotamers))
			return (rotamerIndexOffset.length-1);
		else
			return -1;
	}
	
	//Distributes the different types of energy computation for the entropy calculation
	private void computeEntropyEmatrixMaster(int numRes, String runName, String matrixName, ParamSet sParams, boolean doMinimize, boolean minimizeBB, boolean doBackrubs, String backrubFile,
			float bbEnergies[][], float asasE[][][][], boolean compASASdist, float dist, boolean asDist[][], boolean intraRun){
		
		
		int mutEnerMatrixSize = 0;
		int residueMap[] = null;
		OneMutation mutArray[] = null;
		int numInAS = numRes;
		
		int numMut = 0;
		
		if (compASASdist){ //compute the min distance between any pair of rotamers for each pair of residue positions
			mutArray = new OneMutation[numInAS];
			for (int i=0; i<asDist.length; i++)			
				mutArray[i] = new OneMutation();
			
			numMut = numRes;
		}
		else {
			if (intraRun) { //computing intra energies
				
				System.out.println("Starting intra-rot energy computation..");
				
				mutEnerMatrixSize = 1 + totalNumRotamers*numInAS;
				
				bbEnergies = new float[mutEnerMatrixSize][1];
				for(int i=0; i<mutEnerMatrixSize; i++) {
					for(int j=0; j<1; j++){
						bbEnergies[i][j] = 0.0f;
					}
				}
				
				numInAS = 1;
				residueMap = new int[numInAS];
				
				mutArray = new OneMutation[numRes];
				for (int i=0; i<mutArray.length; i++){
					mutArray[i] = new OneMutation();
					mutArray[i].flagMutType = "INTRA";
				}
				
				numMut = numRes;
			}
			
			else { //AS-AS energies
				
				System.out.println("Starting rot-to-rot energy computation..");
				
				numInAS = 2;
				
				int numPairs = 0;
				for (int i=0; i<asasE.length; i++){
					for (int j=i+1; j<asasE[0].length; j++){
						if (asasE[i][j]!=null)
							numPairs++;
					}
				}
				mutArray = new OneMutation[numPairs];
				
				int curPair = 0;
				for (int i=0; i<asasE.length; i++){
					for (int j=i+1; j<asasE[0].length; j++){
						if (asasE[i][j]!=null){				
							mutArray[curPair] = new OneMutation();
							mutArray[curPair].flagMutType = "AS-AS";
							mutArray[curPair].resMut = new int[2];
							mutArray[curPair].resMut[0] = i; //strand-relative numbering (system strand)
							mutArray[curPair].resMut[1] = j;
							curPair++;
							
							asasE[i][j] = new float[totalNumRotamers][totalNumRotamers];
						}
					}
				}
				
				numMut = numPairs;
			}
		}
		
		
		MutationManager mutMan = new MutationManager(runName,mutArray,true);
		mutMan.setResidueMap(residueMap);
		mutMan.setRotamerIndexOffset(rotamerIndexOffset);
		mutMan.setLigType(null);
		mutMan.setarpFilenameMin(matrixName);
		mutMan.setIntraEntropyMatrixMin(bbEnergies);
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setSoftStericThresh(softStericThresh);
		mutMan.setNumInAS(numInAS);
		mutMan.numTotalRotamers(totalNumRotamers);
		mutMan.numResAllowed(numAAallowed);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setDoBackrubs(doBackrubs);
		mutMan.setBackrubFile(backrubFile);
		mutMan.setCalculateVolumes(false);
		mutMan.setLigPresent(false);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setEntropyComp(true);
		mutMan.setPairEntropyMatrix(asasE);
		mutMan.setASdistMatrix(asDist);
		mutMan.setASdist(dist);
		mutMan.setCompASdist(compASASdist);

		
		try{
			handleDoMPIMaster(mutMan,numMut);
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}
		
		if (compASASdist){
			asDist = mutMan.getASdistMatrix();
			outputObject(asDist,matrixName);
		}
		else {
			if (intraRun){
				bbEnergies = mutMan.getMinEmatrixEntropy();
				
				int numCompEntries = 0;				
				for(int i1=0; i1<mutEnerMatrixSize; i1++){
					if ((bbEnergies[i1][0]!=0.0f))
						numCompEntries++;
				}
				System.out.println("Num computed entries: "+numCompEntries);
				
				outputObject(bbEnergies,matrixName);
			}
			else {
				asasE = mutMan.getPairEntropyEmatrix();		
				
				PrintStream logPS = setupOutputFile(matrixName);
				for (int i=0; i<asasE.length; i++){
					if (asasE[i]!=null){
						String fn = ("peme/pem_entr_"+i);
						logPS.println(i+" "+fn);
						outputObject(asasE[i],fn);
						logPS.flush();
					}
				}
				logPS.close();
			}
		}
	}
	
	private CommucObj handleDoResEntropySlave(CommucObj cObj){
		
		long startTime = System.currentTimeMillis();

		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,cObj.params,false,null); //the ligand is not used here
		
		
		if (cObj.compASdist){ //AS-AS distance computation
			cObj.asDist = getProxAS(m,cObj.mutationNumber,cObj.dist,cObj.asDist);
			cObj.compEE = new SamplingEEntries[0];
		}
		else { //AS-AS or INTRA energy computation
			
			int numInAS = cObj.numInAS;
		
			int resMut[] = new int[numInAS];
			int residueMap[] = new int[cObj.residueMap.length];
			
			boolean shellRun = false; boolean ligPresent = false; boolean intraRun = false; boolean templateOnly = false;
			
			if (cObj.flagMutType.compareTo("AS-AS")==0){ //AS-AS run
				
				for (int i=0; i<numInAS; i++)
					resMut[i] = 1;
				
				m = getASASMolEntropy(m,cObj.residueMap);
				
				residueMap[0] = 0;
				residueMap[1] = 1;
			}
			
			else if (cObj.flagMutType.compareTo("INTRA")==0){
				
				intraRun = true;
				
				m = getASASMolEntropy(m,cObj.residueMap);
				
				residueMap = new int[1];
				residueMap[0] = 0;
				resMut = new int[1];
				resMut[0] = 1;
			}
			
			else {
				System.out.println("ERROR: only AS-AS and INTRA runs allowed for the pairwise entropy matrix precomputation.");
				System.exit(1);
			}
			
			RotamerSearch rs = new RotamerSearch(m, sysStrNum, ligStrNum, hElect, hVDW, hSteric, true,
						true, cObj.epsilon, cObj.stericThresh, cObj.softStericThresh, cObj.distDepDielect, cObj.dielectConst, cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult, rl, grl);
			
			for(int j=0; j<numInAS; j++) {
				for(int q=0;q<resAllowed.length;q++)
					rs.setAllowable(residueMap[j],resAllowed[q]);
			}
			
			//initialize the pairwise energy matrices (partial initialization - only for the residues involved in this computation, e.g., AS-AS)
			PEMHandler pemH = new PEMHandler();
			float minEmatrix[][][][][][] = pemH.initializePairEMatrix(numInAS,ligPresent,resMut,residueMap,rs,ligPresent,null,shellRun,intraRun,rl,grl,numAAallowed,false);
			float maxEmatrix[][][][][][] = pemH.copyMultiDimArray(minEmatrix);			
			
			rs.simplePairwiseMutationAllRotamerSearch(residueMap,numInAS,cObj.rotamerIndexOffset,
					cObj.numTotalRotamers,cObj.doMinimization,ligPresent,shellRun,intraRun,
					resMut,minEmatrix,maxEmatrix,cObj.minimizeBB,cObj.doBackrubs,templateOnly,cObj.backrubFile);			
			
			long stopTime = System.currentTimeMillis();
			cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);
			
			//Store the information in less space to allow the master node to buffer several cObj at once
			cObj.compEE = pemH.generateCompEE(minEmatrix, maxEmatrix);
		}
		
		return cObj;
	}
	
	//Returns a molecule that contains only the residues in the system strand (sysStrNum) of molecule m that are specified by residueMap[];
	//	This function is used for the pairwise energy matrix computation in the residue entropy calculations
	private Molecule getASASMolEntropy (Molecule m, int residueMap[]){
		
		Molecule m1 = new Molecule();
		
		for (int i=0; i<residueMap.length; i++){
			
			Residue oldResidue = m.strand[sysStrNum].residue[residueMap[i]];
			
			Residue newResidue = new Residue();
			newResidue.name = oldResidue.name;
			newResidue.fullName = oldResidue.fullName;
			
			for (int j=0; j<oldResidue.numberOfAtoms; j++){
				
				Atom oldAtom = oldResidue.atom[j];
				
				Atom newAtom = new Atom(oldAtom.name,oldAtom.coord[0],oldAtom.coord[1],oldAtom.coord[2]);
				newAtom.modelAtomNumber = oldAtom.modelAtomNumber;
				newAtom.strandNumber = oldAtom.strandNumber;
				newAtom.elementType = oldAtom.elementType;
				newResidue.addAtom(newAtom);
			}
			
			m1.addResidue(0,newResidue);
		}
		
		//Determine the bonds between the atoms in the molecule
		m1.determineBonds();
		
		// Assign the molecule relative atom numbers
		m1.updateMoleculeAtomNumbers();
		
		m1.strand[0].isProtein = true;
		
		return m1;
	}
	
	//Computes the amino acid reference energies using the intra-rotamer energies from intraEnergies or intraEnergiesEntropy;
	//For each amino acid type, takes the min energy among all rotamers for that amino acid type, for all numRes residues
	private float [] compEref(float intraEnergies[][][][][][], StrandRotamers sysLR, int residueMap[], float intraEnergiesEntropy[][], int numRes){
		
		float bigE = (float)Math.pow(10,38);
		float eRef[] = new float[numAAallowed];
		for (int i=0; i<eRef.length; i++)
			eRef[i] = bigE;
		
		int ind = 1; //skip the entry [0][0], since this is the fixed template energy
		for (int i=0; i<numRes; i++){
			int numAA = numAAallowed;
			if (intraEnergies!=null) { //the six-dimensional, so the energies for only a subset of the amino acid types are available 
				numAA = sysLR.getNumAllowable(residueMap[i]);
			}
			for (int j=0; j<numAA; j++){
				int aaInd = j;
				if (intraEnergies!=null)
					aaInd = sysLR.getIndexOfNthAllowable(residueMap[i],j);
				int numRot = rl.getNumRotForAAtype(aaInd);
				if (numRot==0) //ALA or GLY
					numRot = 1;
				float curMin = bigE;
				for (int k=0; k<numRot; k++){
					if (intraEnergies!=null)
						curMin = Math.min(curMin,intraEnergies[i][aaInd][k][i][0][0]);
					else
						curMin = Math.min(curMin,intraEnergiesEntropy[ind][0]);
					ind++;
				}
				eRef[aaInd] = Math.min(eRef[aaInd],curMin);
			}			
		}
		
		for (int i=0; i<eRef.length; i++){
			if (eRef[i]==bigE)
				eRef[i] = 0.0f;
		}
		return eRef;
	}
//////////////////////////////////////////////////////
// End Compute Residue Entropy Section
//////////////////////////////////////////////////////
	
	private void selectResidues(String s) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Residue search filename (string)		
		
		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters	
		
		String runName = (String)sParams.getValue("RUNNAME");
		int numRes = (new Integer((String)sParams.getValue("NUMRES"))).intValue();
		int pdbRes[] = new int[numRes];
		float dist[] = new float[numRes];
		boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = (String)(sParams.getValue("LIGTYPE"));
		
		String resString = (String)sParams.getValue("RESIDUES");
		String distString = ((String)sParams.getValue("DIST"));
		for (int i=0; i<numRes; i++){
			pdbRes[i] = new Integer((String)getToken(resString,i+1)).intValue();
			dist[i] = new Float((String)getToken(distString,i+1)).floatValue();
		}
		
		
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,ligPresent,ligType);
		
		
		//Map from the pdb residue numbers to the residue index in m.residue[]
		int residues[] = new int[numRes];
		int curRes = 0;
		for (int i=0; i<m.numberOfResidues; i++){
			for (int j=0; j<numRes; j++){
				if (m.residue[i].getResNumber()==pdbRes[j]){
					residues[curRes] = i;
					curRes++;
					break;
				}
			}
		}
		
		
		boolean asProx[] = new boolean[m.numberOfResidues];
		for (int i=0; i<asProx.length; i++)
			asProx[i] = false;
		
		for (int res=0; res<numRes; res++){
			Residue r1 = m.residue[residues[res]];
			
			for (int i=0; i<asProx.length; i++){
				
				if (i!=residues[res]){
					
					Residue r2 = m.residue[i];
					
					if (r1.getDist(r2,true)<=dist[res])
						asProx[i] = true;
				}
				else
					asProx[i] = true;
			}
		}
		
		Molecule m1 = new Molecule(); //add all proximate residues; the connectivity/bonds will not be valid
		for (int i=0; i<asProx.length; i++){
			if (asProx[i]){
				m1.addResidue(0, m.residue[i]);
			}
		}
		saveMolecule(m1,runName+".pdb",0.0f);
	}
	
	//Computes the information necessary to generate a residue interaction graph for the given system;
	//		computes the minimum distance and minimum energy (absolute value) for each residue pair in residueMap[], 
	//		considering all possible unpruned rotamers for the given residues;
	//		ligand interactions are also computed if a ligand is present
	private void genInteractionGraph(int numInAS, RotamerSearch rs, boolean prunedRotAtRes[], boolean ligPresent, String ligType,
			int numLigRotamers, String runName, int residueMap[], float eInteractionCutoff, float distCutoff, Molecule m, 
			boolean usePairSt, float pairSt) {
		
		if (eInteractionCutoff<0.0f) //the cutoff should be non-negative, since we are comparing absolute values of energies against it
			eInteractionCutoff = 0.0f;
		
		float dist[][] = new float[numInAS][numInAS];
		float eInteraction[][] = new float[numInAS][numInAS];		
		for (int i=0; i<numInAS; i++){
			for (int j=0; j<numInAS; j++){
				dist[i][j] = (float)Math.pow(10, 38);
				eInteraction[i][j] = 0.0f;
			}
		}
		
		float ligDist[] = null;
		float ligE[] = null;
		if (ligPresent){
			ligDist = new float[numInAS];
			ligE = new float[numInAS];
			for (int i=0; i<numInAS; i++){
				ligDist[i] = (float)Math.pow(10, 38);
				ligE[i] = 0.0f;
			}
		}
		
		
		int numRes = numInAS;
		if (ligPresent)
			numRes++;
		
		int molResMap[] = new int[numRes];
		for (int i=0; i<numInAS; i++)
			molResMap[i] = m.strand[sysStrNum].residue[residueMap[i]].moleculeResidueNumber;
		
		if (ligPresent)
			molResMap[numInAS] = m.strand[ligStrNum].residue[0].moleculeResidueNumber;
		
		
		for (int i=0; i<numInAS; i++){

			for(int q1=0;q1<rs.sysLR.getNumAllowable(residueMap[i]);q1++) {
				
				int AAindex1 = rs.sysLR.getIndexOfNthAllowable(residueMap[i],q1);
				
				int numRot1 = rl.getNumRotForAAtype(AAindex1);
				if (numRot1==0)
					numRot1 = 1;
				
				for (int r1=0; r1<numRot1; r1++){
					
					if (!prunedRotAtRes[i*totalNumRotamers + rotamerIndexOffset[AAindex1] + r1]){ //rotamer not pruned
				
						for (int j=i+1; j<numInAS; j++){				
							
							for(int q2=0;q2<rs.sysLR.getNumAllowable(residueMap[j]);q2++) {
								
								int AAindex2 = rs.sysLR.getIndexOfNthAllowable(residueMap[j],q2);
								
								int numRot2 = rl.getNumRotForAAtype(AAindex2);
								if (numRot2==0)
									numRot2 = 1;
								
								for (int r2=0; r2<numRot2; r2++){
									
									if (!prunedRotAtRes[j*totalNumRotamers + rotamerIndexOffset[AAindex2] + r2]){
										
										float pairE = rs.getMinMatrix()[i][AAindex1][r1][j][AAindex2][r2];
										
										int smallMolResMap[] = new int[2];
										smallMolResMap[0] = molResMap[i];
										smallMolResMap[1] = molResMap[j];
										
										Molecule m1 = getMolRes(m,smallMolResMap);
										
										StrandRotamers sysLR = new StrandRotamers(rl,m1.strand[sysStrNum]);
										sysLR.setAllowable(0,rl.getAAName(AAindex1));
										sysLR.setAllowable(1,rl.getAAName(AAindex2));
										
										sysLR.changeResidueType(m1,0,rl.getAAName(AAindex1),true,true);
										sysLR.applyRotamer(m1, 0, r1);
									
										sysLR.changeResidueType(m1,1,rl.getAAName(AAindex2),true,true);
										sysLR.applyRotamer(m1, 1, r2);
										
										float d = m1.residue[0].getDist(m1.residue[1],true);
										dist[i][j] = Math.min(dist[i][j],d);
										dist[j][i] = Math.min(dist[j][i],d);
										
										if ( (!usePairSt) || (pairE<=pairSt) ) {
											eInteraction[i][j] = Math.max(eInteraction[i][j],Math.abs(pairE));
											eInteraction[j][i] = Math.max(eInteraction[j][i],Math.abs(pairE));
										}
									}
								}							
							}
						}
						
						if (ligPresent) {
							
							int AAindex2 = grl.getAARotamerIndex(ligType);
							
							int numRot2 = grl.getNumRotForAAtype(AAindex2);
							if (numRot2==0)
								numRot2 = 1;
							
							for (int r2=0; r2<numRot2; r2++){
								
								if (!prunedRotAtRes[numInAS*totalNumRotamers + r2]){
									
									float pairE = rs.getMinMatrix()[numInAS][AAindex2][r2][i][AAindex1][r1];
									
									int smallMolResMap[] = new int[2];
									smallMolResMap[0] = molResMap[i];
									smallMolResMap[1] = molResMap[numInAS];
									
									Molecule m1 = getMolRes(m,smallMolResMap);
									
									StrandRotamers sysLR = new StrandRotamers(rl,m1.strand[sysStrNum]);
									sysLR.setAllowable(0,rl.getAAName(AAindex1));									
									sysLR.changeResidueType(m1,0,rl.getAAName(AAindex1),true,true);
									sysLR.applyRotamer(m1, 0, r1);
									
									StrandRotamers ligROT = new StrandRotamers(grl,m1.strand[ligStrNum]);		
									ligROT.setAllowable(0,m1.strand[ligStrNum].residue[0].name);
									ligROT.applyRotamer(m1, 0, r2);
									
									float d = m1.residue[0].getDist(m1.residue[1],true);
									ligDist[i] = Math.min(ligDist[i],d);
									
									if ( (!usePairSt) || (pairE<=pairSt) )
										ligE[i] = Math.max(ligE[i],Math.abs(pairE));
								}
							}
						}
					}
				}
			}
		}
		
		PrintStream logPS = setupOutputFile(runName+".log");
		PrintStream logPS2 = setupOutputFile(runName);
		
		logPS2.println("PIG:0 "+runName); //output in Pigale-compatible ASCII format
		
		//Output data
		for (int i=0; i<numInAS; i++){
			int pdbResNum1 = m.strand[sysStrNum].residue[residueMap[i]].getResNumber();
			for (int j=i+1; j<numInAS; j++){
				int pdbResNum2 = m.strand[sysStrNum].residue[residueMap[j]].getResNumber();
				
				logPS.println(pdbResNum1+" "+pdbResNum2+" "+dist[i][j]+" "+eInteraction[i][j]);
				if ( (dist[i][j]<=distCutoff) && (eInteraction[i][j]>eInteractionCutoff) ) //these two residues interact
					logPS2.println(pdbResNum1+" "+pdbResNum2);
			}
			if (ligPresent){
				int pdbResNum2 = m.strand[ligStrNum].residue[0].getResNumber();
				
				logPS.println(pdbResNum1+" "+pdbResNum2+" "+ligDist[i]+" "+ligE[i]);
				if ( (ligDist[i]<=distCutoff) && (ligE[i]>eInteractionCutoff) )
					logPS2.println(pdbResNum1+" "+pdbResNum2);
			}
		}
		logPS2.println("0 0");
		
		logPS.flush();logPS.close();
		logPS2.flush();logPS2.close();
		
		outputObject(prunedRotAtRes,runName+"_pruneInfo.obj");
	}
	
	//Returns a molecule m1 that contains only the residues in molecule m that are specified by residueMap[] (molecul-relative residue indexing);
	private Molecule getMolRes (Molecule m, int residueMap[]){
		
		Molecule m1 = new Molecule();
		
		for (int i=0; i<m.numberOfStrands; i++){ //create the same number of strands
			m1.addStrand(m.strand[i].name);
			m1.strand[i].isProtein = m.strand[i].isProtein;
		}
		
		for (int i=0; i<residueMap.length; i++){
			
			Residue oldResidue = m.residue[residueMap[i]];
			
			Residue newResidue = new Residue();
			newResidue.name = oldResidue.name;
			newResidue.fullName = oldResidue.fullName;
			
			for (int j=0; j<oldResidue.numberOfAtoms; j++){
				
				Atom oldAtom = oldResidue.atom[j];
				
				Atom newAtom = new Atom(oldAtom.name,oldAtom.coord[0],oldAtom.coord[1],oldAtom.coord[2]);
				newAtom.modelAtomNumber = oldAtom.modelAtomNumber;
				newAtom.strandNumber = oldAtom.strandNumber;
				newAtom.elementType = oldAtom.elementType;
				newResidue.addAtom(newAtom);
			}
			
			m1.addResidue(oldResidue.strandNumber,newResidue);
		}
		
		//Determine the bonds between the atoms in the molecule
		m1.determineBonds();
		
		// Assign the molecule relative atom numbers
		m1.updateMoleculeAtomNumbers();
		
		return m1;
	}
	public void handleCompBranchDGMEC (String s) { //funciton added for BWM - Swati

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutation search parameter filename (string)
		// 3: Branch decomposition filename (string)
		
		// Read System parameters for the reference structure
		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters
		
		int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();		
		String runNameEMatrixMin = (String)(sParams.getValue("MINENERGYMATRIXNAME"));
		boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = (String)(sParams.getValue("LIGTYPE"));
		boolean useEref = (new Boolean((String)sParams.getValue("USEEREF"))).booleanValue();
		boolean prunedRotAtRes[] = (boolean [])readObject(sParams.getValue("PRUNEDROTFILE"),false);
		String bdFile = sParams.getValue("BRANCHDFILE");
		
		Molecule m = new Molecule();
		int numLigRotamers = setupMolSystem(m,sParams,ligPresent,ligType);
		
		int numLevels = numInAS;
		if (ligPresent)
			numLevels++;
		
		int residueMap[] = new int[numInAS];
		String resMapString = (String)sParams.getValue("RESIDUEMAP");
		String resDefault[] = new String[numInAS];
		int molResMap[] = new int[numLevels];
		int invResMap[] = new int[m.numberOfResidues];
		for (int i=0; i<invResMap.length; i++)
			invResMap[i] = -1;
		System.out.print("Mol residue map:");
		for(int i=0;i<numInAS;i++){
			int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
			residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			molResMap[i] = m.mapPDBresNumToMolResNum(pdbResNum);
			invResMap[molResMap[i]] = i;
			System.out.print(" "+molResMap[i]+"("+m.residue[molResMap[i]].fullName+")");
		}
		if (ligPresent) { //ligand is present
			molResMap[numInAS] = m.strand[ligStrNum].residue[0].moleculeResidueNumber;
			invResMap[molResMap[numInAS]] = numInAS;
			System.out.print(" "+molResMap[numInAS]+"("+m.residue[molResMap[numInAS]].fullName+")");
		}
		System.out.println();
		
		RotamerSearch rs = new
RotamerSearch(m,sysStrNum,ligStrNum,hElect,hVDW,hSteric,true,true,0.0f,stericThresh,softStericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,softvdwMultiplier,rl,grl);
		
		System.out.print("Loading precomputed energy matrix...");
		loadPairwiseEnergyMatrices(sParams,rs,runNameEMatrixMin+".dat",false,null);
		System.out.println("done");
		
		long startTimeCPU = CPUTime();
		long startTimeUser = UserTime();
		long startTimeSystem = SystemTime();
		
		//Set the allowable AAs for each AS residue
		boolean addWT = (new Boolean((String)sParams.getValue("ADDWT"))).booleanValue();
		for (int j=0; j<numInAS; j++){
			setAllowablesHelper(rs, sParams, addWT, j, residueMap, resDefault);
		}
		
		float eRef[] = null;
		if (useEref) { //add the reference energies to the min (and max) intra-energies
			eRef = getResEntropyEmatricesEref(useEref,rs.getMinMatrix(),rs.sysLR,residueMap,null,numInAS);
			rs.addEref(eRef, false, ligPresent, numInAS, residueMap);
		}
		
		int numRotForRes[] = compNumRotForRes(numInAS, rs, numLigRotamers, residueMap);
		int numUnprunedRot[] = new int[numRotForRes.length];
		for (int curRes=0; curRes<numInAS; curRes++){			
			int curPruned = 0;
			for (int curRot=0; curRot<totalNumRotamers; curRot++){
				if (prunedRotAtRes[curRes*totalNumRotamers+curRot]) //cur rot is pruned (pruned rotamers are necessarily in the cur set of allowed AAs)
					curPruned++;
			}			
			numUnprunedRot[curRes] = numRotForRes[curRes] - curPruned;
		}
		if (ligPresent) { //ligand is present
			int curPruned = 0;
			for (int curRot=0; curRot<numLigRotamers; curRot++)
				if (prunedRotAtRes[numInAS*totalNumRotamers+curRot])
					curPruned++;
			numUnprunedRot[numInAS] = numLigRotamers - curPruned;
		}
		
//		BranchTree bt = new BranchTree(bdFile,m,numUnprunedRot,molResMap,invResMap,sysStrNum,numInAS,ligPresent);
//		bt.traverseTree(rs.sysLR, rs.ligROT, m, rl, grl, prunedRotAtRes, totalNumRotamers, rotamerIndexOffset, rs.getMinMatrix());
		
		long stopTimeCPU = CPUTime();
		long stopTimeUser = UserTime();
		long stopTimeSystem = SystemTime();
		System.out.println("GBD done");
		System.out.println("Total execution time: CPU "+((stopTimeCPU-startTimeCPU)/(60.0*1000000000.0)));
		System.out.println("Total execution time: User "+((stopTimeUser-startTimeUser)/(60.0*1000000000.0)));
		System.out.println("Total execution time: System "+((stopTimeSystem-startTimeSystem)/(60.0*1000000000.0)));
	}

//////////////////////////////////////////////////////
// Begin Steric Overlap Check Section
//////////////////////////////////////////////////////
	//Compute the amount of overlap between a set of structures and a reference structure
	public void handleCompStericOverlap (String s) {

		// Takes the following parameters
		// 1: System parameter filename for the reference structure (string)
		// 2: System parameter filename for the set of structures to compare (string)
		// 2: Mutation search parameter filename (string)		
		
		// Read System parameters for the reference structure
		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		Molecule mRef = new Molecule();
		setupMolSystem(mRef,sParams,false,null);
		int numInASref = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();		
		int posMapRef[] = new int[numInASref]; //molecule-relative numbering
		String posMapString = (String)sParams.getValue("RESIDUEMAP");
		for(int i=0;i<numInASref;i++)
			posMapRef[i] = (new Integer(getToken(posMapString,i+1))).intValue();
		
		
		
		// Read System parameters for the set of structures to compare
		sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,3)); //read system parameters
		sParams.addParamsFromFile(getToken(s,4)); //read system parameters
		
		String runName = (String)sParams.getValue("RUNNAME");
		
		String protPDBname = (String)sParams.getValue("PROTPDBNAME");
		int numPDBfiles = (new Integer((String)sParams.getValue("NUMPDBFILES"))).intValue();
		int numInAS2 = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
		
		int posMap2[] = new int[numInAS2]; //molecule-relative numbering
		posMapString = (String)sParams.getValue("RESIDUEMAP");
		for(int i=0;i<numInAS2;i++)
			posMap2[i] = (new Integer(getToken(posMapString,i+1))).intValue();
		
		String pdbFiles[] = getPDBfiles(protPDBname,numPDBfiles); //get the PDB filenames
		numPDBfiles = pdbFiles.length;
		
		double minMaxOverlap[] = new double[numInAS2];
		int minMaxOverlapStruct[] = new int[numInAS2];
		for (int i=0; i<minMaxOverlap.length; i++){
			minMaxOverlap[i] = (float)Math.pow(10, 10);
			minMaxOverlapStruct[i] = -1;
		}
		
		// The idea is: given a structure from the set, for each of the residues in posMap2[], determine the largest
		//		steric overlap with an atom in the posMapRef[] residues from the reference structure;
		//		then, for each residue in posMap2[], find the minimum such largest overlap among all structures in the set
		for (int i=0; i<numPDBfiles; i++){
			
			System.out.println("Starting structure "+pdbFiles[i]+" ("+i+")");
			
			sParams.setValue("PDBNAME",pdbFiles[i]);
		
			//Setup the molecule system
			Molecule m2 = new Molecule();
			setupMolSystem(m2,sParams,false,null);
			
			for (int res2=0; res2<numInAS2; res2++){ //for each included residue in the given structure
				double maxOverlap = 0.0;
				for (int at2=0; at2<m2.residue[posMap2[res2]].numberOfAtoms; at2++){ //for each atom in that residue
					Atom a2 = m2.residue[posMap2[res2]].atom[at2];
					if ( hSteric || (!a2.elementType.equalsIgnoreCase("H"))){
						for (int resRef=0; resRef<numInASref; resRef++){ //for each included residue in the reference structure
							for (int atRef=0; atRef<mRef.residue[posMapRef[resRef]].numberOfAtoms; atRef++){ //for each atom
								Atom aRef = mRef.residue[posMapRef[resRef]].atom[atRef];
								if ( hSteric || (!aRef.elementType.equalsIgnoreCase("H"))){
									double overlap = ((a2.radius + aRef.radius)/100.0) - a2.distance(aRef);
									if (overlap<0.0)
										overlap = 0.0;
									maxOverlap = Math.max(maxOverlap, overlap);
								}
							}
						}
					}
				}
				if (minMaxOverlap[res2]>maxOverlap){
					minMaxOverlap[res2] = maxOverlap;
					minMaxOverlapStruct[res2] = i;
				}
			}
		}
		
		
		//Output the computed distances
		PrintStream logPS = setupOutputFile(runName);		
		for (int i=0; i<numInAS2; i++){
			logPS.println(posMap2[i]+" "+minMaxOverlap[i]+" "+pdbFiles[minMaxOverlapStruct[i]]);
		}		
		logPS.flush();
		logPS.close();
	}
	
	//Reads the pdb filenames
	private String [] getPDBfiles(String fName, int numFiles){
		
		BufferedReader bufread = null;
		try {
			File file = new File(fName);
			FileReader fr = new FileReader(file);
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			System.out.println(" ... pdbs config file not found");
			System.exit(1);
		}
		
		String pdbFiles[] = new String[numFiles];

		boolean done = false;
		String str = null;
		int curFile = 0;
		
		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
			 System.out.println("ERROR: An error occurred while reading input");
			 System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else if (getToken(str,1).charAt(0)=='%') //skip comment lines
				continue;
			else {
				if (curFile>=numFiles)
					break;
				else {
					pdbFiles[curFile] = getToken(str,1);
					curFile++;
				}
			}
		}
		
		if (curFile<numFiles){
			String tmp[] = new String[curFile];
			System.arraycopy(pdbFiles, 0, tmp, 0, tmp.length);
			pdbFiles = tmp;
		}
					
		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}
		
		return pdbFiles;
	}
//////////////////////////////////////////////////////
// End Steric Overlap Check Section
//////////////////////////////////////////////////////
	
//////////////////////////////////////////////////////
// Begin Backrub Precomputation Section
//////////////////////////////////////////////////////
	//Compute the amount of overlap between a set of structures and a reference structure
	public void handlePrecomputeBackrubs (String s) {

		// Takes the following parameters
		// 1: System parameter filename
		// 2: Number of backrub samples in each direction
		// 3: Backrub step size
		// 4: Output file name
		
		// Read System parameters for the reference structure
		ParamSet sParams = new ParamSet();
		sParams.addParamsFromFile(getToken(s,2)); //read system parameters
		int numBackrubSamples = new Integer(getToken(s,3)).intValue();
		float backrubStepSize = new Float(getToken(s,4)).floatValue();
		String backrubFile = getToken(s,5);
		
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,false,null);
		int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();		
		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.getValue("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
			residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
		}
		System.out.println();
		
		Amber96ext a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, softvdwMultiplier);
		a96ff.calculateTypesWithTemplates();
		
		BackrubMinimizer brMin = new BackrubMinimizer();
		brMin.initialize(m, a96ff, residueMap, sysStrNum, backrubFile, hSteric, stericThresh);
		brMin.precomputeBackrubs(numBackrubSamples, backrubStepSize);
		
		System.out.println("DONE: Backrub angle precomputation..");
	}
//////////////////////////////////////////////////////
// End Backrub Precomputation Section
//////////////////////////////////////////////////////
/////////////////////////////////////////////////////
//Begin section and function for CPU time computation - added by Swati
///////////////////////////////////////////////////
	public long CPUTime(){
		ThreadMXBean thread = ManagementFactory.getThreadMXBean();
		if(thread.isCurrentThreadCpuTimeSupported())
			return thread.getCurrentThreadCpuTime();
		else
			return 0L;
	}
	public long UserTime(){
		ThreadMXBean thread = ManagementFactory.getThreadMXBean();
		if(thread.isCurrentThreadCpuTimeSupported())
			return thread.getCurrentThreadUserTime();
		else
			return 0L;
	}
	public long SystemTime(){

		ThreadMXBean thread = ManagementFactory.getThreadMXBean();
		if(thread.isCurrentThreadCpuTimeSupported())
			return thread.getCurrentThreadCpuTime() - thread.getCurrentThreadUserTime();
		else
			return 0L;
	}
	
} // end of KSParser class
