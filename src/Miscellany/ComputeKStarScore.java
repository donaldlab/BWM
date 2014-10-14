package Miscellany;
import java.io.*;
import java.util.HashSet;
import java.util.ArrayList;
import java.util.Set;
import java.util.StringTokenizer;
import java.io.Serializable;
import java.math.BigDecimal;
import java.math.MathContext;
import kstar.ExpFunction;
import kstar.RotamerSearch;
/***
 * This class takes a string of numbers and computes the partition function.
 * @author jj
 *
 */
public class ComputeKStarScore {
	public static final BigDecimal hundred = new BigDecimal(100);
	

	public static void main (String[] args)
	{
		/* Go through all files */
		String path = args[0];
		File directory = new File(path);
		File[] files = directory .listFiles();
		Set<String> BWMStarFileNames = new HashSet<String>();
		for(int i = 0; i < files.length; i++)
		{
			if(files[i].getName().contains("BWM"))
			    BWMStarFileNames.add(files[i].getName());
		}
		/* If it's an AStar conf file, look for the BWM* conf file with the same sequence. */
		for(int i = 0; i < files.length; i++)
		{
		    if(files[i].getName().contains("AStar"))
		    {
		        String AStarFileName = files[i].getName();
		        String BWMStarFileName = AStarFileName.replace("AStar", "BWM");
		        if(BWMStarFileNames.contains(BWMStarFileName))
				{
		            outputPartitionFunction(BWMStarFileName.replace("c_BWM_", ""), path+"/"+AStarFileName, path+"/"+BWMStarFileName);
				}
		    }
		}
		/* compute the scores of both and output to the CSV file. */
	}


	public static String compareKStarScores(String sequence, String AStarFileName, String BWMStarFileName)
	{
		String output = "Error reading "+sequence+", "+AStarFileName+", "+BWMStarFileName;
		/* Open the file */
		try
		{
			BufferedReader BWMReader = new BufferedReader(new FileReader(BWMStarFileName));
			BufferedReader AStarReader = new BufferedReader(new FileReader(AStarFileName));
			
			BigDecimal BWMStarPartitionFunctionScore = new BigDecimal(0.0);
			BigDecimal AStarPartitionFunctionScore = new BigDecimal(0.0);
			ExpFunction ef = new ExpFunction();
			double minEnergy = 10000;
			while(AStarReader.ready() && BWMReader.ready())
			{
				String AStarLine = AStarReader.readLine();
				String[] AStarTokens = AStarLine.split(" ");
				String BWMLine = BWMReader.readLine();
				String[] BWMStarTokens = BWMLine.split(" ");
				/* extract the energy */
				double AStarEnergy = Double.valueOf(AStarTokens[AStarTokens.length -3]);
				minEnergy = Math.min(AStarEnergy, minEnergy);
				BigDecimal AStarConfScore = ef.exp(-AStarEnergy/RotamerSearch.constRT);
				AStarPartitionFunctionScore = AStarPartitionFunctionScore.add(AStarConfScore);
				
				double BWMEnergy = Double.valueOf(BWMStarTokens[BWMStarTokens.length -1]);
				BigDecimal BWMConfScore = ef.exp(-BWMEnergy/RotamerSearch.constRT);
				BWMStarPartitionFunctionScore = BWMStarPartitionFunctionScore.add(BWMConfScore);
			}
			BigDecimal percentError = AStarPartitionFunctionScore.subtract(BWMStarPartitionFunctionScore).abs().divide(AStarPartitionFunctionScore, new MathContext(5)).multiply(hundred);
			output = minEnergy+", "+printBigNum(AStarPartitionFunctionScore,3)
				+", "+printBigNum(BWMStarPartitionFunctionScore,3)  
				+", "+printBigNum(percentError, 3)+"%";
			BWMReader.close();
			AStarReader.close();
		} 
		catch (IOException e)
		{
			e.printStackTrace();
		}
		return output;
		
	}

	public static String outputPartitionFunction(String sequence, String AStarFileName, String BWMStarFileName)
	{
		String output = "Error reading "+sequence+", "+AStarFileName+", "+BWMStarFileName;
		try
		{
			BufferedReader BWMReader = new BufferedReader(new FileReader(BWMStarFileName));
			BufferedReader AStarReader = new BufferedReader(new FileReader(AStarFileName));
			BufferedWriter AStarWriter  = new BufferedWriter(new OutputStreamWriter(
										new FileOutputStream(sequence+"PartitionFunctionAStar.csv"), "utf-8"));
			BufferedWriter BWMStarWriter  = new BufferedWriter(new OutputStreamWriter(
										new FileOutputStream(sequence+"PartitionFunctionBWM.csv"), "utf-8"));
			BufferedWriter BothWriter  = new BufferedWriter(new OutputStreamWriter(
										new FileOutputStream(sequence+"PartitionFunctionBoth.csv"), "utf-8"));
			
			BigDecimal BWMStarPartitionFunctionScore = new BigDecimal(0.0);
			BigDecimal AStarPartitionFunctionScore = new BigDecimal(0.0);
			ExpFunction ef = new ExpFunction();
			double minEnergy = 10000;
			int numConfs = 0;
			ArrayList<BigDecimal> AStarScores = new ArrayList<BigDecimal>();
			ArrayList<Double> AStarEnergies = new ArrayList<Double>();
			ArrayList<BigDecimal> BWMStarScores = new ArrayList<BigDecimal>();
			ArrayList<Double> BWMStarEnergies = new ArrayList<Double>();
			ArrayList<String> lines = new ArrayList<String>();
			ArrayList<String> AStarConformationList = new ArrayList<String>();
			HashSet<String> AStarConformations = new HashSet<String>();
			ArrayList<String> BWMStarConformationList = new ArrayList<String>();
			HashSet<String> BWMStarConformations = new HashSet<String>();
			while(AStarReader.ready() && BWMReader.ready())
			{
				String AStarLine = AStarReader.readLine();
				String[] AStarTokens = AStarLine.split(" ");
				AStarConformationList.add(ConfFromTokens(AStarTokens));
				AStarConformations.add(ConfFromTokens(AStarTokens));
				String BWMLine = BWMReader.readLine();
				String[] BWMStarTokens = BWMLine.split(" ");
				BWMStarConformationList.add(ConfFromTokens(BWMStarTokens));
				BWMStarConformations.add(ConfFromTokens(BWMStarTokens));

				double AStarEnergy = Double.valueOf(AStarTokens[AStarTokens.length -3]);
				BigDecimal AStarConfScore = ef.exp(-AStarEnergy/RotamerSearch.constRT);
				AStarScores.add(AStarConfScore);
				AStarEnergies.add(AStarEnergy);
				AStarPartitionFunctionScore = AStarPartitionFunctionScore.add(AStarConfScore);
				
				double BWMEnergy = Double.valueOf(BWMStarTokens[BWMStarTokens.length -1]);
				BigDecimal BWMConfScore = ef.exp(-BWMEnergy/RotamerSearch.constRT);
				BWMStarScores.add(BWMConfScore);
				BWMStarEnergies.add(BWMEnergy);
				AStarPartitionFunctionScore = AStarPartitionFunctionScore.add(AStarConfScore);
				BWMStarPartitionFunctionScore = BWMStarPartitionFunctionScore.add(BWMConfScore);
				numConfs++;
			}
			BWMReader.close();
			AStarReader.close();
			double[] BWMStarProbabilities = new double[numConfs];
			double[] AStarProbabilities = new double[numConfs];
			for(int i = 0; i <numConfs; i++)
			{
				String AStarConf = AStarConformationList.get(i);
				double AStarEnergy = AStarEnergies.get(i);
				BigDecimal AStarConfScore = AStarScores.get(i);
				double AStarPercent = AStarConfScore.divide(AStarPartitionFunctionScore, MathContext.DECIMAL32).doubleValue(); 
				System.out.println(sequence+", "+AStarEnergy+","+printBigNum(AStarConfScore,3)+","+printBigNum(AStarPartitionFunctionScore,3)+", "+AStarPercent);
				if(BWMStarConformations.contains(AStarConf))
				{
					BothWriter.write(AStarEnergy+","+printBigNum(AStarConfScore,3)+","+printBigNum(AStarPartitionFunctionScore,3)+", "+AStarPercent+"\n");
				}
				else
				{
					AStarWriter.write(AStarEnergy+","+printBigNum(AStarConfScore,3)+","+printBigNum(AStarPartitionFunctionScore,3)+", "+AStarPercent+"\n");
				}
			}
			for(int i = 0; i <numConfs; i++)
			{
				String BWMStarConf = BWMStarConformationList.get(i);
				double BWMStarEnergy = BWMStarEnergies.get(i);
				BigDecimal BWMStarConfScore = BWMStarScores.get(i);
				double BWMStarPercent = BWMStarConfScore.divide(BWMStarPartitionFunctionScore, MathContext.DECIMAL32).doubleValue(); 
				System.out.println(sequence+", "+BWMStarEnergy+","+printBigNum(BWMStarConfScore,3)+","+printBigNum(BWMStarPartitionFunctionScore,3)+", "+BWMStarPercent);
				if(!AStarConformations.contains(BWMStarConf))
				{
					BWMStarWriter.write(BWMStarEnergy+","+printBigNum(BWMStarConfScore,3)+","+printBigNum(BWMStarPartitionFunctionScore,3)+", "+BWMStarPercent+"\n");
				}
			}
			BothWriter.close();
			AStarWriter.close();
			BWMStarWriter.close();
		} 
		catch (IOException e)
		{
			e.printStackTrace();
		}
		return "";
		
	}

	public static String ConfFromTokens(String[] tokens)
	{
		String output = "";
		for(int i = 1; i < tokens.length-1; i++)
		{
			output+=tokens[i]+" ";
		}
		return output;
	}
	
	public static String printBigNum(BigDecimal bd, int sigDigits){
        int numDigits = (bd.toPlainString()).indexOf('.');
        if(numDigits == -1)
                numDigits = (bd.toPlainString().length());

        int newScale = -1*(numDigits - sigDigits);

        return bd.setScale(newScale,BigDecimal.ROUND_HALF_UP).toString();
	}

}
