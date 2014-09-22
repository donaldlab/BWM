package Miscellany;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;
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
		Set<String> BWMFileNames = new HashSet<String>();
		for(int i = 0; i < files.length; i++)
		{
			if(files[i].getName().contains("BWM"))
			    BWMFileNames.add(files[i].getName());
		}
		/* If it's an AStar conf file, look for the BWM* conf file with the same sequence. */
		for(int i = 0; i < files.length; i++)
		{
		    if(files[i].getName().contains("AStar"))
		    {
		        String AStarFileName = files[i].getName();
		        String BWMFileName = AStarFileName.replace("AStar", "BWM");
		        if(BWMFileNames.contains(BWMFileName))
		            System.out.println(compareKStarScores(BWMFileName.replace("c_BWM_", ""), path+"/"+AStarFileName, path+"/"+BWMFileName));
		    }
		}
		/* compute the scores of both and output to the CSV file. */
	}


	public static String compareKStarScores(String sequence, String AStarFileName, String BWMFileName)
	{
		String output = "Error reading "+sequence+", "+AStarFileName+", "+BWMFileName;
		/* Open the file */
		try
		{
			BufferedReader BWMReader = new BufferedReader(new FileReader(BWMFileName));
			BufferedReader AStarReader = new BufferedReader(new FileReader(AStarFileName));
			
			BigDecimal BWMPartitionFunctionScore = new BigDecimal(0.0);
			BigDecimal AStarPartitionFunctionScore = new BigDecimal(0.0);
			ExpFunction ef = new ExpFunction();
			while(AStarReader.ready())
			{
				String AStarLine = AStarReader.readLine();
				String[] AStarTokens = AStarLine.split(" ");
				String BWMLine = AStarReader.readLine();
				String[] BWMTokens = BWMLine.split(" ");
				/* extract the energy */
				double AStarEnergy = Double.valueOf(AStarTokens[AStarTokens.length -3]);
				BigDecimal AStarConfScore = ef.exp(-AStarEnergy/RotamerSearch.constRT);
				AStarPartitionFunctionScore = AStarPartitionFunctionScore.add(AStarConfScore);
				
				double BWMEnergy = Double.valueOf(BWMTokens[BWMTokens.length -1]);
				BigDecimal BWMConfScore = ef.exp(-BWMEnergy/RotamerSearch.constRT);
				BWMPartitionFunctionScore = AStarPartitionFunctionScore.add(BWMConfScore);
			}
			BigDecimal percentError = AStarPartitionFunctionScore.subtract(BWMPartitionFunctionScore).abs().divide(AStarPartitionFunctionScore, new MathContext(5)).multiply(hundred);
			output = sequence+", "+printBigNum(AStarPartitionFunctionScore,3)
				+", "+printBigNum(BWMPartitionFunctionScore,3)  
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

	
	public static String printBigNum(BigDecimal bd, int sigDigits){
        int numDigits = (bd.toPlainString()).indexOf('.');
        if(numDigits == -1)
                numDigits = (bd.toPlainString().length());

        int newScale = -1*(numDigits - sigDigits);

        return bd.setScale(newScale,BigDecimal.ROUND_HALF_UP).toString();
	}

}
