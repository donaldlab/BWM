package Tests;

import kstar.PairwiseEnergyMatrix;
import java.io.*;

public class EnergyMatrixComparison
{
	public static void main(String[] arg)
	{
		System.out.println("Start.");
		String matrixFile1 = arg[0];
		String matrixFile2 = arg[1];
		PairwiseEnergyMatrix eMatrix1 = new PairwiseEnergyMatrix();
		eMatrix1.eMatrix = (float[][][][][][]) readObject(matrixFile1);
		PairwiseEnergyMatrix eMatrix2 = new PairwiseEnergyMatrix();
		eMatrix2.eMatrix = (float[][][][][][]) readObject(matrixFile2);
		compare(eMatrix1.eMatrix, eMatrix2.eMatrix);
		System.out.println("Finished.");

	}

	static void compare(float[][][][][][] m1, float[][][][][][] m2)
	{
		if(m1.length != m2.length)
			return;
		for(int i = 0; i < m1.length; i++)
		{
			if(m1[i] == null || m2[i] == null) 
			{
				if( m1[i] != m2[i])
				{
					System.out.println("NULL ENTRY MISMATCH!!");
					return;
				}
				continue;
			}
			if(m1[i].length != m2[i].length)
				return;
			for(int j = 0; j < m1[i].length; j++)
			{
				if(m1[i][j] == null || m2[i][j] == null) 
				{
					if( m1[i][j] != m2[i][j])
					{
						System.out.println("NULL ENTRY MISMATCH!!");
						return;
					}
					continue;
				}
				if(m1[i][j].length != m2[i][j].length)
					return;
				for(int k = 0; k < m1[i][j].length; k++)
				{
					if(m1[i][j][k] == null || m2[i][j][k] == null) 
					{
						if( m1[i][j][k] != m2[i][j][k])
						{
							System.out.println("NULL ENTRY MISMATCH!!");
							return;
						}
						continue;
					}
					if(m1[i][j][k].length != m2[i][j][k].length)
						return;
					for(int l = 0; l < m1[i][j][k].length; l++)
					{
						if(m1[i][j][k][l] == null || m2[i][j][k][l] == null) 
						{
							if( m1[i][j][k][l] != m2[i][j][k][l])
							{
								System.out.println("NULL ENTRY MISMATCH!!");
								return;
							}
							continue;
						}
						if(m1[i][j][k][l].length != m2[i][j][k][l].length)
							return;
						for(int m = 0; m < m1[i][j][k][l].length; m++)
						{
							if(m1[i][j][k][l][m] == null || m2[i][j][k][l][m] == null) 
							{
								if(m1[i][j][k][l][m] != m2[i][j][k][l][m])
								{
									System.out.println("NULL ENTRY MISMATCH!!");
									return;
								}
								continue;
							}
							if(m1[i][j][k][l][m].length != m2[i][j][k][l][m].length)
								return;
							for(int n = 0; n < m1[i][j][k][l][m].length; n++)
							{
								if(m1[i][j][k][l][m][n] != m2[i][j][k][l][m][n])
								{
									if(m1[i][j][k][l][m][n] > 100 
											|| m1[i][j][k][l][m][n] 
											/ Math.abs(m1[i][j][k][l][m][n] - m2[i][j][k][l][m][n]) > 100
											|| Math.abs(m1[i][j][k][l][m][n] - m2[i][j][k][l][m][n]) < 0.001)
										continue;
									System.out.println("["+i+"]"+"["+j+"]"+"["+k+"]"+"["+l+"]"+"["+m+"]"+"["+n+"]");
									System.out.println("m1:" + m1[i][j][k][l][m][n]+", m2:"+m2[i][j][k][l][m][n]
											+", difference "+(m1[i][j][k][l][m][n] - m2[i][j][k][l][m][n]));
								}
							}
						}
					}
				}
			}
		}
	}

	static Object readObject(String inFile){
		return readObject(inFile,false);
	}

	static Object readObject(String inFile, boolean repeat){
		Object inObj = null;
		boolean done = false;
		while (!done){
			try{
				ObjectInputStream in = new ObjectInputStream(new FileInputStream(inFile));
				inObj = in.readObject();
				in.close();
				done = true;
			}
			catch (ClassNotFoundException e){
				System.out.println(e);
				e.printStackTrace();
				System.out.print(e.getCause());
				System.out.println("Your object is lying!");
				System.out.println("ERROR ERROR ERROR ERROR ERROR");
				System.out.println("ERROR ERROR ERROR ERROR ERROR");
				System.out.println("ERROR ERROR ERROR ERROR ERROR");
				if (repeat)
					done = false;
				else
					done = true;
				System.exit(-1);
			}
			catch (Exception e)
			{
				System.out.println(e);
				e.printStackTrace();
				System.out.println("ERROR: An exception occurred while reading from object file");
				if (repeat)
					done = false;
				else
					done = true;
			}

		}
		return inObj;
	}



}
