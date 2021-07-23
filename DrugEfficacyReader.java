import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Iterator; 
import java.io.*;

public class DrugEfficacyReader
{
	public static void main(String[] args)
	{
		final String ComboDrugFile = "ComboDrugGrowth_Nov2017.csv";
		//append variable to string, read 60 files
		try (PrintWriter writer = new PrintWriter(new FileWriter("ComboDrugMat/786-0-output.csv")))
		{
			try (BufferedReader br = new BufferedReader(new FileReader(ComboDrugFile)))
			{
				int linenum = 0;
				br.readLine(); //Skip first line
				linenum++;
				ArrayList<HashMap<String, Integer>> drugmap = new ArrayList<HashMap<String, Integer>>(); //Maps Drug IDs to ComboMat Indices
				HashMap<String, PrintWriter> writers = new HashMap<String, PrintWriter>(); //Maps NCI60 to individual writers
				HashMap<String, Integer> indexmap = new HashMap<String, Integer>(); //Maps NCI60 cell lines to output matrix
				int[] counters = new int[60];
				int matcount = 0;
				int index = 0;
				double[][][] ComboMat = new double[105][105][60]; //Stores ComboIndex values
				String line = br.readLine();
				linenum++;
				String[] data = line.split(",");
				while(line != null)
				{	
					if(data[28].equals("SF-539"))
					{
						break;
					}	
					if(!writers.containsKey(data[28]))
					{
						String filename = null;
						if(data[28].equals("A549/ATCC")) //Backslash in cell line name makes a subdirectory, have to deal with it
						{
							filename = "ComboDrugMat/A549ATCCMat.csv";
						}
						else if(data[28].equals("NCI/ADR-RES"))
						{
							filename = "ComboDrugMat/NCIADR-RESMat.csv";
						}
						else if(data[28].equals("MDA-MB-231/ATCC"))
						{
							filename = "ComboDrugMat/MDA-MB-231ATCCMat.csv";
						}
						else if(data[28].equals("SF-539"))
						{
							filename = "ComboDrugMat/SF-539Mat.csv";
						}
						else
						{
							filename = "ComboDrugMat/" + data[28] + "Mat.csv";
						}
						//System.out.println(filename);
						PrintWriter w = new PrintWriter(new FileWriter(filename));
						writers.put(data[28], w);
					}
					if(!indexmap.containsKey(data[28]))
					{
						indexmap.put(data[28], matcount);
						drugmap.add(new HashMap<String, Integer>());
						matcount++;
					}
					//If the Drug Id hasn't been added to the hasmap and assigned an index, add it
					int curmat = indexmap.get(data[28]); //Get matrix assigned to this cell line
					HashMap<String, Integer> hmap = drugmap.get(curmat); //Get hashmap of drugs to matrix indices for this cell line
					if(!hmap.containsKey(data[8]) && !data[8].equals(""))
					{
						hmap.put(data[8], counters[curmat]);
						counters[curmat]++;
					}
					if(!hmap.containsKey(data[14]) && !data[14].equals(""))
					{
						hmap.put(data[14], counters[curmat]);
						counters[curmat]++;
					}
					double[][] EfficacyMat = new double[5][5];
					if(data[16].equals("0"))//Case 1: Individual Data, skip next 3 lines
					{
						line = br.readLine();
						linenum++;
						if(line == null)
						{
							break;
						}						
						data = line.split(",");
					}
					else //Otherwise Combo Data, get average of ComboScores
					{
						int count = 0;
						String cellline = data[28];
						String drug1 = data[8];
						String drug2 = data[14];
						int index1 = hmap.get(data[8]);
						int index2 = hmap.get(data[14]);
						if(index2 < index1) //Keep data in order
						{
							int swap = index1;
							index1 = index2;
							index2 = swap;
						}
						int test = 0;
						while(data[8].equals(drug1) && data[14].equals(drug2) && data[28].equals(cellline))//While the lines are interactions of the same drug and cell line
						{
							if(!data[25].equals(""))
							{
								int x = Integer.parseInt(data[10]) - 1;
								int y = Integer.parseInt(data[16]) - 1;
								/*if (data[20].isEmpty())
								{
									System.out.println(Arrays.toString(data));
								}*/
								EfficacyMat[x][y] = Double.parseDouble(data[19]);
								count++;
								line = br.readLine();
								linenum++;
								data = line.split(",");
							}
							else
							{
								line = br.readLine();
								linenum++;
								data = line.split(",");
							}
						}
						double Combined_efficacy = 0;
						if (EfficacyMat[4][4] == 0) //Case 1: matrix is 3x3
						{
							Combined_efficacy = (EfficacyMat[1][2] + EfficacyMat[2][1] + EfficacyMat[2][2])/3.0;
						} 
						else //Case 1: matrix is 5x5
						{
							Combined_efficacy = (EfficacyMat[4][3] + EfficacyMat[3][4] + EfficacyMat[4][4])/3.0;
						}			
						ComboMat[index1][index2][curmat] = Combined_efficacy;
					}
				}
				//Write to excel file
				int count = 0;
				Iterator hmIterator = writers.entrySet().iterator(); 
				while (hmIterator.hasNext())
				{ 
			        Map.Entry entry = (Map.Entry)hmIterator.next();
			        PrintWriter w = ((PrintWriter)entry.getValue());
			        int matID = indexmap.get(entry.getKey());
			        for (int i = 0; i < 104; i++)
			        {
					    for (int j = 0; j < 104; j++)
					    {
					        w.write(ComboMat[i][j][matID] + ",");
					    }
					    w.write("\n");
					}
			        w.close();
        		}
				writer.flush();
				br.close();
		    }
		    catch (IOException e)
		    {
		        e.printStackTrace();
		    }
		    writer.close();
		}
		catch (IOException e)
	    {
	        e.printStackTrace();
	    }
	}
}