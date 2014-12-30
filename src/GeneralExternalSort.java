// adapted from http://www.codeodor.com/index.cfm/2007/5/14/Re-Sorting-really-BIG-files---the-Java-source-code/1208
//  Re: Sorting really BIG files - the Java source code
// Posted by Sam on May 14, 2007 at 08:01 AM UTC - 6 hrs

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GeneralExternalSort {
public static  boolean reverse;

public static int sort_col,name_col;
	public static void main(String[] args) {
		if(args.length == 0) { 
			System.out.println("USAGE: <name_col> <sort_col> <decreasing(0,1)> <num_lines_per_file> <filename>");
			return;
		}
		name_col = Integer.parseInt(args[0]);
		sort_col = Integer.parseInt(args[1]);
		if(args[2].equals("1"))
			reverse = true;
		else
			reverse = false;
		
		int numfiles = Integer.parseInt(args[3]);
		String file = args[4];
		
		externalSort(file,numfiles);
	}
	
	public static void externalSort(String inputfile,int numlines)
	{
		System.out.println("reading file in chunks of " + numlines + ", sorting, and writing them to temp files.");
	     try
	     {
	         FileReader inputFileReader = new FileReader(inputfile);
	         BufferedReader in = new BufferedReader(inputFileReader);
	         ArrayList<GeneralSortElement> rows = new ArrayList<GeneralSortElement>();
	                    String[] row;
	         int numFiles = 0;
	         String line = "";
	         
	         /******** 
	          * SET VARS HERE
	          */
	         System.out.println("SORTING ON COLUMN "+sort_col+ " (decreasing order? " + reverse+")");
	         
	         while (line!=null)
	         {
	        	 System.out.println("  file #" + numFiles);
	             // get 1million rows
	        	 for(int i=0; i<numlines; i++) {
	        		 line = in.readLine();
	        		 if (line==null)
	        			 break;
	        		 row = line.split("\\s+");
	        		 try {
	        			 double sortval = Double.parseDouble(row[sort_col]);
	        			 rows.add(new GeneralSortElement(line,row[name_col],sortval));
	        		 } catch (NumberFormatException e) {
	        			 System.out.println("Skipping " + line);
	        		 }
	        	 }
	        	 // sort the rows
	        	 Collections.sort(rows);

	             // write to disk
	             FileWriter fw = new FileWriter(inputfile + "_chunk" + numFiles + ".txt");
	             BufferedWriter bw = new BufferedWriter(fw);
	             for(int i=0; i<rows.size(); i++) {
	            	 bw.append(rows.get(i).read+"\t"+rows.get(i).perc+"\t"+rows.get(i).row+"\n");
	             }
	             bw.close();
	             fw.close();
	             numFiles++;
	             rows.clear();
	         }
	         in.close();
	         
	         System.out.println("there are " + numFiles + " files.");
	         mergeFiles(inputfile, numFiles);
	        
	     }
	     catch (Exception ex)
	     {
	         ex.printStackTrace();
	         System.exit(-1);
	     }
	    
	    
	}

	private static void mergeFiles(String inputfile, int numFiles)
	{
		System.out.println("Merging files...");
	     try
	     {
	         ArrayList<FileReader> mergefr = new ArrayList<FileReader>();
	         ArrayList<BufferedReader> mergefbr = new ArrayList<BufferedReader>();
	         ArrayList<GeneralSortElement> filerows = new ArrayList<GeneralSortElement>();
	         FileWriter outFileWriter = new FileWriter(inputfile + ".sorted");
	         BufferedWriter out = new BufferedWriter(outFileWriter);
	            
	       String[] splitline;
	         boolean someFileStillHasRows = false;
	        String tmp;
	         for (int i=0; i<numFiles; i++)
	         {
	             mergefr.add(new FileReader(inputfile+"_chunk"+i+".txt"));
	             mergefbr.add(new BufferedReader(mergefr.get(i)));
	           
	             // get the first row
	             String line = mergefbr.get(i).readLine();
	             if (line != null) {
	            	 splitline = line.split("\t");
	            	 tmp = splitline[2];
	            	 for(int j=3;j<splitline.length;j++)
	            		 tmp+="\t"+splitline[j];

	            	 filerows.add(new GeneralSortElement(tmp,splitline[0],Double.parseDouble(splitline[1])));
	                 someFileStillHasRows = true;
	             } else {
	                 filerows.add(null);
	             }
	                
	         }
	        
	         GeneralSortElement row;
	         int cnt = 0;
	         while (someFileStillHasRows)
	         {
	             GeneralSortElement min;
	             int minIndex = 0;
	            
	             row = filerows.get(0);
	             if (row!=null) {
	                 min = row;
	                 minIndex = 0;
	             }
	             else {
	                 min = null;
	                 minIndex = -1;
	             }
	            
	             // check which one is min
	             for(int i=1; i<filerows.size(); i++)
	             {
	                 row = filerows.get(i);
	                 if (min!=null) {
	                    
	                     if(row!=null && row.compareTo(min) < 0)
	                     {
	                         minIndex = i;
	                         min = row;
	                     }
	                 }
	                 else
	                 {
	                     if(row!=null)
	                     {
	                         min = row;
	                         minIndex = i;
	                     }
	                 }
	             }
	            
	             if (minIndex < 0) {
	                 someFileStillHasRows=false;
	             }  else {
	                 // write to the sorted file
	                 out.append(filerows.get(minIndex).row+"\n");
	                
	                 // get another row from the file that had the min
	                 String line = mergefbr.get(minIndex).readLine();
	                 if (line != null)
	                 {
	                	 splitline = line.split("\t");
		            	 tmp = splitline[2];
		            	 for(int j=3;j<splitline.length;j++)
		            		 tmp+="\t"+splitline[j];

		            	 filerows.set(minIndex,new GeneralSortElement(tmp,splitline[0],Double.parseDouble(splitline[1])));
	                 }
	                 else
	                 {
	                     filerows.set(minIndex,null);
	                 }
	             }                                
	             // check if one still has rows
	             for(int i=0; i<filerows.size(); i++)
	             {
	                
	                 someFileStillHasRows = false;
	                 if(filerows.get(i)!=null)
	                 {
	                     if (minIndex < 0)
	                     {
	                         System.out.println("mindex lt 0 and found row not null" + filerows.get(i).row);
	                         System.exit(-1);
	                     }
	                     someFileStillHasRows = true;
	                     break;
	                 }
	             }
	            
	             // check the actual files one more time
	             if (!someFileStillHasRows)
	             {
	                
	                 //write the last one not covered above
	                 for(int i=0; i<filerows.size(); i++)
	                 {
	                     if (filerows.get(i) == null)
	                     {
	                         String line = mergefbr.get(i).readLine();
	                         if (line!=null)
	                         {
	                        	 someFileStillHasRows=true;
	                             
	                             splitline = line.split("\t");
				            	 tmp = splitline[2];
				            	 for(int j=3;j<splitline.length;j++)
				            		 tmp+="\t"+splitline[j];
				            	 
				            	 filerows.set(i,new GeneralSortElement(tmp,splitline[0],Integer.parseInt(splitline[1])));
	                         }
	                     }
	                            
	                 }
	             }
	                
	         }
	        
	        
	        
	         // close all the files
	         out.close();
	         outFileWriter.close();
	         for(int i=0; i<mergefbr.size(); i++)
	             mergefbr.get(i).close();
	         for(int i=0; i<mergefr.size(); i++)
	             mergefr.get(i).close();

	     }
	     catch (Exception ex)
	     {
	         ex.printStackTrace();
	         System.exit(-1);
	     }
	}
	
}
