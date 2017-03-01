/*
 * The MIT License (MIT)
 * Copyright (c) 2015 dnaase <Yaping Liu: lyping1986@gmail.com>

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * RandomMatchedInterval.java
 * Dec 4, 2014
 * 3:54:07 PM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.utils;



import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStreamWriter;
import java.nio.file.OpenOption;

import static java.nio.file.StandardOpenOption.CREATE;
import static java.nio.file.StandardOpenOption.WRITE;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;












import net.sf.javaml.core.kdtree.KDTree;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;








/**
 * only process chromsome by chromosome, so reduce the memory cost
 *  
 */
public class RandomMatchedIntervalByKdtreeByChr {

	/**
	 * @param args
	 */
	

	@Option(name="-percDiff",usage="the percentage of differences allowed between input region and random paired region to search in KD-tree space, default: 10")
	public double percDiff = 10;

	@Option(name="-minDistFromRandom",usage="minimum distance required to find matched random intervals, in order to avoid just find random intervals very closed to existing interval, default: 0")
	public int minDistFromRandom = 0;

	@Option(name="-minDistFracFromRandom",usage="minimum fraction of distance required to find matched random intervals, in order to avoid just find random intervals very closed to existing interval, 0.5 means <=50% overlap with original links. default: 0")
	public double minDistFracFromRandom = 0.0;

	@Option(name="-nearest",usage="when percDiff critera did not find matched case, use the N nearest node in KD-Tree instead, default: 0")
	public int nearest = 0;

	@Option(name="-chrUse",usage="only process chromosome [chrUse]. default: chr22")
	public String chrUse = "chr22";
	
	@Option(name="-excludeRegions",usage="bed files indicated excluded regions. random intervals will be regenerated when there is overlap with the regions. -excludeRegions trackFileName. Default: null")
	public ArrayList<String> excludeRegions = null;

	@Option(name="-excludeCols",usage="indicate which column in the bed file should be excluded. 0-based coordinate. so column 8 will be 7 in this case. allow multiple different column. Default: null")
	public ArrayList<Integer> excludeCols = null;
	
	@Option(name="-bedplus",usage="the first X column is treated as bed item, the rest of the columns are used for KD-tree search values, default: 6")
	public int bedplus = 6;

	
	
	@Option(name="-debugMode",usage="debugMode, default: not enabled")
	public boolean debugMode = false;
	
	final private static String USAGE = "RandomMatchedIntervalByKdtreeByChr [opts] input.bed[.gz] background.bed[.gz] output_random_matched.bed.gz";

	@Argument
	private List<String> arguments = new ArrayList<String>();
	

	private OutputStreamWriter writer = null; 
	private MersenneTwister randomizer = null;
	private static final Logger log = Logger.getLogger(RandomMatchedIntervalByKdtreeByChr.class);
	
	private static long startTime = -1;
	private static long omitRows = 0;
	private static long useNearestRows = 0;
	private static long lineNum=0;
	private static int dimension = -1;
	
	private static double[] meanKey;
	private static double[] meanSquareKey;
	private static long[] countKey;
	private static double[] sdKey;
	
	private static HashSet<Integer> excludeColsSet = null;
	
	
	OpenOption[] options = new OpenOption[] { WRITE, CREATE };
	
	/**
	 * @param args file
	 * @param args bedgraph
	 * @param args matched random bedgraph file
	 */
	public static void main(String[] args)
	throws Exception
	{
			RandomMatchedIntervalByKdtreeByChr rmibk = new RandomMatchedIntervalByKdtreeByChr();
			BasicConfigurator.configure();
			rmibk.doMain(args);
	}

		
	public void doMain(String[] args)
	throws Exception {

			CmdLineParser parser = new CmdLineParser(this);
			// if you have a wider console, you could increase the value;
			// here 80 is also the default
			parser.setUsageWidth(80);
			try
			{
				parser.parseArgument(args);
				if (arguments.size() < 3) throw new CmdLineException(USAGE);
			
			}
			catch (CmdLineException e)
			{
				System.err.println(e.getMessage());
				// print the list of available options
				parser.printUsage(System.err);
				System.err.println();
				return;
			}
			//read input bed file, for each row,
			String inputBedFile = arguments.get(0);
			String backgroundBedFile = arguments.get(1);
			String outputBedFile = arguments.get(2);
			
			initiate(outputBedFile);
			
			if(excludeCols != null && excludeCols.size()>0){
				excludeColsSet = new HashSet<Integer>(excludeCols);
			}
			
			HashMap<String,IntervalTree<Integer>> ignoreLocCollections = null;
			if(excludeRegions!=null && !excludeRegions.isEmpty()){
				log.info("Excluding intervals ... ");
				ignoreLocCollections = new HashMap<String,IntervalTree<Integer>>();
				
				for(String excludeRegion : excludeRegions){
					BufferedReader br = new BufferedReader(new FileReader(excludeRegion));
					String line;
					
					while( (line = br.readLine()) != null){
						if(line.startsWith("#"))
							continue;
						String[] splitLines = line.split("\t");
						if(splitLines.length<3){
							continue;
						}
						String chr = splitLines[0];
						int start = Integer.parseInt(splitLines[1]);
						int end = Integer.parseInt(splitLines[2]);
						IntervalTree<Integer> tree;
						if(ignoreLocCollections.containsKey(chr)){
							tree = ignoreLocCollections.get(chr);
						}else{
							tree = new IntervalTree<Integer>();
						}
						tree.put(start, end, 1);
						ignoreLocCollections.put(chr, tree);
					}
					br.close();
				
				}
			}
			
			log.info("Building KD-Tree for background bed file ...");
			KDTree treeWholeGenome = buildKdtreeWholeGenome(backgroundBedFile, ignoreLocCollections);
 
					
			
			log.info("Parsing input bed file ...");
			BufferedReader br = new BufferedReader(new FileReader(inputBedFile));
			String line;
			
			while( (line = br.readLine()) != null){
				if(line.startsWith("#"))
					continue;
				String[] splitin = line.split("\t");
				
				String chr = splitin[0];
				if(!chr.equalsIgnoreCase(chrUse)){
					continue;
				}
				if(dimension < 0){
					dimension = splitin.length - bedplus + 1;
					if(excludeColsSet != null && excludeColsSet.size()>0){
						dimension -= excludeColsSet.size();
					}
				}
				
				double[] key = new double[dimension];
				double[] lowKey = new double[dimension];
				double[] highKey = new double[dimension];
				key[0] = (double)((int)((Integer.parseInt(splitin[2])-Integer.parseInt(splitin[1]))/1000))*1000;
				key[0]= (key[0]-meanKey[0])/sdKey[0]; 
				lowKey[0] = key[0]*(100-percDiff)/100;
				highKey[0] = key[0]*(100+percDiff)/100;
				//lowKey[0] = key[0];
				//highKey[0] = key[0];
				for(int i = bedplus, index=1; i < splitin.length; i++, index++){
					if(excludeColsSet != null && !excludeColsSet.isEmpty()){
						if(excludeColsSet.contains(i)){
							index--;
							continue;
						}
					}
					if(splitin[i].equalsIgnoreCase("NA")){
						key[index]=Double.NaN;
					}else{
						key[index]=Double.parseDouble(splitin[i]);
						key[index]= (key[index]-meanKey[index])/sdKey[index]; 
					}
					
					lowKey[index]=key[index]*(100-percDiff)/100;
					highKey[index]=key[index]*(100+percDiff)/100;   
				}
				
				
				String radnomMatchedInterval = null;

					
					ArrayList<String> intervalCollection = new ArrayList<String>();
					
						if(nearest > 0){
							try{
								//Object[] os = treeWholeGenome.nearest(key, nearest);
								//System.err.println(os.length);
								//for(Object o : treeWholeGenome.nearest(key, nearest)){
								//	System.err.println(o);
								//}
								for(Object o : treeWholeGenome.nearest(key, nearest)){
									intervalCollection.addAll((ArrayList<String>)o);
								}
								
								useNearestRows++;
							}catch(NullPointerException e){
								log.info("NullPointerException: " + e);
							}
							
							
						}else{
							Object[] target = treeWholeGenome.range(lowKey, highKey);
							if(target == null || target.length == 0){
								
							}else{
								for(Object o : target){
									intervalCollection.addAll((ArrayList<String>)o);
								}
							}
						}
					
					
					if(intervalCollection.size() > 0){
						ArrayList<String> intervalCollectionFinalList = new ArrayList<String>();
						for(int i = 0; i < intervalCollection.size(); i++){
							String intervalString = intervalCollection.get(i);
							String[] intervalStringSeps = intervalString.split("\t");
							if(Math.abs(Integer.parseInt(splitin[2])-Integer.parseInt(intervalStringSeps[2])) > minDistFromRandom &&
									Math.abs(Integer.parseInt(splitin[1])-Integer.parseInt(intervalStringSeps[1])) > minDistFromRandom && 
									((double)Math.abs(Integer.parseInt(splitin[2])-Integer.parseInt(intervalStringSeps[2]))/(double)(Integer.parseInt(splitin[2])-Integer.parseInt(splitin[1]))>minDistFracFromRandom)){
								
								//if(ignoreLocCollections != null && ignoreLocCollections.size() > 0 && ignoreLocCollections.containsKey(intervalStringSeps[0])){
								//	Iterator<Node<Integer>> node = ignoreLocCollections.get(intervalStringSeps[0]).overlappers(Integer.parseInt(intervalStringSeps[1]), Integer.parseInt(intervalStringSeps[2]));
								//	if(node != null && node.hasNext()){
										
								//	}else{
								//		intervalCollectionFinalList.add(intervalString);
								//	}
								//}else{
									intervalCollectionFinalList.add(intervalString);
								//}
								
							}
						}
						if(intervalCollectionFinalList.size()>0){
							radnomMatchedInterval = intervalCollectionFinalList.get(randomizer.nextInt(intervalCollectionFinalList.size()));
						}
						
					}
					
				
				
				if(radnomMatchedInterval != null){
					writer.write(line + "\t" + radnomMatchedInterval + "\n");
					
				}else{
					omitRows++;
				}

				lineNum++;
				if(lineNum % 100000 == 0){
					log.info("Processing line: " + lineNum);
					
				}
				
			}
			br.close();
			
			finished(outputBedFile);
			
	}
	

	//sd: sqrt(N/(N-1))*sqrt(mean(X^2)-(mean(X))^2)
	private KDTree buildKdtreeWholeGenome(String backgroundBedFile, HashMap<String,IntervalTree<Integer>> ignoreLocCollections) throws IOException{
		GZIPInputStream gzipInputStream = null;
		BufferedReader br1;
		if(backgroundBedFile.endsWith(".gz")){
			gzipInputStream = new GZIPInputStream(new FileInputStream(backgroundBedFile));
			br1 = new BufferedReader(new InputStreamReader(gzipInputStream));
			
		}else{
			br1 = new BufferedReader(new FileReader(backgroundBedFile));
		}
		String line;
		while( (line = br1.readLine()) != null){
			if(line.startsWith("#"))
				continue;
			String[] splitin = line.split("\t");
			dimension = splitin.length - bedplus + 1;
			if(excludeColsSet != null && !excludeColsSet.isEmpty()){
				if(excludeColsSet != null && excludeColsSet.size()>0){
					dimension -= excludeColsSet.size();
				}
			}
			break;
		}
		br1.close();
		
		meanKey = new double[dimension];
		meanSquareKey = new double[dimension];
		countKey = new long[dimension];
		for(int i = 0; i < dimension; i++){
			meanKey[i]=0.;
			meanSquareKey[i]=0.;
			countKey[i]=0;
		}
		if(backgroundBedFile.endsWith(".gz")){
			gzipInputStream = new GZIPInputStream(new FileInputStream(backgroundBedFile));
			br1 = new BufferedReader(new InputStreamReader(gzipInputStream));
			
		}else{
			br1 = new BufferedReader(new FileReader(backgroundBedFile));
		}
		while( (line = br1.readLine()) != null){
			if(line.startsWith("#"))
				continue;
			String[] splitin = line.split("\t");
			if(ignoreLocCollections != null && ignoreLocCollections.size() > 0 && ignoreLocCollections.containsKey(splitin[0])){
				Iterator<Node<Integer>> node = ignoreLocCollections.get(splitin[0]).overlappers(Integer.parseInt(splitin[1]), Integer.parseInt(splitin[2]));
				if(node != null && node.hasNext()){
					continue;
				}
			}
			
			
			double len = Double.parseDouble(splitin[2])-Double.parseDouble(splitin[1]);
			countKey[0]++;
			double delta = len - meanKey[0];
			meanKey[0] = meanKey[0] + delta/countKey[0];
			meanSquareKey[0] = meanSquareKey[0] + delta*(len - meanKey[0]);
			
			for(int i = bedplus, index=1; i < splitin.length; i++, index++){
				if(!splitin[i].equalsIgnoreCase("NA")){
					if(excludeColsSet != null && !excludeColsSet.isEmpty()){
						if(excludeColsSet.contains(i)){
							index--;
							continue;
						}
					}
					
					//meanKey[index] += Double.parseDouble(splitin[i]);
					//meanSquareKey[index] += Double.parseDouble(splitin[i])*Double.parseDouble(splitin[i]);
					//countKey[index]++;
					
					countKey[index]++;
					double delta1 = Double.parseDouble(splitin[i]) - meanKey[index];
					meanKey[index] = meanKey[index] + delta1/countKey[index];
					meanSquareKey[index] = meanSquareKey[index] + delta1*(Double.parseDouble(splitin[i]) - meanKey[index]);
					
				}
			}
		}
		br1.close();
		
		//for(int i = 0; i < dimension; i++){
		//	meanKey[i]=meanKey[i]/(double)countKey[i];
		//	meanSquareKey[i]=meanSquareKey[i]/(double)countKey[i];
		//}
		
		sdKey = new double[dimension];
		for(int i = 0; i < dimension; i++){
			//sdKey[i] = Math.sqrt((double)countKey[i]/(double)(countKey[i]-1)) * Math.sqrt(meanSquareKey[i]-meanKey[i]*meanKey[i]);
			sdKey[i] = Math.sqrt(meanSquareKey[i]/(countKey[i]-1));
			log.info("Scaling factor for " + i + "\tmean: " + meanKey[i] + "\tsd:" + sdKey[i] + "\tmeanSquare: " + meanSquareKey[i] + "\tcount:" + countKey[i]);
		}
		
		log.info("Building " + dimension + " dimensional KD-Tree  ...");
		
		KDTree tree = new KDTree(dimension);
		if(backgroundBedFile.endsWith(".gz")){
			gzipInputStream = new GZIPInputStream(new FileInputStream(backgroundBedFile));
			br1 = new BufferedReader(new InputStreamReader(gzipInputStream));
			
		}else{
			br1 = new BufferedReader(new FileReader(backgroundBedFile));
		}
		
		long bgLineNum = 0;
		double[] meanafterKey = new double[dimension];
		double[] meanSquareAfterKey = new double[dimension];
		long[] countAfterKey = new long[dimension];
		double[] sdAfterKey = new double[dimension];
		for(int i = 0; i < dimension; i++){
			meanafterKey[i]=0.;
			meanSquareAfterKey[i]=0.;
			countAfterKey[i]=0;
		}
		
		while( (line = br1.readLine()) != null){
			if(line.startsWith("#"))
				continue;
			String[] splitin = line.split("\t");
			String chr = splitin[0];
			//Interval interval = new Interval(chr, Integer.parseInt(splitin[1]), Integer.parseInt(splitin[2]));
			if(!chr.equalsIgnoreCase(chrUse)){
				continue;
			}
			
			if(ignoreLocCollections != null && ignoreLocCollections.size() > 0 && ignoreLocCollections.containsKey(chr)){
				Iterator<Node<Integer>> node = ignoreLocCollections.get(chr).overlappers(Integer.parseInt(splitin[1]), Integer.parseInt(splitin[2]));
				if(node != null && node.hasNext()){
					continue;
				}
			}
			
			double[] key = new double[dimension];
			key[0] = Double.parseDouble(splitin[2])-Double.parseDouble(splitin[1]);
			key[0] = sdKey[0] == 0 ? (key[0]-meanKey[0]) : (key[0]-meanKey[0])/sdKey[0];
			//meanafterKey[0] += key[0];
			//meanSquareAfterKey[0] += key[0]*key[0];
			//countAfterKey[0]++;
			
			countAfterKey[0]++;
			double delta = key[0] - meanafterKey[0];
			meanafterKey[0] = meanafterKey[0] + delta/countAfterKey[0];
			meanSquareAfterKey[0] = meanSquareAfterKey[0] + delta*(key[0] - meanafterKey[0]);
			
			
			for(int i = bedplus, index=1; i < splitin.length; i++, index++){
				if(excludeColsSet != null && !excludeColsSet.isEmpty()){
					if(excludeColsSet.contains(i)){
						index--;
						continue;
					}
				}
				if(splitin[i].equalsIgnoreCase("NA")){
					key[index]=Double.NaN;
				}else{
					key[index]=Double.parseDouble(splitin[i]);
					key[index]= sdKey[index] == 0 ? (key[index]-meanKey[index]) : (key[index]-meanKey[index])/sdKey[index];  //standardize the result..
					//meanafterKey[index] += key[index];
					//meanSquareAfterKey[index] += key[index]*key[index];
					//countAfterKey[index]++;
					
					countAfterKey[index]++;
					double delta1 = key[index] - meanafterKey[index];
					meanafterKey[index] = meanafterKey[index] + delta1/countAfterKey[index];
					meanSquareAfterKey[index] = meanSquareAfterKey[index] + delta1*(key[index] - meanafterKey[index]);
				}
			}
			
			if(tree.search(key) != null){
				ArrayList<String> intervals = (ArrayList<String>) tree.search(key);
				intervals.add(line);
				tree.insert(key, intervals);
			}else{
				ArrayList<String> intervals = new ArrayList<String>();
				intervals.add(line);
				tree.insert(key, intervals);
			}
			
			bgLineNum++;
			if(bgLineNum % 1000000 == 0){
				log.info("Processing backgrond file line: " + bgLineNum);
				
			}
			
		}
		
		br1.close();
		
		//for(int i = 0; i < dimension; i++){
		//	meanafterKey[i]=meanafterKey[i]/(double)countAfterKey[i];
		//	meanSquareAfterKey[i]=meanSquareAfterKey[i]/(double)countAfterKey[i];
		//}

		for(int i = 0; i < dimension; i++){
			//sdAfterKey[i] = Math.sqrt((double)countAfterKey[i]/(double)(countAfterKey[i]-1)) * Math.sqrt(meanSquareAfterKey[i]-meanafterKey[i]*meanafterKey[i]);
			sdAfterKey[i] = Math.sqrt(meanSquareAfterKey[i]/(countAfterKey[i]-1));
			log.info("After standization, scaling factor for " + i + "\tmean: " + meanafterKey[i] + "\tsd:" + sdAfterKey[i] + "\tmeanSquare: " + meanSquareAfterKey[i] + "\tcount:" + countAfterKey[i]);
		}
		//System.exit(1);
		return tree;
	}
	


	
	private void initiate(String outputFile) throws IOException{
		startTime = System.currentTimeMillis();

		randomizer = new MersenneTwister();
		writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile)), "UTF-8");
		
	}
	
	private void finished(String outputFile) throws IOException{
		
		writer.close();
		
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;


		log.info( omitRows + " rows out of " + lineNum + " total rows ("+ String.format("%.2f", 100*(double)omitRows/(double)lineNum) +"%) can't find matched region");
		log.info( useNearestRows + " rows could only be found by nearest point");
		log.info("RandomMatchedInterval's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");

	}
	
	//just 
	

}
