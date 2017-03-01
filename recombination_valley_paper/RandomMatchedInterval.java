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
package edu.mit.compbio.qrf.java.utils;



import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.OpenOption;


import static java.nio.file.StandardOpenOption.CREATE;
import static java.nio.file.StandardOpenOption.WRITE;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.ml.distance.EuclideanDistance;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;


import edu.unc.genomics.Interval;

import edu.unc.genomics.io.BigWigFileReader;
import edu.unc.genomics.io.WigFileException;
import edu.usc.epigenome.uecgatk.bissnp.writer.bedObject;




/**
 *
 */
public class RandomMatchedInterval {

	/**
	 * @param args
	 */
	
	@Option(name="-bigWig",usage="input big wig file for matching, default: not enabled")
	public List<String> bigWigFiles = null;
	
	@Option(name="-eps",usage="the minimum distance allowed between input region and random paired region, default: 0.05")
	public double eps = 0.05;
	
	@Option(name="-iteration",usage="the maximum iteration allowed to search the best matched random reigon, default: 1000")
	public int iteration = 1000;
	
	@Option(name="-treatNaAs",usage="treat NA as , default: 0")
	public double treatNaAs = 0.0;

	@Option(name="-randomWholeGenome",usage="randomize the region in the whole genome , default is just random within the same chromosome.. default: not enabled")
	public boolean randomWholeGenome = false;
	
	@Option(name="-regionMode",usage="find the matched value in the specified region rather than the two ends' value , default: not enabled")
	public boolean regionMode = false;
	
	@Option(name="-useMean0",usage="0 means useMean0, 1 or others means not useMean0.  average over bases with non-covered bases counting as zeroes, good for motif density, only work for regionMode,should be matched to the bigWig files' order. default: not enabled")
	public List<Integer> useMean0 = null;

	@Option(name="-includeNA",usage="includeNA, when did not find matched intervals, output NA.  default: not enabled")
	public boolean includeNA = false;
	
	@Option(name="-debugMode",usage="debugMode, default: not enabled")
	public boolean debugMode = false;
	
	final private static String USAGE = "RandomMatchedInterval [opts] hg19.fa.fai input.bed random_matched.bed";

	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	private File genomeIndex = null;
	private PrintWriter bedWriter = null; 
	private List<BigWigFileReader> bigWigFileReaders = null;
	private RandomGenomicRegionGenerator randomGenomicGenerator = null;
	private static final Logger log = Logger.getLogger(RandomMatchedInterval.class);
	
	private static long startTime = -1;
	private static long omitRows = 0;
	private static long lineNum=0;
	
	OpenOption[] options = new OpenOption[] { WRITE, CREATE };
	
	/**
	 * @param args file
	 * @param args bedgraph
	 * @param args matched random bedgraph file
	 */
	public static void main(String[] args)
	throws Exception
	{
			RandomMatchedInterval rmi = new RandomMatchedInterval();
			rmi.doMain(args);
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
			genomeIndex = new File(arguments.get(0));
			String inputBedFile = arguments.get(1);
			String outputBedFile = arguments.get(2);
			
			initiate(outputBedFile);
			log.info("Parsing input bed file ...");
			BufferedReader br = new BufferedReader(new FileReader(inputBedFile));
			String line;
			
			while( (line = br.readLine()) != null){
				if(line.startsWith("#"))
					continue;
				String[] splitin = line.split("\t");
				bedObject randomRegion = matchedRegion(splitin[0],Integer.parseInt(splitin[1]),Integer.parseInt(splitin[2]));
				
				if(randomRegion != null){
					bedWriter.print(randomRegion.getChr() + "\t" + randomRegion.getStart() + "\t" + randomRegion.getEnd() + "\t.\t.");
					for(Object o : randomRegion.getValueObject()){
						bedWriter.print("\t" + o);
					}
					bedWriter.println();
				}else{
					omitRows++;
				}
				lineNum++;
				if(lineNum % 5000 == 0){
					log.info("Processing line: " + lineNum);
				}
			}
			br.close();
			
			finished(arguments.get(2));
			
	}
	
	
	
	private bedObject matchedRegion(String chr, int start, int end) throws IOException, WigFileException{

		Interval randomRegion;
		if(randomWholeGenome){
			randomRegion = randomGenomicGenerator.next(end-start);
		}else{
			randomRegion = randomGenomicGenerator.next(chr, end-start);
		}
				
		//if(bigWigFiles != null && bigWigFiles.length>0){		
		if(bigWigFiles != null && !bigWigFiles.isEmpty()){// when big wig file is provided
			int it = 0;
			if(regionMode){ // measure the value in the whole region
				double[] valuesInInputRegion = new double[bigWigFileReaders.size()];
				int j = 0;
				for(BigWigFileReader bbreader : bigWigFileReaders){				 
					 SummaryStatistics valueInInputRegion = bbreader.queryStats(chr, start, end);
					 if(useMean0.get(j) == 0){
						 valuesInInputRegion[j] = (double)valueInInputRegion.getSum()/(double)(end-start);
						 
					 }else{
						 valuesInInputRegion[j] = valueInInputRegion.getMean();
					 }
					 
					 if(Double.isNaN(valuesInInputRegion[j]))
						 valuesInInputRegion[j]=treatNaAs;
					 j++;
			 	}
				
				while(it < iteration){
					double[] valuesInRandomRegion = new double[bigWigFileReaders.size()];
					j = 0;
					for(BigWigFileReader bbreader : bigWigFileReaders){				 
						 SummaryStatistics valueInRandomRegion = bbreader.queryStats(randomRegion);
						 if(useMean0.get(j) == 0){
							 valuesInRandomRegion[j] = (double)valueInRandomRegion.getSum()/(double)(end-start);
						 }else{
							 valuesInRandomRegion[j] = valueInRandomRegion.getMean();
						 }
						 
						 if(Double.isNaN(valuesInRandomRegion[j]))
							 valuesInRandomRegion[j]=treatNaAs;
						 j++;
				 	}
					EuclideanDistance euclideanDistance = new EuclideanDistance();
					if(valuesInRandomRegion.length > 0 && valuesInInputRegion.length > 0 && euclideanDistance.compute(valuesInInputRegion, valuesInRandomRegion) < eps){
						bedObject bo = UtilsFuns.intervalToBedObject(randomRegion);
						for(double v : valuesInRandomRegion){
							bo.addValue(v);
						}
						bo.addValue((new Interval(chr, start, end)).toString());
						for(double v : valuesInInputRegion){
							bo.addValue(v);
						}
						return bo;
					}
					if(randomWholeGenome){
						randomRegion = randomGenomicGenerator.next(end-start);
					}else{
						randomRegion = randomGenomicGenerator.next(chr, end-start);
					}
					it++;
				}
				
			}else{ // measure the value near start and end seperately
				double[] valuesInInputRegion = new double[bigWigFileReaders.size()*2];
				int j = 0;
				for(BigWigFileReader bbreader : bigWigFileReaders){				 
					 SummaryStatistics valueInInputRegion = bbreader.queryStats(chr, start, start);
					 valuesInInputRegion[j] = valueInInputRegion.getMean();
					 if(Double.isNaN(valuesInInputRegion[j]))
						 valuesInInputRegion[j]=treatNaAs;
					 j++;
					 valueInInputRegion = bbreader.queryStats(chr, end, end);
					 valuesInInputRegion[j] = valueInInputRegion.getMean();
					 if(Double.isNaN(valuesInInputRegion[j]))
						 valuesInInputRegion[j]=treatNaAs;
					 j++;
			 	}
				while(it < iteration){
					double[] valuesInRandomRegion = new double[bigWigFileReaders.size()*2];
					j = 0;
					for(BigWigFileReader bbreader : bigWigFileReaders){				 
						 SummaryStatistics valueInRandomRegion = bbreader.queryStats(randomRegion.getChr(),randomRegion.getStart(),randomRegion.getStart());
						 valuesInRandomRegion[j] = valueInRandomRegion.getMean();
						 if(Double.isNaN(valuesInRandomRegion[j]))
							 valuesInRandomRegion[j]=treatNaAs;
						 j++;
						 valueInRandomRegion = bbreader.queryStats(randomRegion.getChr(),randomRegion.getStop(),randomRegion.getStop());
						 valuesInRandomRegion[j] = valueInRandomRegion.getMean();
						 if(Double.isNaN(valuesInRandomRegion[j]))
							 valuesInRandomRegion[j]=treatNaAs;
						 j++;
				 	}
					EuclideanDistance euclideanDistance = new EuclideanDistance();
					if(valuesInRandomRegion.length > 0 && valuesInInputRegion.length > 0 && euclideanDistance.compute(valuesInInputRegion, valuesInRandomRegion) < eps){
						bedObject bo = UtilsFuns.intervalToBedObject(randomRegion);
						for(double v : valuesInRandomRegion){
							bo.addValue(v);
						}
						bo.addValue((new Interval(chr, start, end)).toString());
						for(double v : valuesInInputRegion){
							bo.addValue(v);
						}
						return bo;
					}
					if(randomWholeGenome){
						randomRegion = randomGenomicGenerator.next(end-start);
					}else{
						randomRegion = randomGenomicGenerator.next(chr, end-start);
					}
					it++;
				}
			}
			if(it >= iteration){
				if(includeNA){
					bedObject bo = UtilsFuns.intervalToBedObject(randomRegion);
					for(BigWigFileReader bbreader : bigWigFileReaders){
						bo.addValue("NA");
					}
					bo.addValue((new Interval(chr, start, end)).toString());
					for(BigWigFileReader bbreader : bigWigFileReaders){
						bo.addValue("NA");
					}
					return bo;
				}else{
					return null;
				}
				
			}
			
			
		}//if big wig is not provided, just return the length matched region. 
			 
		
		bedObject bo = UtilsFuns.intervalToBedObject(randomRegion);
		bo.addValue((new Interval(chr, start, end)).toString());
		return bo;
	}
	

	
	private void initiate(String outputFile) throws IOException{
		startTime = System.currentTimeMillis();
		//System.out.println("HmmHunter started at : " + startTime);
		//System.out.println();
		randomGenomicGenerator  = new RandomGenomicRegionGenerator(genomeIndex);

		bedWriter = new PrintWriter(new File(outputFile));
		
		bigWigFileReaders = new ArrayList<BigWigFileReader>();
		if(bigWigFiles != null && !bigWigFiles.isEmpty()){
			for(int i = 0; i < bigWigFiles.size(); i++){
				//System.err.println(bigWigFiles.get(i) + "\t" + i + "\t" + (new File(bigWigFiles.get(i)).toPath().toString()));
				BigWigFileReader bbreader = new BigWigFileReader((new File(bigWigFiles.get(i))).toPath());
				 //if(!bbreader.isBigWig(p)){
				//	 throw new IOException(bigWigFiles.get(i) + " is not a valid big wig file format!!");
				// }
				bigWigFileReaders.add(bbreader);
			}
		}
		
		
	}
	
	private void finished(String outputFile) throws IOException{
		if(bigWigFiles != null){
			for(int i = 0; i < bigWigFiles.size(); i++){
				bigWigFileReaders.get(i).close();
			}
		}
		
		bedWriter.close();
		
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;


		log.info( omitRows + " rows out of " + lineNum + " total rows ("+ String.format("%.2f", 100*(double)omitRows/(double)lineNum) +"%) can't find matched region");
		log.info("RandomMatchedInterval's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");

	}
	
	//just 
	

}
