/**
 * AlignWig2BedMatrix2CorMatrix.java
 * Dec 16, 2014
 * 10:05:56 AM
 * yaping    lyping1986@gmail.com
 */
package edu.usc.epigenome.dmntools.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.ujmp.core.doublematrix.impl.DefaultDenseDoubleMatrix2D;
import org.ujmp.core.doublematrix.impl.DefaultSparseDoubleMatrix;



/**
 *
 */
public class AlignWig2BedMatrix2CorMatrix {

	/**
	 * @param args
	 */
	
	//@Option(name="-step",usage="matrix step, default:20")
	//public int step = 20;
	@Option(name="-matrix",usage="matrix step, default:2401")
	public int matrix = 2401;
	
	@Option(name="-start",usage="summary start position,for DKO1 paper F2 means -100bp to TSS. default:1100")
	public int start = 1100;
	
	@Option(name="-end",usage="summary end position,for DKO1 paper F2  means +100bp to TSS.  default:1300")
	public int end = 1300;
	
	@Option(name="-extend",usage="summary plot region, so the correlation matrix will be calculated from 1100+500 to 1300+1500, for DKO1 paper F2 means -100(-100 to +400) to +100 (+100 to +500) of TSS. default:500")
	public int extend = 500;
	
	@Option(name="-randomMatch",usage="use random region in the genome for the background correction, default:false")
	public boolean randomMatch = false;

	final private static String USAGE = "AlignWig2BedMatrix2CorMatrix [opts] input.matrix.txt output.correlation.matrix.txt";

	@Argument
	private List<String> arguments = new ArrayList<String>();
	//private File genomeIndex = null;
	
	public static void main(String[] args)
	throws Exception
	{
		AlignWig2BedMatrix2CorMatrix awm = new AlignWig2BedMatrix2CorMatrix();
		awm.doMain(args);
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
				if (arguments.size() < 2) throw new CmdLineException(USAGE);
			
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
			//genomeIndex = new File(arguments.get(0));
			File inputMatrix = new File(arguments.get(0));
			File corMatrix = new File(arguments.get(1));
			PrintWriter writer = new PrintWriter(corMatrix);
			
			BufferedReader br = new BufferedReader(new FileReader(inputMatrix));
			String line;
			ArrayList<Double[]> dataMatrix = new ArrayList<Double[]>();
			long lineNum = 0;
			while( (line = br.readLine()) != null){
				if(line.startsWith("#"))
					continue;
				String[] splitin = line.split("\t");
				Double[] listTmp = new Double[end+extend-start+1];
				for(int i = start+4, j=0; i<= end+extend+4; i++, j++){
					if(splitin[i].equalsIgnoreCase("NA")){
						listTmp[j] = Double.NaN;
					}else{
						listTmp[j] = Double.parseDouble(splitin[i]);
						
					}
					
				}
				dataMatrix.add(listTmp);
				lineNum++;
				if(lineNum % 1000 == 0){
					System.err.println("Processing line: " + lineNum);
				}
			}
			br.close();
		
			DefaultDenseDoubleMatrix2D doubleMatrix = listToMatrix(dataMatrix);
			
			double[][] pearsonScoreMatrix = doubleMatrix.corrcoef(org.ujmp.core.calculation.Calculation.Ret.NEW, true).toDoubleArray();
			//double[][] pearsonScoreMatrix = pearsonScore.computeCorrelationMatrix(doubleMatrix).getData();
			//from start to end, 0 to end-start+1
			for(int i = 0; i < end-start+1; i++){
				for(int j = i+1; j < extend+i; j++){
					writer.print(pearsonScoreMatrix[i][j] + "\t");
				}
				writer.println();
			}
			
			writer.close();
			System.err.println("Finished ! Total lines: " + lineNum);
	}
	
	private DefaultDenseDoubleMatrix2D listToMatrix(ArrayList<Double[]> dataMatrix){
		DefaultDenseDoubleMatrix2D doubleMatrix = new DefaultDenseDoubleMatrix2D(dataMatrix.size(), end+extend-start+1);
		for(int i = 0; i < dataMatrix.size(); i++){
			Double[] t = dataMatrix.get(i);
			for(int j = 0; j < t.length; j++){
				doubleMatrix.setDouble(t[j], i, j);
				
			}
				
			
		}
		
		return doubleMatrix;
	}
	
}
