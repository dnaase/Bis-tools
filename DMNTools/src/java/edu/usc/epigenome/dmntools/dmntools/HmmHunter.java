/**
 * 
 */
package edu.usc.epigenome.dmntools.dmntools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistributionImpl;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.dmntools.distribution.OpdfBeta;
import edu.usc.epigenome.dmntools.distribution.OpdfBetaFactory;
import edu.usc.epigenome.dmntools.hmm.BbForwardBackwardScaledCalculator;
import edu.usc.epigenome.dmntools.hmm.BbViterbiCalculator;
import edu.usc.epigenome.dmntools.hmm.NdrForwardBackwardScaledCalculator;
import edu.usc.epigenome.dmntools.hmm.ObservationMethy;
import edu.usc.epigenome.dmntools.utils.BisulfiteGenomicLocHmm;
import edu.usc.epigenome.dmntools.utils.FisherExactTest;
import edu.usc.epigenome.dmntools.utils.GenomeLocus;
import edu.usc.epigenome.uecgatk.bissnp.writer.bedObject;
import edu.usc.epigenome.uecgatk.bissnp.writer.bedObjectWriterImp;


import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardScaledCalculator;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.learn.KMeansLearner;


/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Oct 16, 2013 6:26:03 PM
 * 
 */
public abstract class HmmHunter {
	
	/**
	 * @param args
	 */
	@Option(name="-minCT",usage="minimum CT read coverage required to be count in methylation level, default: 1")
	public int minCT = 1;

	//options for HMM model
	@Option(name="-onlyTrain",usage="only enable the training mode and output HMM parameter, default: not enabled")
	public boolean onlyTrain = false;
	
	@Option(name="-onlyDecode",usage="only enable the decoding step and output segments, default: not enabled")
	public boolean onlyDecode = false;

	@Option(name="-states",usage="number of states in HMM model, default: 2")
	public int states = 2;
	
	@Option(name="-hmmFile",usage="read/write HMM model from/to a file, default: read HMM parameters from a file", required=true)
	public String hmmFile = null;
	
	@Option(name="-tol",usage="tolerence level for the converge, default: 1e-5")
	public double tol = 1e-5;
	
	@Option(name="-iteration",usage="maximum number of iteration for HMM model converge, default: 20")
	public int iteration = 20;

	@Option(name="-beta",usage="using beta model rather than beta binomial model in the decoding step and output segments, default: not enabled")
	public boolean beta = false;

	@Option(name="-kmeans",usage="use k means method rather than randomly to initiate the best initial parameters, when number of GCH/HCG is very large, it is not efficient to do it. default: not enabled")
	public boolean kmeans = false;
	
	@Option(name="-bedFormat",usage="input bed format. 1: 6plus2 bed format. 2: standard bed format, default: 1")
	public int bedFormat = 1;
	
	//options to calculate p value for each segment	
	@Option(name="-sigTestMode",usage="The mode to calculate segment p value [1:permutation, 2:binomial, 3:fisher, 4:betaDiff]. default: 2")
	public int sigTest = 2;
	
	@Option(name="-gap",usage="max gap size between two GCH/HCG. If two GCH distrance are more than -gap, it will restart a new segment, rather than estimate a transition probability, default: 10000")
	public int gap = 10000;
	
	@Option(name="-dataP",usage="minimum number of data points in each segment, default: 5")
	public int dataP = 5;

	//options for random permutation
	@Option(name="-permutNum",usage="number of random permutation to calculate p value, default: 100")
	public int permutNum = 100;
	
	//option for halfLocal decoding (simple binomial test)
	@Option(name="-adjWindow",usage="in the binomial mode, define the adjacent window size to compare with (100 means +/-100 bp around segment), default: 1000000")
	public int window = 1000000;
	
	@Option(name="-ctInAdjWindow",usage="in the binomial mode, minimum of CT reads in the adjacent window, default: 100")
	public int ctInWindow = 100;
	
	@Option(name="-cInAdjWindow",usage="in the binomial mode, minimum of Cytosine position in the adjacent window, default: 10")
	public int cInAdjWindow = 10;
	
	@Option(name = "-gchInAdjWindow", usage = "in the binomial mode, minimum number of GCH in the window, default: 10")
	public int gchInWindow = 10;
	
	private final static int MINIMUM_DATA_POINTS = 2;
	final private static int PURGE_INTERVAL = 1000000; // interval to show the progress.
	
	private static long startTime = -1;
	
	public ArrayList<ArrayList<ObservationMethy>> methyObj;
	public ArrayList<ArrayList<ObservationReal>> methyValue;
	public ArrayList<ArrayList<GenomeLocus>> position;
	
	
	protected bedObjectWriterImp stateWriter = null;
	protected bedObjectWriterImp segWriter = null;  //mar and mpr are all output into the same file

	protected FileWriter hmmWriter = null;
	
	protected void initiate(String prefix) throws IOException{
		startTime = System.currentTimeMillis();
		//System.out.println("HmmHunter started at : " + startTime);
		//System.out.println();
		
		String hmmStateSeqFile = prefix + ".HMMstates.txt";
		String marFile = prefix + ".segment.txt";

		
		if(!onlyTrain){
			stateWriter = new bedObjectWriterImp(new File(hmmStateSeqFile));
			stateWriter.setBedgraph();
			stateWriter.addHeader("#chr\tstart\tend\tstate\tnum_C\tnum_CT\n");
			segWriter = new bedObjectWriterImp(new File(marFile));
			segWriter.setBedgraph();
			String header = "#chr\tstart\tend\tstate_of_segment\tave_methy(%)\tnum_GCH\ttotal_num_methy_C_in_segment\ttotal_num_CT_in_segment\ttotal_num_methy_C_in_adj_segment\ttotal_num_CT_in_adj_segment\tstart_bound_dist_to_pre_GCH\t" +
					"start_bound_dist_to_next_GCH\tend_bound_dist_to_pre_GCH\tend_bound_dist_to_next_GCH\ttransition_probability_to_segment\ttransition_probability_out_segment\tp_value\n";
			segWriter.addHeader(header);

		}
		if(!onlyDecode){
			hmmWriter = new FileWriter(hmmFile);
		}
		
		methyObj = new ArrayList<ArrayList<ObservationMethy>>();
		methyValue = new ArrayList<ArrayList<ObservationReal>>();
		position = new ArrayList<ArrayList<GenomeLocus>>();
	}
	
	
	
	protected void finished() throws IOException{
		if(!onlyTrain){
			stateWriter.close();
			segWriter.close();
		}
		long endTime   = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		//System.out.println();
		//System.out.println("HmmHunter finshed at : " + endTime);
		System.out.println();
		System.out.println("HmmHunter's running time is: " + TimeUnit.DAYS.toHours(totalTime));
		
	}
	
	/*
	 * @param:bedFile provide input bed file
	 * it parse 6plus2 bed file, and push all position and methylation value into the array list 
	 */
	protected void parseBedFile(String bedFile) throws Exception,FileNotFoundException, IOException{
		System.out.println("Parsing input bed file ...");
		BufferedReader br = new BufferedReader(new FileReader(bedFile));
		String line;
		ArrayList<ObservationMethy> methyObjIsland = new ArrayList<ObservationMethy>();
		ArrayList<ObservationReal> methyValueIsland = new ArrayList<ObservationReal>();
		ArrayList<GenomeLocus> positionIsland = new ArrayList<GenomeLocus>();

		while( (line = br.readLine()) != null){
			if(!line.startsWith("chr"))
				continue;
			String[] splitin = line.split("\t");
			double tmp = Double.NaN;
			int coverage = Integer.MIN_VALUE;
			if(bedFormat == 1){
				if(splitin.length != 8)
					throw new Exception("Not 6plus2 bed format! 6plus2 bed file format should be like: \n chr\tstart\tend\tname\tscore\tstrand\tmethylation(%)\tnum_CT_coverage\nchr1\t1001\t1002\t.\t750\t+\t75.0\t10\nchr1\t1005\t1006\t.\t200\t-\t20.0\t2\n");
				tmp = Double.parseDouble(splitin[6])/100.0;
				coverage = Integer.parseInt(splitin[7]);
				
			}else if(bedFormat == 2){
				if(splitin.length < 6)
					throw new Exception("Not standard bed format! 6plus2 bed file format should be like: \n chr\tstart\tend\tname\tscore\tstrand\n");
				tmp = Double.parseDouble(splitin[3])/100.0;
				coverage = Integer.parseInt(splitin[4]);
			}else{
				throw new Exception("Not such a bed format allowed\n");
			}
			if(tmp == 1.0){
				tmp-=(Math.random())/1000000000;
			}
			if(tmp == 0.0){
				tmp+=(Math.random())/1000000000;
			}
			
			GenomeLocus loc = new GenomeLocus(splitin[0], Integer.parseInt(splitin[1]), Integer.parseInt(splitin[2]));
			ObservationReal methyTmp = new ObservationReal(tmp);
			ObservationMethy methyObjTmp = new ObservationMethy(tmp);
			methyObjTmp.setCoverage(coverage);
			int dist = 0;
			if(!positionIsland.isEmpty()){
				dist =positionIsland.get((positionIsland.size()-1)).distance(loc);
			}
			 
			methyObjTmp.setDistance(dist);  //TODO: need to check whether this should be set 0 in the beginning of each island/chromsome?
			if(!positionIsland.isEmpty() && positionIsland.get((positionIsland.size()-1)).onSameContig(loc) && dist <= gap){
				positionIsland.add(loc);
				methyValueIsland.add(methyTmp);
				methyObjIsland.add(methyObjTmp);
				
			}else{
				if(positionIsland.size() >= MINIMUM_DATA_POINTS){
					methyObj.add( new ArrayList<ObservationMethy>(methyObjIsland));
					methyValue.add( new ArrayList<ObservationReal>(methyValueIsland));
					position.add( new ArrayList<GenomeLocus>(positionIsland));
				}
				methyObjIsland.clear();
				methyValueIsland.clear();
				positionIsland.clear();
				
				positionIsland.add(loc);
				methyValueIsland.add(methyTmp);
				methyObjIsland.add(methyObjTmp);
			}
			//System.err.println(methyValue.size() + "\t" + positionIsland.size());
			
		}
		methyObjIsland.clear();
		methyValueIsland.clear();
		positionIsland.clear();
		System.out.println("Training island total number: " + position.size());
		System.out.println("Training island size\tChr\tDistance to the last island");
		int preEnd = 0;
		String preChr = null;
		for( ArrayList<GenomeLocus> it : position){
			if(preChr == null){
				preChr = it.get(it.size()-1).getChr();
				preEnd = it.get(it.size()-1).getStart();
				System.out.println(it.size() + "\t" + preChr + "\t" + 0);
			}else if(!it.get(it.size()-1).getChr().equalsIgnoreCase(preChr)){
				preChr = it.get(it.size()-1).getChr();
				preEnd = it.get(it.size()-1).getStart();
				System.out.println(it.size() + "\t" + preChr + "\t" + 0);
			}else{
				System.out.println(it.size() + "\t" + preChr + "\t" + (it.get(0).getStart()-preEnd));
				preEnd = it.get(it.size()-1).getStart();
			}
			
		}
		br.close();
	}

	
	protected Hmm<ObservationReal> buildInitHmmRandomlyByBeta(ArrayList<ArrayList<ObservationReal>> seqs)
	{	
		Hmm<ObservationReal> hmm = 
				new Hmm<ObservationReal>(2,new OpdfBetaFactory());
		double tmp3 = Math.random();
			hmm.setPi(0, tmp3);
			hmm.setPi(1, 1-tmp3);
			
			hmm.setOpdf(0, new OpdfBeta(0.03,1.1));
			hmm.setOpdf(1, new OpdfBeta(1.0,0.3));
			double tmp1 = Math.random();
			
			hmm.setAij(0, 1, tmp1);
			hmm.setAij(0, 0, 1-tmp1);
			double tmp2 = Math.random();
			hmm.setAij(1, 0, tmp2);
			hmm.setAij(1, 1, 1-tmp2);
		return hmm;

	}
	
	protected Hmm<ObservationReal> buildInitHmmByBeta(ArrayList<ArrayList<ObservationReal>> seqs)
	{	

		KMeansLearner<ObservationReal> kl = new KMeansLearner<ObservationReal>(states, new OpdfBetaFactory(),
				seqs);
		System.out.println("KMeansLearner...");
		Hmm<ObservationReal> hmm = kl.learn();
		return hmm;
	}
	
	abstract protected void trainHmm() throws IOException;
	
	abstract protected void decodeHmm() throws Exception;
	
	protected void segmentHmmState(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, NdrForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm) throws Exception{

		
		if(sigTest == 1){
			segmentHmmStateByRandomPermutation(loci, methyState, hiddenState, nfbsc,hmm);
		}else if(sigTest == 2){
			segmentHmmStateByBinomialTestWithFBSC(loci, methyState, hiddenState, nfbsc,hmm);
		}else if(sigTest == 3){
			segmentHmmStateByFisherTest(loci, methyState, hiddenState, nfbsc,hmm);
		}else if(sigTest == 4){
			segmentHmmStateByBetaDiffTest(loci, methyState, hiddenState, nfbsc,hmm);
		}else{
			throw new Exception("Not such a sigTestMode");
		}
	}
	
	protected void segmentHmmState(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, BbForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm) throws Exception{

		
		if(sigTest == 1){
			segmentHmmStateByRandomPermutation(loci, methyState, hiddenState, nfbsc,hmm);
		}else if(sigTest == 2){
			segmentHmmStateByBinomialTestWithFBSC(loci, methyState, hiddenState, nfbsc,hmm);
		}else if(sigTest == 3){
			segmentHmmStateByFisherTest(loci, methyState, hiddenState, nfbsc,hmm);
		}else if(sigTest == 4){
			segmentHmmStateByBetaDiffTest(loci, methyState, hiddenState, nfbsc,hmm);
		}else{
			throw new Exception("Not such a sigTestMode");
		}
	}
	
	
	
	abstract protected void segmentHmmStateByRandomPermutation(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, ForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm);

	abstract protected void segmentHmmStateByBinomialTestWithFBSC(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, NdrForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm);
	
	abstract protected void segmentHmmStateByFisherTest(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, NdrForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm);
	
	abstract protected void segmentHmmStateByBetaDiffTest(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, ForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm);
	
	abstract protected void segmentHmmStateByRandomPermutation(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, BbForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm);

	abstract protected void segmentHmmStateByBinomialTestWithFBSC(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, BbForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm);
	
	abstract protected void segmentHmmStateByFisherTest(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, BbForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm);
	
	abstract protected void segmentHmmStateByBetaDiffTest(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, BbForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm);

	
	protected double getRandomPermutatedPvalue(double value, double[] randomScoreDistribution, boolean reverse){
		int rank = 0;
		for(int i = 0; i < randomScoreDistribution.length; i++){
			if(randomScoreDistribution[i] >= value){
				rank++;
			}
		}
		if(reverse)
			return (double)rank/(double)randomScoreDistribution.length;
		else
			return 1 - (double)rank/(double)randomScoreDistribution.length;
	}
	
	protected double getBinomialSigTest(int k, int n, double pSucess, boolean reverseP){
		
		BinomialDistributionImpl binomial = new BinomialDistributionImpl(n, pSucess);
		
		 double p = Double.NaN;
		try {
			if(reverseP){
				p = 1-binomial.cumulativeProbability(k)+binomial.probability(k); //p(X>= x)
			}
			else{
				p = binomial.cumulativeProbability(k); //p(X<= x)
			}
			
		} catch (MathException e) {
			
			e.printStackTrace();
		}
		return p;
		
	}
	
	protected double getFisherPvalue(int methy, int unmethy, int methyAdj, int unmethyAdj, boolean rightTail){
		
		double pValue = Double.NaN;
		
		if(methy+unmethy == 0 || unmethyAdj + methyAdj == 0){
			return pValue;
		}
		
		int sum = methy + methyAdj + unmethy + unmethyAdj + 100;
			FisherExactTest fisherExact = new FisherExactTest(sum);
			if(rightTail){
				pValue = fisherExact.getRightTailedP(methy,methyAdj, unmethy,  unmethyAdj);
			}
			else{
				pValue = fisherExact.getLeftTailedP(methy,methyAdj, unmethy,  unmethyAdj);
			}
			 
		 return pValue;
	}
	
	protected enum sigTestMode{
		permutation,
		binomial,
		fisher,
		betaDiff //default
	}
	
}
