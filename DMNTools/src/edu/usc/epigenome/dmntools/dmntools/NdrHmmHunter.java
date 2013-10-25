/**
 * 
 */
package edu.usc.epigenome.dmntools.dmntools;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import org.broad.tribble.annotation.Strand;


import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;


import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardScaledCalculator;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;

import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchScaledLearner;
import be.ac.ulg.montefiore.run.jahmm.toolbox.KullbackLeiblerDistanceCalculator;


import edu.usc.epigenome.dmntools.hmm.BbForwardBackwardScaledCalculator;
import edu.usc.epigenome.dmntools.hmm.BbViterbiCalculator;
import edu.usc.epigenome.dmntools.hmm.HmmReader;
import edu.usc.epigenome.dmntools.hmm.HmmWriter;
import edu.usc.epigenome.dmntools.hmm.NdrForwardBackwardScaledCalculator;
import edu.usc.epigenome.dmntools.hmm.ObservationMethy;
import edu.usc.epigenome.dmntools.hmm.OpdfBetaReader;
import edu.usc.epigenome.dmntools.hmm.OpdfBetaWriter;
import edu.usc.epigenome.dmntools.utils.BisulfiteGenomicLocHmm;
import edu.usc.epigenome.dmntools.utils.BisulfiteSlidingWindow;
import edu.usc.epigenome.dmntools.utils.GenomeLocus;
import edu.usc.epigenome.dmntools.distribution.OpdfBeta;





import edu.usc.epigenome.uecgatk.bissnp.writer.bedObject;
import gnu.trove.map.hash.THashMap;




/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Oct 16, 2013 5:12:03 PM
 * 
 */
public class NdrHmmHunter extends HmmHunter{

	/**
	 * @param args
	 */
		
		final private static String USAGE = "NdrHmmHunter [opts] outputPrefix input.6plus2.bed";

		@Argument
		private List<String> arguments = new ArrayList<String>();
		
		//option for local decoding (Fisher exact test)
		
		//option for local decoding (BetaDiff  test)

		
		private Hmm<ObservationReal> hmm = null;
		

		public static void main(String[] args)
		throws Exception
		{
			NdrHmmHunter nhh = new NdrHmmHunter();
			nhh.doMain(args);
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
			
			//initiate global list and writer
			initiate(arguments.get(0));
			
			//read the list of input signal
			parseBedFile(arguments.get(1));			
			
			//train HMM
			if(!onlyDecode){
				trainHmm();
			}
			
			
			//decoding by HMM model, also segment HMM state into different states' segment
			if(!onlyTrain){
				decodeHmm();
				//close writer
				finished();
			}

		}

		
		
		//initiate HMM && training HMM
		protected void trainHmm() throws IOException{
			System.out.println("HMM training....");
			if(kmeans){
				hmm = buildInitHmmByBeta(methyValue);
			}
			else{
				hmm = buildInitHmmRandomlyByBeta(methyValue);
			}
			
			BaumWelchScaledLearner bwl = new BaumWelchScaledLearner();

			Hmm<ObservationReal> prevHmm = null;
			try {
				prevHmm = hmm.clone();
			} catch (CloneNotSupportedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			
			// This object measures the distance between two HMMs
			KullbackLeiblerDistanceCalculator klc = 
				new KullbackLeiblerDistanceCalculator();
			
			double distance = Double.MAX_VALUE;
			int i = 0;
			// Incrementally improve the solution
			while(Math.abs(distance) >= tol){
				i++;
				
				try {
					prevHmm = hmm.clone();
				} catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				System.out.println("HMM pre:\n" + prevHmm);
				
				//hmm = bwl.iterate(hmm, result.value);
				hmm = bwl.iterate(hmm, methyValue);
				distance = klc.distance(prevHmm, hmm);
				System.out.println("Distance at iteration " + i + ": " +
						distance);
			}  

			System.out.println("Resulting HMM:\n" + hmm);

			HmmWriter.write(hmmWriter, new OpdfBetaWriter(), hmm);
			if(!onlyDecode){
				hmmWriter.close();
			}
		}
		
		//decoding HMM 
		protected void decodeHmm() throws Exception{
			System.out.println("\nDecoding ...\n");	
			
			hmm = HmmReader.read(new FileReader(hmmFile), new OpdfBetaReader());

				for(int j=0; j < methyObj.size(); j++){
					System.out.println("Decoding island " + (j+1) + " in total " + methyObj.size() + " of island. Now in the island beginning at " + position.get(j).get(0).toString() + ". Island size (# of GCH) is: " + methyObj.get(j).size());	
					
					int[] hiddenState = null;
					if(beta){
						hiddenState = hmm.mostLikelyStateSequence(methyValue.get(j));
					}
					else{
						hiddenState = (new BbViterbiCalculator(methyObj.get(j), hmm)).stateSequence();
					}

					GenomeLocus[] loci = position.get(j).toArray(new GenomeLocus[position.get(j).size()]);

					ObservationMethy[] methyState = methyObj.get(j).toArray(new ObservationMethy[methyObj.get(j).size()]);

					
					if(beta){
						NdrForwardBackwardScaledCalculator nfbsc = new NdrForwardBackwardScaledCalculator(methyValue.get(j), hmm);
						segmentHmmState(loci, methyState, hiddenState, nfbsc, hmm);
					}
					else{
						BbForwardBackwardScaledCalculator nfbsc = new BbForwardBackwardScaledCalculator(methyObj.get(j), hmm);
						segmentHmmState(loci, methyState, hiddenState, nfbsc, hmm);
					}
					
					//segment HMM state by different significant test mode
					
					
					for(int i = 0; i < hiddenState.length; i++){
						List<Object> tmp = new ArrayList<Object>();
						tmp.add(hiddenState[i]);
						tmp.add((int)(methyState[i].value*methyState[i].coverage));
						tmp.add(methyState[i].coverage);
						bedObject bedLine = new bedObject(loci[i].getChr(), loci[i].getStart(), loci[i].getEnd(), (List)tmp);
						stateWriter.add(bedLine);
					}
					//j++;
				}
		}

		
		//segment 
		
		/*******************Binomial Test to do segmentation *************************/
		@Override
		protected void segmentHmmStateByBinomialTestWithFBSC(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, NdrForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm){
			
			//TODO: find the way for multi-state
			
			double maxMean = Math.max(((OpdfBeta)hmm.getOpdf(1)).mean(), ((OpdfBeta)hmm.getOpdf(0)).mean());

			for(int i=0; i<hmm.nbStates(); i++){
				if(((OpdfBeta)hmm.getOpdf(i)).mean() >= maxMean){
					segmentHmmStateByBinomialTestWithFBSCInSingleState(loci, methyState, hiddenState, nfbsc, hmm, i, true);
				}else{
					segmentHmmStateByBinomialTestWithFBSCInSingleState(loci, methyState, hiddenState, nfbsc, hmm, i, false);
				}
				
			}
			
			
		}
		
		protected void segmentHmmStateByBinomialTestWithFBSCInSingleState(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, NdrForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm, int hmmState,boolean reverseP){
			THashMap<GenomeLocus, Integer> numCtLeftBound = new THashMap<GenomeLocus, Integer>();
			THashMap<GenomeLocus, Integer> numCLeftBound = new THashMap<GenomeLocus, Integer>();
			THashMap<GenomeLocus, Integer> numCtRightBound = new THashMap<GenomeLocus, Integer>();
			THashMap<GenomeLocus, Integer> numCRightBound = new THashMap<GenomeLocus, Integer>();
			
			getAdjacentWindowCtStatusAtEachPos(numCtLeftBound, numCLeftBound, numCtRightBound, numCRightBound, hmmState, loci, methyState, hiddenState,nfbsc);
			getSegmentAndPvalue(numCtLeftBound, numCLeftBound, numCtRightBound, numCRightBound, hmmState, loci, methyState, hiddenState,nfbsc, reverseP);
		}
		
		private void getAdjacentWindowCtStatusAtEachPos(THashMap<GenomeLocus, Integer> numCtLeftBound, THashMap<GenomeLocus, Integer> numCLeftBound, THashMap<GenomeLocus, Integer> numCtRightBound, THashMap<GenomeLocus, Integer> numCRightBound, 
				int hmmState, GenomeLocus[] loci, ObservationMethy[] observationMethyState, int[] hiddenState, NdrForwardBackwardScaledCalculator nfbsc){
			BisulfiteSlidingWindow slidingWindowLeft = new BisulfiteSlidingWindow(ctInWindow, gchInWindow, window, hmmState);
			BisulfiteSlidingWindow slidingWindowRight = new BisulfiteSlidingWindow(ctInWindow, gchInWindow, window, hmmState);
			double[] methyState = new double[observationMethyState.length];
			int[] numCTState = new int[observationMethyState.length];
			for(int i=0; i< observationMethyState.length; i++ ){
				methyState[i] = observationMethyState[i].value;
				numCTState[i] = observationMethyState[i].coverage;
			}
			
			for(int z=0; z < hiddenState.length; z++){
				BisulfiteGenomicLocHmm data = new BisulfiteGenomicLocHmm(loci[z].getChr(),loci[z].getStart(), loci[z].getEnd(), new ObservationReal(methyState[z]),numCTState[z], (int)(numCTState[z]*methyState[z]),hiddenState[z]);
				if(slidingWindowLeft.getLength() < window || slidingWindowLeft.getGchNum() < gchInWindow || slidingWindowLeft.getCtReadsNum() < ctInWindow){ // in the beginning of the chromosome
					slidingWindowLeft.addLast(data);
					if(z==hiddenState.length-1){
						int numC = slidingWindowLeft.getCReadsNum();
						int numCt = slidingWindowLeft.getCtReadsNum();
						for(int j = 0; j <= z; j++){
							numCLeftBound.put(loci[j], numC);
							numCtLeftBound.put(loci[j], numCt);
							numCRightBound.put(loci[j], 0);
							numCtRightBound.put(loci[j], 0);
						}
					}
					
				}
				else if(slidingWindowRight.windowList.isEmpty()){
					int numC = slidingWindowLeft.getCReadsNum();
					int numCt = slidingWindowLeft.getCtReadsNum();
					for(int j = 0; j < z; j++){
						numCLeftBound.put(loci[j], numC);
						numCtLeftBound.put(loci[j], numCt);
						numCRightBound.put(loci[j], 0);
						numCtRightBound.put(loci[j], 0);
					}
					slidingWindowRight.addLast(data);
					numCLeftBound.put(data.position, numC);
					numCtLeftBound.put(data.position, numCt);
					numCRightBound.put(data.position, 0);
					numCtRightBound.put(data.position, 0);

				}
				else if(slidingWindowRight.getLength() < window || slidingWindowRight.getGchNum() < gchInWindow || slidingWindowRight.getCtReadsNum() < ctInWindow){
					slidingWindowRight.addLast(data);
					
					numCLeftBound.put(data.position, slidingWindowRight.getCReadsNum());
					numCtLeftBound.put(data.position, slidingWindowRight.getCtReadsNum());
					numCRightBound.put(data.position, 0);
					numCtRightBound.put(data.position, 0);
			
					if(z==hiddenState.length-1){
						int numC = slidingWindowRight.getCReadsNum();
						int numCt = slidingWindowRight.getCtReadsNum();
						for(int j = z - slidingWindowRight.windowList.size(); j <= z; j++){
							numCRightBound.put(loci[j], numC);
							numCtRightBound.put(loci[j], numCt);
						}
					}
					
				}
				else{
					int numC = slidingWindowRight.getCReadsNum();
					int numCt = slidingWindowRight.getCtReadsNum();
					LinkedList<BisulfiteGenomicLocHmm> tmpValue = slidingWindowRight.addLast(data, true);
					numCLeftBound.put(data.position, slidingWindowRight.getCReadsNum());
					numCtLeftBound.put(data.position, slidingWindowRight.getCtReadsNum());
					numCRightBound.put(data.position, 0);
					numCtRightBound.put(data.position, 0);
					for(BisulfiteGenomicLocHmm tmpData: tmpValue){
						
						numCLeftBound.put(tmpData.position, slidingWindowLeft.getCReadsNum());
						numCtLeftBound.put(tmpData.position, slidingWindowLeft.getCtReadsNum());
						slidingWindowLeft.addLast(tmpData, true);
						numCRightBound.put(tmpData.position, numC);
						numCtRightBound.put(tmpData.position, numCt);

					}
					// first add CT reads number when new value add into slidingWindowRight, it will be updated when they are popped out.

				}
				
			}
			
		}
		
		private void getSegmentAndPvalue(THashMap<GenomeLocus, Integer> numCtLeftBound, THashMap<GenomeLocus, Integer> numCLeftBound, THashMap<GenomeLocus, Integer> numCtRightBound, THashMap<GenomeLocus, Integer> numCRightBound,
				int hmmState, GenomeLocus[] loci, ObservationMethy[] observationMethyState, int[] hiddenState, NdrForwardBackwardScaledCalculator nfbsc, boolean reverseP){
			
			double[] methyState = new double[observationMethyState.length];
			int[] numCTState = new int[observationMethyState.length];
			int[] numCState = new int[observationMethyState.length];
			for(int i=0; i< observationMethyState.length; i++ ){
				numCTState[i] = observationMethyState[i].coverage;
				numCState[i] = (int)(numCTState[i] * observationMethyState[i].value);
				methyState[i] = (double)numCState[i]/(double)numCTState[i];
			}
			
			
			String chr = null;
			int start = -1;
			int end = -1;
			double score = 0;
			int dataPoint = 0; //number of GCH in the segments
			int numC_npr = 0;
			int numCT_npr = 0;

			GenomeLocus startLoc = null;
			GenomeLocus oneBeforeStartLoc = null;
			GenomeLocus preLoc = null;
			int preState = -1; //record the previous loci's state
			
			double preWeight = 0;
			double postWeight = 0;
			
			for(int i=0; i < hiddenState.length; i++){
				if(preState == -1){
					preState = hiddenState[i]; 
					if(preState == hmmState){
						chr = loci[i].getChr();
						start = loci[i].getStart();
						startLoc = loci[i];
						oneBeforeStartLoc = loci[i];
						score = methyState[i];
						dataPoint=1;
						numC_npr = numCState[i];
						numCT_npr = numCTState[i];

					}	
					
				}
				else{
					if(preState != hmmState){
						if(hiddenState[i] == hmmState){
								chr = loci[i].getChr();
								preWeight = 1 - nfbsc.getAlpha(hmm, i, preState, hiddenState[i]);

								start =(int)( (loci[i].getStart() - preLoc.getStart()) * preWeight + preLoc.getStart()) ;

								startLoc = loci[i];
								oneBeforeStartLoc = preLoc;
								numC_npr = numCState[i];
								numCT_npr = numCTState[i];
								score = methyState[i];
								dataPoint = 1;
						}

					}
					else{
						if(hiddenState[i] == hmmState){
							numC_npr += numCState[i];
							numCT_npr += numCTState[i];
							score += methyState[i];
							dataPoint++;
						}
						else{

							postWeight = 1 - nfbsc.getAlpha(hmm, i, preState, hiddenState[i]);

							end =(int)( (loci[i].getEnd() - preLoc.getEnd()) * postWeight + preLoc.getEnd()) ;

							int numC_back = numCLeftBound.get(startLoc) + numCRightBound.get(preLoc);
							int numCT_back = numCtLeftBound.get(startLoc) + numCtRightBound.get(preLoc);
							double pValue = getBinomialSigTest(numC_npr, numCT_npr, (double)numC_back/(double)numCT_back, reverseP);
							List<Object> tmp = new LinkedList<Object>();
							tmp.add(hmmState);
							tmp.add(String.format("%.2f",100*score/dataPoint));
							tmp.add(dataPoint);
							tmp.add(numC_npr);
							tmp.add(numCT_npr);
							tmp.add(numC_back);
							tmp.add(numCT_back);
							tmp.add(start - oneBeforeStartLoc.getStart());//distance of segment start point to the previous GCH, 
							tmp.add(startLoc.getStart() - start);//distance of segment start point to the next GCH, 
							tmp.add(end - preLoc.getEnd());//distance of segment end point to the previous GCH, 
							tmp.add(loci[i].getEnd() - end);//distance of segment end point to the next GCH, 
							tmp.add(String.format("%.3f",1-preWeight));
							tmp.add(String.format("%.3f",1-postWeight));
							tmp.add(String.format("%.6f",pValue));

							if(numCT_npr >= minCT && numCT_back >= minCT){
								bedObject bedLine = new bedObject(chr, start, end, (List)tmp);
								segWriter.add(bedLine);
							}

						}
					}
					preState = hiddenState[i];
				}
				preLoc = loci[i];
			}
		}

		/* (non-Javadoc)
		 * @see edu.usc.epigenome.dmntools.dmntools.HmmHunter#segmentHmmStateByRandomPermutation(edu.usc.epigenome.dmntools.utils.GenomeLocus[], edu.usc.epigenome.dmntools.hmm.ObservationMethy[], int[], be.ac.ulg.montefiore.run.jahmm.ForwardBackwardScaledCalculator, be.ac.ulg.montefiore.run.jahmm.Hmm)
		 */
		@Override
		protected void segmentHmmStateByRandomPermutation(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, ForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm) {
			// TODO Auto-generated method stub
			
		}


		/* (non-Javadoc)
		 * @see edu.usc.epigenome.dmntools.dmntools.HmmHunter#segmentHmmStateByFisherTest(edu.usc.epigenome.dmntools.utils.GenomeLocus[], edu.usc.epigenome.dmntools.hmm.ObservationMethy[], int[], be.ac.ulg.montefiore.run.jahmm.ForwardBackwardScaledCalculator, be.ac.ulg.montefiore.run.jahmm.Hmm)
		 */
		@Override
		protected void segmentHmmStateByFisherTest(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, ForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm) {
			// TODO Auto-generated method stub
			
		}

		/* (non-Javadoc)
		 * @see edu.usc.epigenome.dmntools.dmntools.HmmHunter#segmentHmmStateByBetaDiffTest(edu.usc.epigenome.dmntools.utils.GenomeLocus[], edu.usc.epigenome.dmntools.hmm.ObservationMethy[], int[], be.ac.ulg.montefiore.run.jahmm.ForwardBackwardScaledCalculator, be.ac.ulg.montefiore.run.jahmm.Hmm)
		 */
		@Override
		protected void segmentHmmStateByBetaDiffTest(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, ForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm) {
			// TODO Auto-generated method stub
			
		}

		/* (non-Javadoc)
		 * @see edu.usc.epigenome.dmntools.dmntools.HmmHunter#segmentHmmStateByRandomPermutation(edu.usc.epigenome.dmntools.utils.GenomeLocus[], edu.usc.epigenome.dmntools.hmm.ObservationMethy[], int[], edu.usc.epigenome.dmntools.hmm.BbForwardBackwardScaledCalculator, be.ac.ulg.montefiore.run.jahmm.Hmm)
		 */
		@Override
		protected void segmentHmmStateByRandomPermutation(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, BbForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm) {
			// TODO Auto-generated method stub
			
		}

		/* (non-Javadoc)
		 * @see edu.usc.epigenome.dmntools.dmntools.HmmHunter#segmentHmmStateByBinomialTestWithFBSC(edu.usc.epigenome.dmntools.utils.GenomeLocus[], edu.usc.epigenome.dmntools.hmm.ObservationMethy[], int[], edu.usc.epigenome.dmntools.hmm.BbForwardBackwardScaledCalculator, be.ac.ulg.montefiore.run.jahmm.Hmm)
		 */
		@Override
		protected void segmentHmmStateByBinomialTestWithFBSC(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, BbForwardBackwardScaledCalculator nfbsc,
				Hmm<? extends Observation> hmm) {
			// TODO Auto-generated method stub
			
		}

		/* (non-Javadoc)
		 * @see edu.usc.epigenome.dmntools.dmntools.HmmHunter#segmentHmmStateByFisherTest(edu.usc.epigenome.dmntools.utils.GenomeLocus[], edu.usc.epigenome.dmntools.hmm.ObservationMethy[], int[], edu.usc.epigenome.dmntools.hmm.BbForwardBackwardScaledCalculator, be.ac.ulg.montefiore.run.jahmm.Hmm)
		 */
		@Override
		protected void segmentHmmStateByFisherTest(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, BbForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm) {
			// TODO Auto-generated method stub
			
		}

		/* (non-Javadoc)
		 * @see edu.usc.epigenome.dmntools.dmntools.HmmHunter#segmentHmmStateByBetaDiffTest(edu.usc.epigenome.dmntools.utils.GenomeLocus[], edu.usc.epigenome.dmntools.hmm.ObservationMethy[], int[], edu.usc.epigenome.dmntools.hmm.BbForwardBackwardScaledCalculator, be.ac.ulg.montefiore.run.jahmm.Hmm)
		 */
		@Override
		protected void segmentHmmStateByBetaDiffTest(GenomeLocus[] loci, ObservationMethy[] methyState, int[] hiddenState, BbForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm) {
			// TODO Auto-generated method stub
			
		}
		
		/*******************Fisher exact Test to do segmentation *************************/
		
		
		/*******************Beta Diff Test to do segmentation *************************/
		
		

		/*******************random permutation to do segmentation *************************/
		/****only good for HMR detection, for nucleosome positioning, it is not a good idea..***/
		/*
		protected void segmentHmmStateByRandomPermutation(GenomeLocus[] loci, Double[] methyState, int[] hiddenState, ForwardBackwardScaledCalculator nfbsc, Hmm<? extends Observation> hmm){
			double[] randomSeqs = getPermutatedRandomSeqs(result.value, hmm);
			getNDRSegment(hiddenState, methyState, loci, randomSeqs, ndrState);
		}
		
		
		private double[] getPermutatedRandomSeqs(ArrayList<LinkedList<ObservationReal>> seqs, Hmm<ObservationReal> hmm, int ndrState){
			ArrayList<ObservationReal> randomValueList = new ArrayList<ObservationReal>();
			
			Iterator<LinkedList<ObservationReal>> it = seqs.iterator();
			while(it.hasNext()){
				randomValueList.addAll(it.next());
			}
			Collections.shuffle(randomValueList);
			int[] hiddenState = hmm.mostLikelyStateSequence(randomValueList);
			double[] score = getScoreInEachStateByRandomPermutation(hiddenState, randomValueList, ndrState);
			return score;
		}
		
		private double[] getScoreInEachStateByRandomPermutation(int[] hiddenState, ArrayList<ObservationReal> randomValueList){
			ArrayList<Double> scoreList = new ArrayList<Double>();
			double score = 0;
			int preState = -1;
			int i = 0;
			Iterator<ObservationReal> it = randomValueList.iterator();
			while(it.hasNext()){
				double methyState = it.next().value;
				if(preState == -1){
					preState = hiddenState[i];
					if(preState == ndrState){
						
						score += methyState;
					}
					i++;
					continue;
				}
				if(preState != ndrState && hiddenState[i] == ndrState){

					score += methyState;
					preState = hiddenState[i];
				}
				else if(preState == ndrState){
					if(hiddenState[i] != ndrState){
						preState = hiddenState[i];
						scoreList.add(score);
						score = 0;
					}
					else{
						preState = hiddenState[i];
						score += methyState;
					}
					
				}
				i++;
			}
			double[] scoreCollection = new double[scoreList.size()];
			Iterator<Double> it2 = scoreList.iterator();
			i = 0;
			while(it2.hasNext()){
				scoreCollection[i++] = it2.next();
			}
			return scoreCollection;
		}
		
		private void segmentHmmStateByRandomPermutation(int[] hiddenState, double[] methyState, GenomeLoc[] loci, THashMap<Integer, double[]> randomSeqsScoreInDiffStates){
			String chr = null;
			int start = -1;
			int end = -1;
			double score = 0;
			int preState = -1;
			GenomeLoc preLoc = null;
			int dataPoint = 0;
			for(int i=0; i < hiddenState.length; i++){
				if(preState == -1){
					preState = hiddenState[i];
					preLoc = loci[i]; 
					if(preState == ndrState){
						chr = loci[i].getContig();
						start = loci[i].getStart()-1;
						score += methyState[i];
						dataPoint++;
					}	
					continue;
				}
				if(preState != ndrState && hiddenState[i] == ndrState){
					chr = loci[i].getContig();
					start = loci[i].getStart()-1;
					score += methyState[i];
					dataPoint++;
					preState = hiddenState[i];
					preLoc = loci[i];
				}
				else if(preState == ndrState){
					if(hiddenState[i] != ndrState){
						end = preLoc.getStart();
						preState = hiddenState[i];
						List<Object> tmp = new LinkedList<Object>();
						tmp.add(score);
						tmp.add(dataPoint);
						tmp.add(end-start);
						tmp.add(getPvalue(score, randomSeqs));
						bedObject bedLine = new bedObject(chr, start, end, Strand.NONE, (List)tmp);
						segWriter.add(bedLine);
						score = 0;
						dataPoint = 0;
						preLoc = loci[i];
					}
					else{
						preState = hiddenState[i];
						preLoc = loci[i];
						score += methyState[i];
						dataPoint++;
					}
					
				}
				
			}
		}
		*/
		

		
}
