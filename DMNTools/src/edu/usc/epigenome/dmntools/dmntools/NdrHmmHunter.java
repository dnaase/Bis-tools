/**
 * 
 */
package edu.usc.epigenome.dmntools.dmntools;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.broad.tribble.annotation.Strand;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardScaledCalculator;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.io.OpdfReader;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchScaledLearner;
import be.ac.ulg.montefiore.run.jahmm.toolbox.KullbackLeiblerDistanceCalculator;


import edu.usc.epigenome.dmntools.hmm.BbForwardBackwardScaledCalculator;
import edu.usc.epigenome.dmntools.hmm.BbViterbiCalculator;
import edu.usc.epigenome.dmntools.hmm.HmmReader;
import edu.usc.epigenome.dmntools.hmm.HmmWriter;
import edu.usc.epigenome.dmntools.hmm.NdrForwardBackwardScaledCalculator;
import edu.usc.epigenome.dmntools.hmm.OpdfBetaReader;
import edu.usc.epigenome.dmntools.hmm.OpdfBetaWriter;
import edu.usc.epigenome.dmntools.utils.GenomeLocus;
import edu.usc.epigenome.dmntools.distribution.OpdfBeta;
import edu.usc.epigenome.dmntools.dmntools.HmmHunter.sigTestMode;


import edu.usc.epigenome.uecgatk.bissnp.writer.bedObject;
import edu.usc.epigenome.uecgatk.bissnp.writer.bedObjectWriterImp;
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
				if (args.length < 2) throw new CmdLineException(USAGE);
			
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
			initiate(args[0]);
			
			//read the list of input signal
			parseBedFile(args[1]);			
			
			//train HMM
			trainHmm();
			
			//decoding by HMM model, also segment HMM state into different states' segment
			decodeHmm();
			
			//close writer
			finished();
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

		}
		
		//decoding HMM 
		protected void decodeHmm() throws Exception{
				hmm = HmmReader.read(new FileReader(hmmFile), new OpdfBetaReader());

				for(int j=0; j < methyObj.size(); j++){
					int[] hiddenState = null;
					if(beta){
						hiddenState = hmm.mostLikelyStateSequence(methyValue.get(j));
					}
					else{
						hiddenState = (new BbViterbiCalculator(methyObj.get(j), hmm)).stateSequence();
					}

					GenomeLocus[] loci = position.get(j).toArray(new GenomeLocus[position.get(j).size()]);

					Double[] methyState = methyObj.get(j).toArray(new Double[methyObj.get(j).size()]);
					
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
						bedObject bedLine = new bedObject(loci[i].getChr(), loci[i].getStart()-1, loci[i].getEnd(), Strand.NONE, (List)tmp);
						stateWriter.add(bedLine);
					}
					//j++;
				}
		}

		
		//segment 
		
		/*******************random permutation to do segmentation *************************/
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
		
		
		
		/*******************Binomial Test to do segmentation *************************/
		
		
		
		/*******************Fisher exact Test to do segmentation *************************/
		
		
		/*******************Beta Diff Test to do segmentation *************************/
		
		//calculate p value
		
		//calculate FDR
		
}
