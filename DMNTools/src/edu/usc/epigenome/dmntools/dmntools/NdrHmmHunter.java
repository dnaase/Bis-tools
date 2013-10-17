/**
 * 
 */
package edu.usc.epigenome.dmntools.dmntools;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.broad.tribble.annotation.Strand;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
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


		// when readCov<10, a>=1 will be filt out, if readCov>=10, a>10% will be filt out
		
		/**
		 * @param args
		 */
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

		
		
		
		//initiate HMM
		
		//training HMM
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
		
		protected void decodeHmm(){
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
						
					}
					else{
						BbForwardBackwardScaledCalculator nfbsc = new BbForwardBackwardScaledCalculator(methyObj.get(j), hmm);
						
					}
					
					if(sigTest == sigTestMode.permutation){
						
					}else if(sigTest == sigTestMode.binomial){
						
					}else if(sigTest == sigTestMode.fisher){
						
					}else if(sigTest == sigTestMode.betaDiff){
						
					}else{
						throw new Exception("Not such a sigTestMode");
					}
					
					getNPRSegmentByHalfLocalComparFBSC(hiddenState, methyState, loci, numCTState, numCState, ndrState, segWriter, true, hmm, nfbsc);
					getNPRSegmentByHalfLocalComparFBSC(hiddenState, methyState, loci, numCTState, numCState, nprState, nprSegWriter, false, hmm, nfbsc);
					
					//getNDRSegment(hiddenState, methyState, loci, numCTState, randomSeqs, ndrState);
					//getNPRSegmentByLocalCompar(hiddenState, methyState, loci, numCTState, ndrState, segWriter, true);
					//getNPRSegmentByLocalCompar(hiddenState, methyState, loci, numCTState, nprState, nprSegWriter, false);
					
					//getNPRSegmentByHalfLocalCompar(hiddenState, methyState, loci, numCTState, ndrState, segWriter, true);
					//getNPRSegmentByHalfLocalCompar(hiddenState, methyState, loci, numCTState, nprState, nprSegWriter, false);
					
				//	getNPRSegmentByHalfLocalCompar(hiddenState, methyState, loci, numCTState, numCState, ndrState2, segWriter2, true);
					
					
					
					for(int i = 0; i < hiddenState.length; i++){
						List<Object> tmp = new ArrayList<Object>();
						tmp.add(hiddenState[i]);
						tmp.add(numCState[i]);
						tmp.add(numCTState[i]);
						//tmp.add(result.numCLeftBound.get(loci));
						//tmp.add(result.numCtLeftBound.get(loci));
						//tmp.add(result.numCRightBound.get(loci));
						//tmp.add(result.numCtRightBound.get(loci));
						bedObject bedLine = new bedObject(loci[i].getChr(), loci[i].getStart()-1, loci[i].getEnd(), Strand.NONE, (List)tmp);
						segWriter.add(bedLine);
					}
					//j++;
				}
		}

		/* (non-Javadoc)
		 * @see edu.usc.epigenome.dmntools.dmntools.HmmHunter#segmentHmmState()
		 */
		@Override
		protected void segmentHmmState() {
			// TODO Auto-generated method stub
			
		}
		
		//initiate HMM state file writer
		
		//decoding
		
		//initiate writer
		
		//segment 
		
		//calculate p value
		
		//calculate FDR
		
		
		

}
