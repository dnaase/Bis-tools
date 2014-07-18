/**
 * 
 */
package edu.usc.epigenome.dmntools.hmm;

import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardScaledCalculator;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.toolbox.KullbackLeiblerDistanceCalculator;
import be.ac.ulg.montefiore.run.jahmm.toolbox.MarkovGenerator;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Dec 12, 2012 2:17:10 PM
 * 
 */
public class NhmmKullbackLeiblerDistanceCalculator extends KullbackLeiblerDistanceCalculator {

	/**
	 * 
	 */
	public NhmmKullbackLeiblerDistanceCalculator() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * Computes the Kullback-Leibler distance between two HMMs.
	 *
	 * @param hmm1 The first HMM against which the distance is computed.
	 *             The distance is mesured with regard to this HMM (this must
	 *             be defined since the Kullback-Leibler distance is not
	 *             symetric).
	 * @param hmm2 The second HMM against which the distance is computed.
	 * @return The distance between <code>hmm1</code> and <code>hmm2</code> with
	 *      regard to <code>hmm1</code>
	 */
	@Override
	public <O extends Observation> double 
	distance(Hmm<O> hmm1, Hmm<? super O> hmm2)
	{			
		double distance = 0.;
		
		for (int i = 0; i < getNbSequences(); i++) {
			
			List<O> oseq = new MarkovGenerator<O>(hmm1).
			observationSequence(getSequencesLength());
			
			distance += (new NhmmForwardBackwardScaledCalculator(oseq, hmm1).
					lnProbability() -
					new NhmmForwardBackwardScaledCalculator(oseq, hmm2).
					lnProbability()) / getSequencesLength();
			//System.err.println("distance:" + distance + new ForwardBackwardScaledCalculator(oseq, hmm1).
		//			lnProbability() + "\t" + new ForwardBackwardScaledCalculator(oseq, hmm2).
			//		lnProbability() + "\tlength:" + sequencesLength);
		}
		
		return distance / getNbSequences();
	}

	
}
