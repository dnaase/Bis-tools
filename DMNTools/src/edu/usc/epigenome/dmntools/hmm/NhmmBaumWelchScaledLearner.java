/**
 * 
 */
package edu.usc.epigenome.dmntools.hmm;

import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardCalculator;
import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardScaledCalculator;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchScaledLearner;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Dec 12, 2012 1:55:49 PM
 * 
 */
public class NhmmBaumWelchScaledLearner extends BaumWelchScaledLearner {

	/**
	 * 
	 */
	public NhmmBaumWelchScaledLearner() {
		// TODO Auto-generated constructor stub
	}

	@Override
	protected <O extends Observation> ForwardBackwardCalculator
	generateForwardBackwardCalculator(List<? extends O> sequence,
			Hmm<O> hmm)
	{
		return new NhmmForwardBackwardScaledCalculator(sequence, hmm, 
				EnumSet.allOf(NhmmForwardBackwardScaledCalculator.Computation.class));
	}
	
	
	/* Here, the xi (and, thus, gamma) values are not divided by the
	 probability of the sequence because this probability might be
	 too small and induce an underflow. xi[t][i][j] still can be
	 interpreted as P[q_t = i and q_(t+1) = j | obsSeq, hmm] because
	 we assume that the scaling factors are such that their product
	 is equal to the inverse of the probability of the sequence. */
	@Override
	protected <O extends Observation> double[][][]
	estimateXi(List<? extends O> sequence, ForwardBackwardCalculator fbc,
			Hmm<O> hmm)
	{	
		if (sequence.size() <= 1)
			throw new IllegalArgumentException("Observation sequence too " + 
			"short");
		
		double xi[][][] = 
			new double[sequence.size() - 1][hmm.nbStates()][hmm.nbStates()];
		
		Iterator<? extends O> seqIterator = sequence.iterator();
		seqIterator.next();
		
		for (int t = 0; t < sequence.size() - 1; t++) {
			O observation = seqIterator.next();
			
			for (int i = 0; i < hmm.nbStates(); i++)
				for (int j = 0; j < hmm.nbStates(); j++)
					xi[t][i][j] = fbc.alphaElement(t, i) *
					((Nhmm<ObservationMethy>)hmm).getAijx(i, j, ((ObservationMethy)observation).distance) * 
					hmm.getOpdf(j).probability(observation) *
					fbc.betaElement(t + 1, j);
		}
		
		return xi;
	}
}
