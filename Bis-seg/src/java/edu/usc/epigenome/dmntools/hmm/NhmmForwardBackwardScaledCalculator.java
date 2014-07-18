/**
 * 
 */
package edu.usc.epigenome.dmntools.hmm;

import java.util.EnumSet;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardScaledCalculator;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Dec 12, 2012 2:03:01 PM
 * 
 */
public class NhmmForwardBackwardScaledCalculator extends ForwardBackwardScaledCalculator {

	/**
	 * @param oseq
	 * @param hmm
	 * @param flags
	 */
	public <O extends Observation> NhmmForwardBackwardScaledCalculator(List<? extends O> oseq, Hmm<O> hmm, EnumSet<Computation> flags) {
		super(oseq, hmm, flags);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param oseq
	 * @param hmm
	 */
	public <O extends Observation> NhmmForwardBackwardScaledCalculator(List<? extends O> oseq, Hmm<O> hmm) {
		super(oseq, hmm);
		// TODO Auto-generated constructor stub
	}
	
	@Override
	protected <O extends Observation> void 
	computeAlphaStep(Hmm<? super O> hmm, O o, int t, int j)
	{
		double sum = 0.;
		
		for (int i = 0; i < hmm.nbStates(); i++)
			sum += alpha[t-1][i] * ((Nhmm<ObservationMethy>)hmm).getAijx(i, j, ((ObservationMethy)o).distance);		

		alpha[t][j] = sum * hmm.getOpdf(j).probability(o);
	}
	
	/* Computes beta[t][i] (t < obs. seq.le length - 1) */
	@Override
	protected <O extends Observation> void 
	computeBetaStep(Hmm<? super O> hmm, O o, int t, int i)
	{
		double sum = 0.;
		
		for (int j = 0; j < hmm.nbStates(); j++){
			sum += beta[t+1][j] * ((Nhmm<ObservationMethy>)hmm).getAijx(i, j, ((ObservationMethy)o).distance) * 
					hmm.getOpdf(j).probability(o);
		//	System.err.print(",sum:" + sum + ",beta:" + beta[t+1][j] + ",getAij:" + hmm.getAij(i, j) + ",getOpdf:" + hmm.getOpdf(j).probability(o) + ",t:" + t + ",i:" + i + ",j:" + j + ",oberve:" + ((ObservationReal)o).value
		//			+ ",mean:" + ((OpdfGaussianMixture)hmm.getOpdf(j)).means()[j] + ",var:" + ((OpdfGaussianMixture)hmm.getOpdf(j)).variances()[j] + ",proportion:" + ((OpdfGaussianMixture)hmm.getOpdf(j)).proportions()[j] );
		//	if(hmm.getOpdf(j).probability(o) > 1)
			//	System.exit(1);
		}
		//System.err.println();
		
		beta[t][i] = sum;
	}


}
