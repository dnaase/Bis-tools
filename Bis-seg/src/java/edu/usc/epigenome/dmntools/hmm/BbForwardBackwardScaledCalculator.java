/**
 * 
 */
package edu.usc.epigenome.dmntools.hmm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;

import edu.usc.epigenome.dmntools.distribution.OpdfBeta;
import edu.usc.epigenome.dmntools.distribution.OpdfBetaBinomial;


import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardScaledCalculator;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardCalculator.Computation;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Feb 1, 2013 10:45:11 PM
 * 
 */
public class BbForwardBackwardScaledCalculator{

	private double[] ctFactors;
	private double lnProbability;
	
	/* alpha[t][i] = P(O(1), O(2),..., O(t+1), i(t+1) = i+1 | hmm), that is the
	 probability of the beginning of the state sequence (up to time t+1)
	 with the (t+1)th state being i+1. */
	protected double[][] alpha = null;
	protected double[][] beta = null;
	protected double probability;

	
	/**
	 * @param arrayList
	 * @param hmm
	 */
	public BbForwardBackwardScaledCalculator(ArrayList<ObservationMethy> oseq, Hmm<ObservationMethy> hmm) {
		// TODO Auto-generated constructor stub
		if (oseq.isEmpty())
			throw new IllegalArgumentException();
		
		ctFactors = new double[oseq.size()];
		Arrays.fill(ctFactors, 0.);
		
		computeAlpha(hmm, oseq);
		
		
		computeProbability(oseq, hmm);
	}

	/* Computes alpha[0][i] */

	protected <O extends Observation> void
	computeAlphaInit(Hmm<ObservationMethy> hmm, ObservationMethy observationMethy, int i)
	{
		OpdfBeta ob = (OpdfBeta) hmm.getOpdf(i);
		int coverage = ((ObservationMethy)(observationMethy)).coverage;
		OpdfBetaBinomial obb = new OpdfBetaBinomial(coverage, ob.alpha(), ob.beta());
		alpha[0][i] = hmm.getPi(i) * obb.probability(((ObservationMethy)(observationMethy)));
		//System.err.println(observationMethy.value + "\t" + obb.probability(observationMethy) + "\t" + hmm.getPi(i) + "\t" + alpha[0][i] + "\t" + i);
	}
	
	
	/* Computes alpha[t][j] (t > 0) */
	
	protected <O extends Observation> void 
	computeAlphaStep(Hmm<ObservationMethy> hmm, O observation, int t, int j)
	{
		double sum = 0.;
		
		for (int i = 0; i < hmm.nbStates(); i++)
			sum += alpha[t-1][i] * hmm.getAij(i, j);		

		OpdfBeta ob = (OpdfBeta) hmm.getOpdf(j);
		int coverage = ((ObservationMethy)(observation)).coverage;
		OpdfBetaBinomial obb = new OpdfBetaBinomial(coverage, ob.alpha(), ob.beta());
		//System.err.println(((ObservationMethy)observation).value + "\t" + obb.probability(((ObservationMethy)(observation))));
		alpha[t][j] = sum * obb.probability(((ObservationMethy)(observation)));
	}
	
	/* Computes beta[t][i] (t < obs. seq.le length - 1) */
	
	protected <O extends Observation> void 
	computeBetaStep(Hmm<? super O> hmm, O o, int t, int i)
	{
		double sum = 0.;
		
		for (int j = 0; j < hmm.nbStates(); j++){
			OpdfBeta ob = (OpdfBeta) hmm.getOpdf(j);
			int coverage = ((ObservationMethy)(o)).coverage;
			OpdfBetaBinomial obb = new OpdfBetaBinomial(coverage, ob.alpha(), ob.beta());
			
			sum += beta[t+1][j] * hmm.getAij(i, j) * 
					obb.probability(((ObservationMethy)(o)));

		}

		
		beta[t][i] = sum;
	}
	
	/* Computes the content of the scaled alpha array */
	protected <O extends Observation> void
	computeAlpha(Hmm<ObservationMethy> hmm, ArrayList<ObservationMethy> oseq)
	{	
		alpha = new double[oseq.size()][hmm.nbStates()];
		
		for (int i = 0; i < hmm.nbStates(); i++)
			computeAlphaInit(hmm, oseq.get(0), i);
		scale(ctFactors, alpha, 0);
		
		Iterator<ObservationMethy> seqIterator = oseq.iterator();
		if (seqIterator.hasNext())
			seqIterator.next();
		
		for (int t = 1; t < oseq.size(); t++) {
			ObservationMethy observation = seqIterator.next();
			
			for (int i = 0; i < hmm.nbStates(); i++)
				computeAlphaStep(hmm, observation, t, i);
			scale(ctFactors, alpha, t);
		}
	}
	
	
	/* Computes the content of the scaled beta array.  The scaling factors are
	 those computed for alpha. */
	protected <O extends Observation> void 
	computeBeta(Hmm<? super O> hmm, List<O> oseq)
	{	
		beta = new double[oseq.size()][hmm.nbStates()];
		
		for (int i = 0; i < hmm.nbStates(); i++)
			beta[oseq.size()-1][i] = 1. / ctFactors[oseq.size()-1];
		
		for (int t = oseq.size() - 2; t >= 0; t--)
			for (int i = 0; i < hmm.nbStates(); i++) {
				computeBetaStep(hmm, oseq.get(t+1), t, i);
				beta[t][i] /= ctFactors[t];
			}
	}
	
	
	/* Normalize alpha[t] and put the normalization factor in ctFactors[t] */
	private void scale(double[] ctFactors, double[][] array, int t)
	{
		double[] table = array[t];
		double sum = 0.;
		
		for (int i = 0; i < table.length; i++)
			sum += table[i];
		
		ctFactors[t] = sum;
		for (int i = 0; i < table.length; i++) 
			table[i] /= sum;
	}
	
	
	private <O extends Observation> void
	computeProbability(ArrayList<ObservationMethy> oseq, Hmm<ObservationMethy> hmm)
	{	
		lnProbability = 0.;
		
		for (int t = 0; t < oseq.size(); t++)
			lnProbability += Math.log(ctFactors[t]);
		
		probability = Math.exp(lnProbability);
	}
	
	
	/**
	 * Return the neperian logarithm of the probability of the sequence that
	 * generated this object.
	 *
	 * @return The probability of the sequence of interest's neperian logarithm.
	 */
	public double lnProbability()
	{
		return lnProbability;
	}
	
	public <O> double getAlpha(Hmm<? extends O> hmm, int t, int Prestate, int state){
		//if(t == 0){
		//	return alpha[0][state];
		//}
		//else{
			//return alpha[t-1][Prestate] * hmm.getAij(Prestate, state);
			return alpha[t][state];
		//}
		
	}
	public <O> double getBeta(Hmm<? extends O> hmm, int t, int Prestate, int state){
		
			//return alpha[t-1][Prestate] * hmm.getAij(Prestate, state);
			return beta[t][Prestate];

		
	}
	
	public double getAlpha( int t,int state){

			return alpha[t][state];

		
	}
	public double getBeta(int t, int Prestate){

			return beta[t][Prestate];

		
	}
	
	public double[] getProbFactor(){
		return ctFactors;
	}

}
