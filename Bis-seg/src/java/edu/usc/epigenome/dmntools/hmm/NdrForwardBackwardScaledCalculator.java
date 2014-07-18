/**
 * 
 */
package edu.usc.epigenome.dmntools.hmm;

import java.util.EnumSet;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardScaledCalculator;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardCalculator.Computation;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jan 28, 2013 2:37:26 PM
 * 
 */
public class NdrForwardBackwardScaledCalculator extends ForwardBackwardScaledCalculator {

	/**
	 * @param oseq
	 * @param hmm
	 * @param flags
	 */
	public <O extends Observation> NdrForwardBackwardScaledCalculator(List<? extends O> oseq, Hmm<O> hmm, EnumSet<Computation> flags) {
		super(oseq, hmm, flags);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param oseq
	 * @param hmm
	 */
	public <O extends Observation> NdrForwardBackwardScaledCalculator(List<? extends O> oseq, Hmm<O> hmm) {
		super(oseq, hmm);
		// TODO Auto-generated constructor stub
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
}
