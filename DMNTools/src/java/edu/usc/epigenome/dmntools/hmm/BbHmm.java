/**
 * 
 */
package edu.usc.epigenome.dmntools.hmm;

import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfFactory;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Dec 11, 2012 4:09:33 PM
 * 
 */
public class BbHmm<O extends Observation> extends Hmm {

	/**
	 * @param nbStates
	 * @param opdfFactory
	 */
	public BbHmm(int nbStates, OpdfFactory opdfFactory) {
		super(nbStates, opdfFactory);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param pi
	 * @param a
	 * @param opdfs
	 */
	public BbHmm(double[] pi, double[][] a, List opdfs) {
		super(pi, a, opdfs);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param nbStates
	 */
	public BbHmm(int nbStates) {
		super(nbStates);
		// TODO Auto-generated constructor stub
	}

	public Opdf<O> getOpdf(int stateNb, int trials)
	{
		
		return super.getOpdf(stateNb);
	}
	
}
