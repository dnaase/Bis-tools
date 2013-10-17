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
 * @time Dec 12, 2012 12:28:04 AM
 * 
 */
public class Nhmm<O extends Observation> extends Hmm {

	/**
	 * @param nbStates
	 * @param opdfFactory
	 */
	public Nhmm(int nbStates, OpdfFactory opdfFactory) {
		super(nbStates, opdfFactory);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param pi
	 * @param a
	 * @param opdfs
	 */
	public Nhmm(double[] pi, double[][] a, List<? extends Opdf<O>> opdfs) {
		super(pi, a, opdfs);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param nbStates
	 */
	public Nhmm(int nbStates) {
		super(nbStates);
		// TODO Auto-generated constructor stub
	}

	//x is normalized distance(bp)
	public double getAijx(int i, int j, double x)
	{
		double sum = 0;
		for(int z=0; z < nbStates(); z++){
			sum += Math.exp(-getAij(i,z) + getAij(i,z)*x);
		}
		return Math.exp(-getAij(i,j) + getAij(i,j)*x)/sum;
	}
	
	@SuppressWarnings("unchecked")
	public Nhmm<O> clone()
	throws CloneNotSupportedException
	{
		Nhmm<O> hmm = new Nhmm<O>(nbStates());
		for(int i = 0; i < nbStates(); i++){
			hmm.setPi(i, this.getPi(i));
			for(int j = 0; j < nbStates(); j++){
				hmm.setAij(i, j, this.getAij(i, j));
			}
				hmm.setOpdf(i, this.getOpdf(i).clone());
		}	
		return hmm;
	}
	
}
