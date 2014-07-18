/**
 * 
 */
package edu.usc.epigenome.dmntools.distribution;

import java.util.Collection;

import edu.usc.epigenome.dmntools.hmm.ObservationMethy;
import be.ac.ulg.montefiore.run.jahmm.OpdfFactory;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 14, 2012 7:24:51 PM
 * 
 */
public class OpdfBetaBinomialFactory implements OpdfFactory<OpdfBetaBinomial> {

	public Collection<? extends ObservationMethy> seqs;
	/**
	 * 
	 */
	public OpdfBetaBinomialFactory() {
	}
	
	//public OpdfBetaBinomialFactory(Collection<? extends ObservationMethy> co) {
	//	seqs = co;
	//}

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.jahmm.OpdfFactory#factor()
	 */
	@Override
	public OpdfBetaBinomial factor() {
		// TODO Auto-generated method stub
		return new OpdfBetaBinomial();
	}
	
	public OpdfBetaBinomial factor(int nState) {
		switch(nState){
		case 0:
			return new OpdfBetaBinomial(10, 0.3,2.8);
		case 1:
			return new OpdfBetaBinomial(10, 1.5,1.5);
		default:
			return new OpdfBetaBinomial(10, 0.75,0.42);
		}
		
	}

}
