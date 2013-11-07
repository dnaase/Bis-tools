/**
 * 
 */
package edu.usc.epigenome.dmntools.hmm;

import java.io.IOException;
import java.io.StreamTokenizer;

import be.ac.ulg.montefiore.run.jahmm.io.FileFormatException;

import edu.usc.epigenome.dmntools.distribution.OpdfBetaBinomial;


/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Dec 10, 2012 9:44:25 PM
 * 
 */
public class OpdfBetaBinomReader extends OpdfReader<OpdfBetaBinomial> {

	private int trials;
	
	public OpdfBetaBinomReader() {
		// TODO Auto-generated constructor stub
	}
	
	public OpdfBetaBinomReader(int trials) {
		this.trials = trials;
	}
	
	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.jahmm.io.OpdfReader#keyword()
	 */
	@Override
	public String keyword() {
		return "BetaOPDF"; // just because our special case..
	}

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.jahmm.io.OpdfReader#read(java.io.StreamTokenizer)
	 */
	@Override
	public OpdfBetaBinomial read(StreamTokenizer st) throws IOException, FileFormatException {
		HmmReader.readWords(st, keyword());
		
		double[] alphaBeta = OpdfReader.read(st, 2);
		
		return new OpdfBetaBinomial(trials, alphaBeta[0], alphaBeta[1]);
	}
	
	public void setTrials(int num){
		trials = num;
	}

}
