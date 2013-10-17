/**
 * 
 */
package edu.usc.epigenome.dmntools.hmm;

import java.io.IOException;
import java.io.StreamTokenizer;

import edu.usc.epigenome.dmntools.distribution.OpdfBeta;

import be.ac.ulg.montefiore.run.jahmm.io.FileFormatException;


/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 15, 2012 7:22:43 PM
 * 
 */
public class OpdfBetaReader extends OpdfReader<OpdfBeta> {

	/**
	 * 
	 */
	public OpdfBetaReader() {
		// TODO Auto-generated constructor stub
	}

	public String keyword()
	{
		return "BetaOPDF";
	}

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.jahmm.io.OpdfReader#read(java.io.StreamTokenizer)
	 */
	@Override
	public OpdfBeta read(StreamTokenizer st) throws IOException, FileFormatException {
		HmmReader.readWords(st, keyword());
		
		double[] alphaBeta = OpdfReader.read(st, 2);
		
		return new OpdfBeta(alphaBeta[0], alphaBeta[1]);
	}

}
