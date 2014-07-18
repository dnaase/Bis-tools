/**
 * 
 */
package edu.usc.epigenome.dmntools.hmm;

import java.io.IOException;
import java.io.Writer;

import edu.usc.epigenome.dmntools.distribution.OpdfBeta;


import be.ac.ulg.montefiore.run.jahmm.io.OpdfWriter;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 15, 2012 7:12:28 PM
 * 
 */
public class OpdfBetaWriter extends OpdfWriter<OpdfBeta> {

	
	@Override
	public void write(Writer writer, OpdfBeta opdf) throws IOException {
		String s = "BetaOPDF [";
		
		s += opdf.alpha() + " " + opdf.beta();
			
		writer.write(s + "]\n");

	}

}
