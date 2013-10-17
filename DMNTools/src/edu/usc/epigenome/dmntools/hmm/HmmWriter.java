/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package edu.usc.epigenome.dmntools.hmm;

import java.io.IOException;
import java.io.Writer;
import java.text.DecimalFormat;

import be.ac.ulg.montefiore.run.jahmm.*;
import be.ac.ulg.montefiore.run.jahmm.io.OpdfWriter;

/**
 * Writes a HMM to a text file compatible with {@link HmmReader}.
 */
public class HmmWriter
{
	/**
	 * Writes a HMM description.
	 * 
	 * @param writer The writer to write the HMM to.
	 * @param opdfWriter The writer used to convert the observation's
	 *        distributions of the HMMs.
	 * @param hmm The HMM to write.
	 */
	static public <O extends Observation> void 
	write(Writer writer, OpdfWriter<? extends Opdf<O>> opdfWriter, Hmm<O> hmm)
	throws IOException
	{
    	writer.write("Hmm v1.0\n\nNbStates " + hmm.nbStates() + "\n\n");
    	for (int i = 0; i < hmm.nbStates(); i++)
    		writeState(writer, opdfWriter, hmm, i);
	}
    
    
	@SuppressWarnings("unchecked") // Cannot guarantee type safety
	static private <O extends Observation, D extends Opdf<O>> void 
	writeState(Writer writer, OpdfWriter<D> opdfWriter,
			Hmm<O> hmm, int stateNb)
    throws IOException
    {
		DecimalFormat formatter = new DecimalFormat("#0.######");
		
    	writer.write("State\nPi " + formatter.format(hmm.getPi(stateNb)));
    	
    	writer.write("\nA ");
    	for (int i = 0; i < hmm.nbStates(); i++)
    		writer.write(formatter.format(hmm.getAij(stateNb, i)) + " ");
    	writer.write("\n");
    	
    	D opdf = (D) hmm.getOpdf(stateNb);
    	opdfWriter.write(writer, opdf);
    	writer.write("\n\n");
    }
}
