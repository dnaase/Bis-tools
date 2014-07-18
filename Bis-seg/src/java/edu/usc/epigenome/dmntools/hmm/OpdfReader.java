/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package edu.usc.epigenome.dmntools.hmm;

import java.io.IOException;
import java.io.StreamTokenizer;
import java.util.ArrayList;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.io.FileFormatException;


/**
 * Reads an observation distribution textual description.
 */
public abstract class OpdfReader<O extends Opdf<?>>
{	
	/**
	 * Returns the keyword identifying the distribution read.
	 * It must be the word beginning the distribution's description.
	 * 
	 * @return The keyword.
	 */
	public abstract String keyword();
	
	
	/**
	 * Reads an
	 * {@link be.ac.ulg.montefiore.run.jahmm.Opdf Opdf} out of a
	 * {@link java.io.StreamTokenizer}.
	 * <p>
	 * The stream tokenizer syntax table must be set according to
	 * of <code>HmmReader.initSyntaxTable(StreamTokenizer)
	 * </code> before the call to this method and reset to this state if
	 * modified before it returns.
	 *
	 * @param st A stream tokenizer.
	 * @return An Opdf.
	 */
	public abstract O read(StreamTokenizer st)
	throws IOException, FileFormatException;
	
	
	/**
	 * Reads a sequence of numbers.  The sequence is between brackets
	 * and numbers are separated by spaces.  Empty array are not allowed.
	 * 
	 * @param st The tokenizer to read the sequence from.
	 * @param length The expected length of the sequence or a strictly negative
	 *        number if it must not be checked.
	 * @return The array read.
	 */
	static protected double[] read(StreamTokenizer st, int length)
	throws IOException, FileFormatException
	{
		List<Double> l = new ArrayList<Double>();
		HmmReader.readWords(st, "[");
		while (st.nextToken() == StreamTokenizer.TT_NUMBER)
			l.add(st.nval);
		st.pushBack();
		HmmReader.readWords(st, "]");

		if (length >= 0 && l.size() != length)
			throw new FileFormatException(st.lineno(),
					"Wrong length of number sequence");
		
		if (l.size() == 0)
			throw new FileFormatException(st.lineno(),
					"Invalid empty sequence");
		
		double[] a = new double[l.size()];
		for (int i = 0; i < a.length; i++)
			a[i] = l.get(i);
		
		return a;
	}
}
