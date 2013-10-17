/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package edu.usc.epigenome.dmntools.hmm;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.*;
import be.ac.ulg.montefiore.run.jahmm.io.FileFormatException;



/**
 * This class can read Hidden Markov Models represented as text files.
 * The file syntax is as follows.
 * <p>
 * A '#' character induces a comment ; the rest of the line is skipped.
 * Words must be separated with a white space (space, tab or new line).
 * The file is case-sensitive.
 * The file must begin with the words <tt>Hmm</tt> and <tt>v1.0</tt>.
 * <p>
 * The next word must be <tt>NbStates</tt> followed by a number. This
 * number is the HMM's number of states.
 * <p>
 * Then comes a description of each state.  The n-th description is
 * related to the n-th state.  A state description begins with the keywords
 * <tt>State</tt> and </tt>Pi</tt>, followed by the initial probability of
 * this state.  Then comes the letter <tt>A</tt> followed by the state
 * transition probabilities separated by a space, in the right order.
 * Then comes a description of an observation distribution which depends on
 * the type of observation handled by the HMM.
 * <p>
 * The opdfs associated with all the states must have the same type.
 * A HMM description file thus looks like this:
 * <pre>
 * # A simple Hmm
 * Hmm
 * v1.0
 * NbStates 2
 * 
 * State
 * Pi 0.7
 * A 0.1 0.9
 * IntegerOPDF [ .2 .3 .4 .1 ]
 *
 * State
 * Pi 0.3
 * A 0.4 0.6
 * IntegerOPDF [ .1 .1 .1 .7 ]
 * </pre>
 * The lines starting with 'IntegerOPDF' are distributions descriptions.
 */
public class HmmReader
{	
	/**
	 * Reads a HMM from a text file.
	 * 
	 * @param reader The reader to read the HMM description from.
	 * @param opdfReader The {@link OpdfReader} used to read the observation
	 *        distributions.
	 */
	public static <O extends Observation> Hmm<O>
	read(Reader reader, OpdfReader<? extends Opdf<O>> opdfReader)
	throws IOException, FileFormatException
	{
		StreamTokenizer st = new StreamTokenizer(reader);
		initSyntaxTable(st);
		
		readWords(st, "Hmm", "v1.0", "NbStates");
		int nbStates = (int) readNumber(st);
		
		double[] pi = new double[nbStates];
		double[][] a = new double[nbStates][nbStates];
		List<Opdf<O>> opdfs = new ArrayList<Opdf<O>>(nbStates); 
		
		for (int i = 0; i < nbStates; i++) 
			readState(st, nbStates, i, pi, a, opdfs, opdfReader);
		
		return new Hmm<O>(pi, a, opdfs);
	}
	
	
	static private <O extends Observation> void
	readState(StreamTokenizer st, int nbStates,	int stateNb, double[] pi,
			double[][] a,
			List<Opdf<O>> opdfs, OpdfReader<? extends Opdf<O>> opdfReader)
	throws IOException, FileFormatException
	{
		readWords(st, "State", "Pi");
		pi[stateNb] = readNumber(st);
		
		readWords(st, "A");
		for (int i = 0; i < nbStates; i++)
			a[stateNb][i] = readNumber(st);
		
		opdfs.add(opdfReader.read(st));
	}
	
	
	/**
	 * Reads some keywords out of a {@link StreamTokenizer}. 
	 * 
	 * @param st A stream tokenizer.
	 * @param words The words to read, in the right order.
	 */
	static void readWords(StreamTokenizer st, String... words)
	throws IOException, FileFormatException
	{
		for (String word : words) {
			st.nextToken();
			
			if (st.ttype == StreamTokenizer.TT_WORD)
				if (st.sval.equals(word))
					continue;
				else
					throw new FileFormatException(st.lineno(),
							"Syntax error: unexpected token '" +
							st.sval + "', ('" + word + "' expected)");
			
			if (st.ttype > 0) // Single character token
				if (word.length() == 1 && st.ttype == (int) word.charAt(0))
					continue;
				else
					throw new FileFormatException(st.lineno(),
							"Syntax error: unexpected token '" +
							(char) st.ttype + "' (" + word + "' expected)");
			
			throw new FileFormatException(st.lineno(), "Syntax error: '" +
					word + "' expected");
		}
	}
	

	static double readNumber(StreamTokenizer st)
	throws IOException, FileFormatException
	{
		st.nextToken();
		
		if (st.ttype != StreamTokenizer.TT_NUMBER)
			throw new FileFormatException(st.lineno(),
					"Syntax error: number expected");
		
		return st.nval;
	}

	
	/* Initialize the syntax table of a stream tokenizer */
	static void initSyntaxTable(StreamTokenizer st)
	{
		st.resetSyntax();
		st.parseNumbers();
		st.wordChars('a', 'z');
		st.wordChars('A', 'Z');
		st.whitespaceChars(0, (int) ' ');
		st.whitespaceChars((int) '\t', (int) '\t');
		st.eolIsSignificant(false);
		st.commentChar((int) '#');
	}
}
