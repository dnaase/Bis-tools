/*
 * The MIT License (MIT)
 * Copyright (c) 2015 dnaase <Yaping Liu: lyping1986@gmail.com>

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * RandomGenomicRegionGenerater.java
 * Dec 5, 2014
 * 10:20:08 AM
 * yaping    lyping1986@gmail.com
 */
package edu.mit.compbio.qrf.java.utils;


import java.io.File;
import java.util.Iterator;

import org.apache.commons.math3.random.MersenneTwister;

import edu.unc.genomics.Interval;



/**
 *
 */
public class RandomGenomicRegionGenerator {

	private MersenneTwister generator = null;
	private FastaIndexFileReader fastaIndex = null;
	public long randomRange = Integer.MAX_VALUE;;

	
	public RandomGenomicRegionGenerator(){
		generator = new MersenneTwister();
	}
	
	public RandomGenomicRegionGenerator(int seed){
		generator = new MersenneTwister(seed);
		setRange();
	}
	
	public RandomGenomicRegionGenerator(File genomeIndex){
		generator = new MersenneTwister();
		fastaIndex =  new FastaIndexFileReader(genomeIndex);
		setRange();
	}
	 
	public RandomGenomicRegionGenerator(File genomeIndex, int seed){
		generator = new MersenneTwister(seed);
		fastaIndex =  new FastaIndexFileReader(genomeIndex);
		setRange();
	}
	
	//random single point region in the whole genome
	public Interval next(){
		long randomStart = generator.nextLong(randomRange);
		return fastaIndex.getContigAndLocation(randomStart);
		
	}
	
	//random bed region with defined length in the whole genome
	public Interval next(int length){
			long randomStart = generator.nextLong(randomRange);
			Interval bed = fastaIndex.getContigAndLocation(randomStart);
			if(bed.getStop()-1+length > fastaIndex.getContigSize(bed.getChr())){
				bed = this.next(length);
			}else{
				bed.setStop(bed.getStop()-1+length);
				
			}
			return bed;
			
	}
	
	//random single point region in the whole chromosome
	public Interval next(String chr){
		setRange(chr);
		long randomStart = generator.nextLong(randomRange);
		long randomEnd = randomStart;
		return new Interval(chr, (int)randomStart, (int)randomEnd);
	}
	
	//random bed region with defined length in the whole chromosome
	public Interval next(String chr, int length){
		setRange(chr);
		long randomStart = generator.nextLong(randomRange);
		long randomEnd = randomStart+length-1;
		if(randomEnd>randomRange) randomEnd = randomRange;
		
		return new Interval(chr, (int)randomStart, (int)randomEnd);
	}
	
	
	public void setRange(){
		randomRange = fastaIndex.genomeSize();
		//System.err.println(fastaIndex);
	}
	
	public void setRange(String chr){
		randomRange = fastaIndex.getContigSize(chr);
	}

}
