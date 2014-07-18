/**
 * 
 */
package edu.usc.epigenome.dmntools.utils;


import edu.usc.epigenome.dmntools.hmm.ObservationMethy;

import be.ac.ulg.montefiore.run.jahmm.ObservationReal;


/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Oct 16, 2013 6:49:36 PM
 * 
 */
public class BisulfiteGenomicLocHmm{
	public ObservationMethy value;
	public GenomeLocus position;
	public String contig;
	public int start;
	public int end;
	public ObservationMethy methy;
	public int numCT;
	public int numC;
	public int hmmState;

	public BisulfiteGenomicLocHmm(String contig, int start, int end, ObservationMethy value, int numCT){
		this.contig = contig;
		this.start = start;
		this.end = end;
		this.value = value;
		this.numCT = numCT;
		
	}
	
	public BisulfiteGenomicLocHmm(String contig, int start, int end, ObservationMethy value, int numCT, int numC){
		this.contig = contig;
		this.start = start;
		this.end = end;
		this.value = value;
		this.numCT = numCT;
		this.numC = numC;
		this.methy = new ObservationMethy(value.value);
		methy.setCoverage(numCT);
	}
	
	public BisulfiteGenomicLocHmm(String contig, int start, int end, ObservationMethy value, int numCT, int numC, int hmmState){
		this.contig = contig;
		this.start = start;
		this.end = end;
		this.value = value;
		this.numCT = numCT;
		this.numC = numC;
		this.hmmState = hmmState;
		this.position = new GenomeLocus(contig, start,end);
	}
	
	public BisulfiteGenomicLocHmm(String contig, int start, int end, int numCT, int numC){
		ObservationMethy v = new ObservationMethy((double)numC/(double)numCT);
		new BisulfiteGenomicLocHmm(contig,  start,  end,  v,  numCT,  numC);
	}

}
