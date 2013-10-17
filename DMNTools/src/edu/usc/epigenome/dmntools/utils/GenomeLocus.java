/**
 * 
 */
package edu.usc.epigenome.dmntools.utils;

import java.util.List;

import org.broad.tribble.annotation.Strand;
import org.broadinstitute.sting.utils.GenomeLoc;

import edu.usc.epigenome.uecgatk.bissnp.writer.genomeObject;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Oct 16, 2013 9:04:19 PM
 * 
 */
public class GenomeLocus implements genomeObject {

	private String chr;
	private int end;
	private int start;
	private Strand strand = Strand.NONE;

	public GenomeLocus(String chr, int start, int end) {
		// TODO Auto-generated constructor stub
		this.chr = chr;
		this.start = start;
		this.end = end;

	}

	public GenomeLocus(String chr, int start, int end, Strand strand) {
		// TODO Auto-generated constructor stub
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.strand = strand;
	}

	@Override
	public String getChr() {
		// TODO Auto-generated method stub
		return chr;
	}

	public int getEnd() {
		// TODO Auto-generated method stub
		return end;
	}

	@Override
	public int getStart() {
		// TODO Auto-generated method stub
		return start;
	}

	public Strand getStrand() {
		// TODO Auto-generated method stub
		return strand;
	}

	public String getStrandAsString() {
		// TODO Auto-generated method stub
		return strand.toString();
	}

	 public final boolean onSameContig(GenomeLocus that) {
	        return (this.chr.equalsIgnoreCase(that.getChr()));
	   }
	 
	 public final int distance( final GenomeLocus that ) {
	        if ( this.onSameContig(that) )
	            return Math.abs(this.getStart() - that.getStart());
	        else
	            return Integer.MAX_VALUE;
	   }

}
