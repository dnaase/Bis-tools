/**
 * 
 */
package edu.usc.epigenome.dmntools.utils;


import org.broad.tribble.annotation.Strand;

import edu.usc.epigenome.uecgatk.bissnp.writer.genomeObject;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Oct 16, 2013 9:04:19 PM
 * 0 based genomic location object
 */
public class GenomeLocus implements genomeObject,Comparable<GenomeLocus>, Cloneable {

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
	 
	 public String toString(){
		 return new String(chr + "\t" + start + "\t" + end);
	 }

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(GenomeLocus that) {
		int result = 0;

        if ( this == that ) {
            result = 0;
        }
        else if(this.getChr().equalsIgnoreCase(that.getChr())){ //TODO: find a safe way to include strand information
        	if ( this.getStart() < that.getStart() ) result = -1;
            if ( this.getStart() > that.getStart() ) result = 1;
        }else{
        	result = 1;
        }

        return result;
	}
	
	@Override
    public int hashCode() {
        return start << 16 | end << 4 | chr.hashCode();
    }
	
	@Override
    public boolean equals(Object other) {
        if(other == null)
            return false;
        if(other instanceof GenomeLocus) {
        	GenomeLocus otherGenomeLoc = (GenomeLocus)other;
        	if(this.compareTo(otherGenomeLoc)==0){
        		return true;
        	}
        }
        return false;
    }
	
	@Override
    public GenomeLocus clone() {
        return new GenomeLocus(chr,start,end,strand);
    }
}
