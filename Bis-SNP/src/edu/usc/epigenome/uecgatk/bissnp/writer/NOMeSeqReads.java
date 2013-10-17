/**
 * 
 */
package edu.usc.epigenome.uecgatk.bissnp.writer;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 14, 2012 4:31:38 PM
 * 
 */
public class NOMeSeqReads extends cpgReads {

	private String refContext;
	private String sampleContext;

	/**
	 * @param chr
	 * @param genomeLoc
	 * @param methyStatus
	 * @param baseQ
	 * @param strand
	 * @param readID
	 */
	public NOMeSeqReads(String chr, int genomeLoc, char methyStatus, byte baseQ, char strand, String readID) {
		super(chr, genomeLoc, methyStatus, baseQ, strand, readID);
		// TODO Auto-generated constructor stub
	}

	public NOMeSeqReads(String chr, int genomeLoc, String sampleContext, String refContext, char methyStatus, byte baseQ, char strand, String readID) {
		super(chr, genomeLoc, methyStatus, baseQ, strand, readID);
		this.sampleContext = sampleContext;
		this.refContext = refContext;
		// TODO Auto-generated constructor stub
	}

	public String getRefContext() {
		return this.refContext;
	}

	public String getSampleContext() {
		return this.sampleContext;
	}

}
