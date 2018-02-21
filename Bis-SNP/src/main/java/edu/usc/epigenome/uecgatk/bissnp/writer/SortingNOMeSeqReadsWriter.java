/**
 * 
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp.writer;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 14, 2012 6:20:58 PM
 * 
 */
public class SortingNOMeSeqReadsWriter extends SortingFormatWriterBase {

	private int maxCachingStartDistance;

	public SortingNOMeSeqReadsWriter(FormatWriterBase readsWriter, int maxCachingStartDistance) {
		this(readsWriter, maxCachingStartDistance, false);
		// TODO Auto-generated constructor stub
	}


	public SortingNOMeSeqReadsWriter(FormatWriterBase innerWriter, int maxCachingStartDistance, boolean takeOwnershipOfInner) {
		super(innerWriter, takeOwnershipOfInner);
		this.maxCachingStartDistance = maxCachingStartDistance;
		// TODO Auto-generated constructor stub
	}

	@Override
	protected void noteCurrentRecord(genomeObject obj) {
		super.noteCurrentRecord(obj); // first, check for errors

		// then, update mostUpstreamWritableLoc:
		int mostUpstreamWritableIndex = obj.getStart() - maxCachingStartDistance;
		this.mostUpstreamWritableLoc = Math.max(BEFORE_MOST_UPSTREAM_LOC, mostUpstreamWritableIndex);
	}

}
