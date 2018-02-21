/**
 * 
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp.writer;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time May 10, 2012 11:46:00 PM
 * 
 */
public class SortingBedObjectWriter extends SortingFormatWriterBase {

	private int maxCachingStartDistance;

	/**
	 * @param innerWriter
	 * @param takeOwnershipOfInner
	 */
	public SortingBedObjectWriter(bedObjectWriterImp innerWriter, int maxCachingStartDistance, boolean takeOwnershipOfInner) {
		super(innerWriter, takeOwnershipOfInner);
		this.maxCachingStartDistance = maxCachingStartDistance;
		// TODO Auto-generated constructor stub
	}

	public SortingBedObjectWriter(FormatWriterBase readsWriter, int maxCachingStartDistance) {
		this((bedObjectWriterImp) readsWriter, maxCachingStartDistance, false);
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
