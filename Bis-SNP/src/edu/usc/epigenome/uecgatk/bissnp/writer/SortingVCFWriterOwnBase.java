package edu.usc.epigenome.uecgatk.bissnp.writer;

import java.util.Comparator;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.PriorityBlockingQueue;

import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import edu.usc.epigenome.uecgatk.bissnp.BaseUtilsMore;

public abstract class SortingVCFWriterOwnBase implements VCFWriter {

	protected static final int BEFORE_MOST_UPSTREAM_LOC = 0; // No real locus

	// index is <= 0
	protected String cachedCurrentChr = null;

	protected boolean enableDiscreteLoci = false;
	// The locus until which we are permitted to write out (inclusive)
	protected Integer mostUpstreamWritableLoc;

	// The set of chromosomes already passed over and to which it is forbidden
	// to return
	private Set<String> finishedChromosomes = null;

	// The VCFWriter to which to actually write the sorted VCF records
	private VCFWriter innerWriter = null;

	// the current queue of un-emitted records
	private PriorityBlockingQueue<VCFRecord> queue = null;

	// Should we call innerWriter.close() in close()
	private boolean takeOwnershipOfInner;

	public SortingVCFWriterOwnBase(VCFWriter innerWriter) {
		this(innerWriter, false); // by default, don't own inner
	}

	/**
	 * create a local-sorting VCF writer, given an inner VCF writer to write to
	 * 
	 * @param innerWriter
	 *            the VCFWriter to write to
	 * @param takeOwnershipOfInner
	 *            Should this Writer close innerWriter when it's done with it
	 */
	public SortingVCFWriterOwnBase(VCFWriter innerWriter, boolean takeOwnershipOfInner) {
		this.innerWriter = innerWriter;
		this.queue = new PriorityBlockingQueue<VCFRecord>(100000000, new VariantContextComparator());
		this.mostUpstreamWritableLoc = BEFORE_MOST_UPSTREAM_LOC;
		this.finishedChromosomes = new TreeSet<String>();
		this.takeOwnershipOfInner = takeOwnershipOfInner;
	}

	/**
	 * add a record to the file
	 * 
	 * @param vc
	 *            the Variant Context object
	 * @param refBase
	 *            the ref base
	 */
	@Override
	public synchronized void add(VariantContext vc) {
		/*
		 * Note that the code below does not prevent the successive add()-ing of: (chr1, 10),
		 * (chr20, 200), (chr15, 100) since there is no implicit ordering of chromosomes:
		 */
		if (cachedCurrentChr == null) {
			cachedCurrentChr = vc.getChr();
		}
		VCFRecord firstRec = queue.peek();
		if (firstRec != null && !vc.getChr().equals(firstRec.vc.getChr())) {

			finishedChromosomes.add(firstRec.vc.getChr());
			cachedCurrentChr = vc.getChr();
			stopWaitingToSort();
		}

		noteCurrentRecord(vc); // overwritten

		queue.add(new VCFRecord(vc));
	//	if(firstRec != null && mostUpstreamWritableLoc != null)
	//		System.err.println(queue.size() + "\t" + firstRec.vc.getStart() + "\t" + mostUpstreamWritableLoc + "\t" + vc.getStart());
		emitSafeRecords();
	}

	/**
	 * attempt to close the VCF file; we need to flush the queue first
	 */
	@Override
	public void close() {
		stopWaitingToSort();

		if (takeOwnershipOfInner)
			innerWriter.close();
	}

	public void emitSafeRecords() {
		emitRecords(false);
	}

	/**
	 * Gets a string representation of this object.
	 * 
	 * @return
	 */
	@Override
	public String toString() {
		return getClass().getName();
	}

	@Override
	public void writeHeader(VCFHeader header) {
		innerWriter.writeHeader(header);
	}

	protected void noteCurrentRecord(VariantContext vc) {
		// did the user break the contract by giving a record too late?
		if (mostUpstreamWritableLoc != null && (cachedCurrentChr != null && cachedCurrentChr.equalsIgnoreCase(vc.getChr()) && vc.getStart() < mostUpstreamWritableLoc) && !enableDiscreteLoci)
			throw new IllegalArgumentException("Permitted to write any record upstream of position : " + cachedCurrentChr + ": " + mostUpstreamWritableLoc + ", but a record at " + vc.getChr() + ":"
					+ vc.getStart() + " was just added.");
	}

	public void emitRecords(boolean emitUnsafe) {
		//System.err.println(queue.size() + "\t" + mostUpstreamWritableLoc);
		//if(!queue.isEmpty() && mostUpstreamWritableLoc != null){
		//	VCFRecord firstRec = queue.peek();
		//	System.err.println(queue.size() + "\t" + firstRec.vc.getStart() + "\t" + mostUpstreamWritableLoc);
		//}

		while (!queue.isEmpty()) {
			VCFRecord firstRec = queue.peek();

			if (emitUnsafe || mostUpstreamWritableLoc == null || firstRec.vc.getStart() <= mostUpstreamWritableLoc || !cachedCurrentChr.equalsIgnoreCase(firstRec.vc.getChr()) || enableDiscreteLoci) {
				queue.poll();
				innerWriter.add(firstRec.vc);
			} else {
				break;
			}
		}
	}

	private void stopWaitingToSort() {
		emitRecords(true);
		mostUpstreamWritableLoc = BEFORE_MOST_UPSTREAM_LOC;
	}

	private static class VariantContextComparator implements Comparator<VCFRecord> {
		@Override
		public int compare(VCFRecord r1, VCFRecord r2) {
		//	if (BaseUtilsMore.contigToIndex(r1.vc.getChr()) - BaseUtilsMore.contigToIndex(r1.vc.getChr()) != 0)
		//		return BaseUtilsMore.contigToIndex(r1.vc.getChr()) - BaseUtilsMore.contigToIndex(r1.vc.getChr());
			return r1.vc.getStart() - r2.vc.getStart();
		}
	}

	private static class VCFRecord {
		public VariantContext vc;

		public VCFRecord(VariantContext vc) {
			this.vc = vc;
		}
	}

}
