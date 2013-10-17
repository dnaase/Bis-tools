package edu.usc.epigenome.uecgatk.bissnp.writer;

import java.util.Comparator;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.PriorityBlockingQueue;

public abstract class SortingFormatWriterBase {
	protected static final int BEFORE_MOST_UPSTREAM_LOC = 0; // No real locus
																// index is <= 0

	// The locus until which we are permitted to write out (inclusive)
	protected Integer mostUpstreamWritableLoc;

	// The set of chromosomes already passed over and to which it is forbidden
	// to return
	private Set<String> finishedChromosomes = null;
	private FormatWriterBase innerWriter = null;

	// the current queue of un-emitted records
	private PriorityBlockingQueue<genomeObject> queue = null;

	// Should we call innerWriter.close() in close()
	private boolean takeOwnershipOfInner;

	public SortingFormatWriterBase(FormatWriterBase innerWriter) {
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
	public SortingFormatWriterBase(FormatWriterBase innerWriter, boolean takeOwnershipOfInner) {
		this.innerWriter = innerWriter;
		this.queue = new PriorityBlockingQueue<genomeObject>(1000, new GenomeObjectComparator());
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
	public synchronized void add(genomeObject go) {
		/*
		 * Note that the code below does not prevent the successive add()-ing of: (chr1, 10),
		 * (chr20, 200), (chr15, 100) since there is no implicit ordering of chromosomes:
		 */
		genomeObject firstRec = queue.peek();
		if (firstRec != null && !go.getChr().equals(firstRec.getChr())) {
			if (finishedChromosomes.contains(go.getChr()))
				throw new IllegalArgumentException("Added a record at " + go.getChr() + ":" + go.getStart() + ", but already finished with chromosome" + go.getStart());

			finishedChromosomes.add(firstRec.getChr());
			stopWaitingToSort();
		}

		noteCurrentRecord(go); // possibly overwritten

		queue.add(go);
		emitSafeRecords();
	}

	/**
	 * attempt to close the VCF file; we need to flush the queue first
	 */
	public void close() {
		stopWaitingToSort();

		if (takeOwnershipOfInner)
			innerWriter.close();
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

	public void writeHeader(Object header) {
		innerWriter.addHeader(header);
	}

	public void writerFlush() {
		this.innerWriter.writerFlush();
	}

	protected void emitSafeRecords() {
		emitRecords(false);
	}

	protected void noteCurrentRecord(genomeObject obj) {
		// did the user break the contract by giving a record too late?
		if (mostUpstreamWritableLoc != null && obj.getStart() < mostUpstreamWritableLoc)
			throw new IllegalArgumentException("Permitted to write any record upstream of position " + mostUpstreamWritableLoc + ", but a record at " + obj.getChr() + ":" + obj.getStart()
					+ " was just added.");
	}

	private void emitRecords(boolean emitUnsafe) {
		while (!queue.isEmpty()) {
			genomeObject firstRec = queue.peek();

			// No need to wait, waiting for nothing, or before what we're
			// waiting for:
			if (emitUnsafe || mostUpstreamWritableLoc == null || firstRec.getStart() <= mostUpstreamWritableLoc) {
				queue.poll();

				innerWriter.add(firstRec);

			} else {
				break;
			}
		}
	}

	private void stopWaitingToSort() {
		emitRecords(true);
		mostUpstreamWritableLoc = BEFORE_MOST_UPSTREAM_LOC;
	}

	private static class GenomeObjectComparator implements Comparator<genomeObject> {

		@Override
		public int compare(genomeObject go1, genomeObject go2) {
			// TODO Auto-generated method stub
			return go1.getStart() - go2.getStart();
		}
	}
}
