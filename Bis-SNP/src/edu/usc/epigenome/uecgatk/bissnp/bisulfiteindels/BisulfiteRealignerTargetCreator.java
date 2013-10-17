/**
 * just simply adapted from GATK's RealignerTargetCreator
 */
package edu.usc.epigenome.uecgatk.bissnp.bisulfiteindels;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadCigarFilter;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroFilter;
import org.broadinstitute.sting.gatk.filters.Platform454Filter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Allows;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.ExtendedEventPileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import edu.usc.epigenome.uecgatk.bissnp.BaseUtilsMore;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time May 3, 2012 2:58:02 PM
 * 
 */
@ReadFilters({ Platform454Filter.class, MappingQualityZeroFilter.class, BadCigarFilter.class })
@Reference(window = @Window(start = -1, stop = 50))
@Allows(value = { DataSource.READS, DataSource.REFERENCE })
@By(DataSource.REFERENCE)
@BAQMode(ApplicationTime = BAQ.ApplicationTime.FORBIDDEN)
public class BisulfiteRealignerTargetCreator extends RodWalker<BisulfiteRealignerTargetCreator.Event, BisulfiteRealignerTargetCreator.EventPair> implements
		TreeReducible<BisulfiteRealignerTargetCreator.EventPair> {

	/**
	 * Any number of VCF files representing known SNPs and/or indels. Could be e.g. dbSNP and/or
	 * official 1000 Genomes indel calls. SNPs in these files will be ignored unless the
	 * --mismatchFraction argument is used.
	 */
	@Input(fullName = "known", shortName = "known", doc = "Input VCF file with known indels", required = false)
	public List<RodBinding<VariantContext>> known = Collections.emptyList();

	/**
	 * Because the realignment algorithm is N^2, allowing too large an interval might take too long
	 * to completely realign.
	 */
	@Argument(fullName = "maxIntervalSize", shortName = "maxInterval", doc = "maximum interval size", required = false)
	protected int maxIntervalSize = 500;

	@Argument(fullName = "minReadsAtLocus", shortName = "minReads", doc = "minimum reads at a locus to enable using the entropy calculation", required = false)
	protected int minReadsAtLocus = 4;

	/**
	 * To disable this behavior, set this value to <= 0 or > 1. This feature is really only
	 * necessary when using an ungapped aligner (e.g. MAQ in the case of single-end read data) and
	 * should be used in conjunction with '--model USE_SW' in the IndelRealigner.
	 */
	@Argument(fullName = "mismatchFraction", shortName = "mismatch", doc = "fraction of base qualities needing to mismatch for a position to have high entropy", required = false)
	protected double mismatchThreshold = 0.0;

	/**
	 * The target intervals for realignment.
	 */
	@Output
	protected PrintStream out;

	/**
	 * Any two SNP calls and/or high entropy positions are considered clustered when they occur no
	 * more than this many basepairs apart.
	 */
	@Argument(fullName = "windowSize", shortName = "window", doc = "window size for calculating entropy or SNP clusters", required = false)
	protected int windowSize = 10;

	static private boolean canBeMerged(Event left, Event right) {
		return left.loc.getContigIndex() == right.loc.getContigIndex() && left.furthestStopPos >= right.loc.getStart();
	}

	// @com.google.java.contract.Requires({"left != null", "right != null"})
	static private Event mergeEvents(Event left, Event right) {
		left.merge(right);
		return left;
	}

	@Override
	public boolean generateExtendedEvents() {
		return true;
	}

	@Override
	public boolean includeReadsWithDeletionAtLoci() {
		return true;
	}

	@Override
	public void initialize() {
		if (windowSize < 2)
			throw new UserException.BadArgumentValue("windowSize", "Window Size must be an integer greater than 1");
	}

	@Override
	public Event map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

		boolean hasIndel = false;
		boolean hasInsertion = false;
		boolean hasPointEvent = false;

		int furthestStopPos = -1;

		// look for insertions in the extended context (we'll get deletions from
		// the normal context)
		if (context.hasExtendedEventPileup()) {
			ReadBackedExtendedEventPileup pileup = context.getExtendedEventPileup();
			if (pileup.getNumberOfInsertions() > 0) {
				hasIndel = hasInsertion = true;
				// check the ends of the reads to see how far they extend
				for (ExtendedEventPileupElement p : pileup.toExtendedIterable())
					furthestStopPos = Math.max(furthestStopPos, p.getRead().getAlignmentEnd());
			}
		}

		// look at the rods for indels or SNPs
		if (tracker != null) {
			for (VariantContext vc : tracker.getValues(known)) {
				switch (vc.getType()) {
				case INDEL:
					hasIndel = true;
					if (vc.isSimpleInsertion())
						hasInsertion = true;
					break;
				case SNP:
					hasPointEvent = true;
					break;
				case MIXED:
					hasPointEvent = true;
					hasIndel = true;
					if (vc.isSimpleInsertion())
						hasInsertion = true;
					break;
				default:
					break;
				}
				if (hasIndel)
					furthestStopPos = vc.getEnd();
			}
		}

		// look at the normal context to get deletions and positions with high
		// entropy
		if (context.hasBasePileup()) {
			ReadBackedPileup pileup = context.getBasePileup();

			int mismatchQualities = 0, totalQualities = 0;
			byte refBase = ref.getBase();
			for (PileupElement p : pileup) {
				if (BadMateFilter.hasBadMate(p.getRead()))
					continue;

				// check the ends of the reads to see how far they extend
				furthestStopPos = Math.max(furthestStopPos, p.getRead().getAlignmentEnd());

				// is it a deletion? (sanity check in case extended event missed
				// it)
				if (p.isDeletion()) {
					hasIndel = true;
				}

				// look for mismatches
				else {
					if (p.getBase() != refBase) {
						boolean secondPair = false;
						if (p.getRead().getReadPairedFlag())
							secondPair = p.getRead().getSecondOfPairFlag();
						if (BaseUtilsMore.isBisulfiteMismatch(refBase, p.getBase(), p.getRead().getReadNegativeStrandFlag(), secondPair)) {
							mismatchQualities += p.getQual();
						}
					}

					totalQualities += p.getQual();
				}
			}

			// make sure we're supposed to look for high entropy
			if (mismatchThreshold > 0.0 && mismatchThreshold <= 1.0 && pileup.getNumberOfElements() >= minReadsAtLocus && (double) mismatchQualities / (double) totalQualities >= mismatchThreshold)
				hasPointEvent = true;
		}

		// return null if no event occurred
		if (!hasIndel && !hasPointEvent)
			return null;

		// return null if we didn't find any usable reads/rods associated with
		// the event
		if (furthestStopPos == -1)
			return null;

		GenomeLoc eventLoc = context.getLocation();
		if (hasInsertion)
			eventLoc = getToolkit().getGenomeLocParser().createGenomeLoc(eventLoc.getContig(), eventLoc.getStart(), eventLoc.getStart() + 1);
		else if (hasIndel && !context.hasBasePileup())
			eventLoc = getToolkit().getGenomeLocParser().createGenomeLoc(eventLoc.getContig(), eventLoc.getStart(), furthestStopPos);

		EVENT_TYPE eventType = (hasIndel ? (hasPointEvent ? EVENT_TYPE.BOTH : EVENT_TYPE.INDEL_EVENT) : EVENT_TYPE.POINT_EVENT);

		return new Event(eventLoc, furthestStopPos, eventType);
	}

	@Override
	public void onTraversalDone(EventPair sum) {
		if (sum.left != null && sum.left.isReportableEvent())
			sum.intervals.add(sum.left.getLoc());
		if (sum.right != null && sum.right.isReportableEvent())
			sum.intervals.add(sum.right.getLoc());

		for (GenomeLoc loc : sum.intervals)
			out.println(loc);
	}

	@Override
	public EventPair reduce(Event value, EventPair sum) {
		if (value == null) {
			; // do nothing
		} else if (sum.left == null) {
			sum.left = value;
		} else if (sum.right == null) {
			if (canBeMerged(sum.left, value))
				sum.left = mergeEvents(sum.left, value);
			else
				sum.right = value;
		} else {
			if (canBeMerged(sum.right, value))
				sum.right = mergeEvents(sum.right, value);
			else {
				if (sum.right.isReportableEvent())
					sum.intervals.add(sum.right.getLoc());
				sum.right = value;
			}
		}

		return sum;
	}

	@Override
	public EventPair reduceInit() {
		return new EventPair(null, null);
	}

	@Override
	public EventPair treeReduce(EventPair lhs, EventPair rhs) {
		EventPair result;

		if (lhs.left == null) {
			result = rhs;
		} else if (rhs.left == null) {
			result = lhs;
		} else if (lhs.right == null) {
			if (rhs.right == null) {
				if (canBeMerged(lhs.left, rhs.left))
					result = new EventPair(mergeEvents(lhs.left, rhs.left), null, lhs.intervals, rhs.intervals);
				else
					result = new EventPair(lhs.left, rhs.left, lhs.intervals, rhs.intervals);
			} else {
				if (canBeMerged(lhs.left, rhs.left))
					result = new EventPair(mergeEvents(lhs.left, rhs.left), rhs.right, lhs.intervals, rhs.intervals);
				else {
					if (rhs.left.isReportableEvent())
						rhs.intervals.add(rhs.left.getLoc());
					result = new EventPair(lhs.left, rhs.right, lhs.intervals, rhs.intervals);
				}
			}
		} else if (rhs.right == null) {
			if (canBeMerged(lhs.right, rhs.left))
				result = new EventPair(lhs.left, mergeEvents(lhs.right, rhs.left), lhs.intervals, rhs.intervals);
			else {
				if (lhs.right.isReportableEvent())
					lhs.intervals.add(lhs.right.getLoc());
				result = new EventPair(lhs.left, rhs.left, lhs.intervals, rhs.intervals);
			}
		} else {
			if (canBeMerged(lhs.right, rhs.left)) {
				Event merge = mergeEvents(lhs.right, rhs.left);
				if (merge.isReportableEvent())
					lhs.intervals.add(merge.getLoc());
			} else {
				if (lhs.right.isReportableEvent())
					lhs.intervals.add(lhs.right.getLoc());
				if (rhs.left.isReportableEvent())
					rhs.intervals.add(rhs.left.getLoc());
			}

			result = new EventPair(lhs.left, rhs.right, lhs.intervals, rhs.intervals);
		}

		return result;
	}

	private enum EVENT_TYPE {
		BOTH, INDEL_EVENT, POINT_EVENT
	}

	class Event {
		public int furthestStopPos;

		private int eventStartPos;
		private int eventStopPos;
		private GenomeLoc loc;
		private ArrayList<Integer> pointEvents = new ArrayList<Integer>();
		private EVENT_TYPE type;

		public Event(GenomeLoc loc, int furthestStopPos, EVENT_TYPE type) {
			this.loc = loc;
			this.furthestStopPos = furthestStopPos;
			this.type = type;

			if (type == EVENT_TYPE.INDEL_EVENT || type == EVENT_TYPE.BOTH) {
				eventStartPos = loc.getStart();
				eventStopPos = loc.getStop();
			} else {
				eventStartPos = -1;
				eventStopPos = -1;
			}

			if (type == EVENT_TYPE.POINT_EVENT || type == EVENT_TYPE.BOTH) {
				pointEvents.add(loc.getStart());
			}
		}

		public GenomeLoc getLoc() {
			return getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(), eventStartPos, eventStopPos);
		}

		public boolean isReportableEvent() {
			return getToolkit().getGenomeLocParser().isValidGenomeLoc(loc.getContig(), eventStartPos, eventStopPos, true) && eventStopPos >= 0 && eventStopPos - eventStartPos < maxIntervalSize;
		}

		public void merge(Event e) {

			// merges only get called for events with certain types
			if (e.type == EVENT_TYPE.INDEL_EVENT || e.type == EVENT_TYPE.BOTH) {
				if (eventStartPos == -1)
					eventStartPos = e.eventStartPos;
				eventStopPos = e.eventStopPos;
				furthestStopPos = e.furthestStopPos;
			}

			if (e.type == EVENT_TYPE.POINT_EVENT || e.type == EVENT_TYPE.BOTH) {
				int newPosition = e.pointEvents.get(0);
				if (pointEvents.size() > 0) {
					int lastPosition = pointEvents.get(pointEvents.size() - 1);
					if (newPosition - lastPosition < windowSize) {
						eventStopPos = Math.max(eventStopPos, newPosition);
						furthestStopPos = e.furthestStopPos;

						if (eventStartPos == -1)
							eventStartPos = lastPosition;
						else
							eventStartPos = Math.min(eventStartPos, lastPosition);
					} else if (eventStartPos == -1 && e.eventStartPos != -1) {
						eventStartPos = e.eventStartPos;
						eventStopPos = e.eventStopPos;
						furthestStopPos = e.furthestStopPos;
					}
				}
				pointEvents.add(newPosition);
			}
		}
	}

	class EventPair {
		public TreeSet<GenomeLoc> intervals = new TreeSet<GenomeLoc>();
		public Event left, right;

		public EventPair(Event left, Event right) {
			this.left = left;
			this.right = right;
		}

		public EventPair(Event left, Event right, TreeSet<GenomeLoc> set1, TreeSet<GenomeLoc> set2) {
			this.left = left;
			this.right = right;
			intervals.addAll(set1);
			intervals.addAll(set2);
		}
	}

}
