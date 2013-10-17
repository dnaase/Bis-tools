package edu.usc.epigenome.uecgatk.bissnp;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;

import net.sf.samtools.SAMFileWriter;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Mar 20, 2012 3:39:46 PM
 * 
 */
public class DownsamplingBAM {

	protected SAMFileWriter samWriter = null;
	private BisulfiteArgumentCollection BAC;
	private int SAMPLE_READS_MEAN_COVERAGE = 30;

	/**
	 * 
	 */
	public DownsamplingBAM(BisulfiteArgumentCollection BAC, SAMFileWriter samWriter) {
		// TODO Auto-generated constructor stub
		this.BAC = BAC;
		this.samWriter = samWriter;
	}
	
	public DownsamplingBAM(BisulfiteArgumentCollection BAC, SAMFileWriter samWriter, int coverage) {
		// TODO Auto-generated constructor stub
		this.BAC = BAC;
		this.samWriter = samWriter;
		this.SAMPLE_READS_MEAN_COVERAGE = coverage;
	}

	/*
	 * we use different downsampling strategy. e.g. when downsampling 10X, it randomly pick up
	 * s*10/r reads.(r is the mean coverage of the sample, we use 30 here for our sample, s is the
	 * total reads covered in this position) it does NOT treat reads like GATK, when reads number
	 * more than 10, cut off the reads number to 10, and keep the same when reads number lower than
	 * 10.
	 */
	public int downsamplingBamFile(AlignmentContext rawContext) {
		int downsampledCoverage = 0;
		if (rawContext.hasReads()) {
			String tag = "Xi";
			Integer coverageMarked = 0;

			int covergaeLimit = BAC.orcad;
			//covergaeLimit = Math.max((covergaeLimit * rawContext.getBasePileup().depthOfCoverage()) / SAMPLE_READS_MEAN_COVERAGE, 1);
			covergaeLimit = (int)((covergaeLimit * rawContext.getBasePileup().depthOfCoverage()) / SAMPLE_READS_MEAN_COVERAGE);
			
			ReadBackedPileup downsampledPileup = getDownsampledPileup(rawContext.getBasePileup(), covergaeLimit);
			downsampledCoverage = downsampledPileup.depthOfCoverage();
			for (PileupElement p : rawContext.getBasePileup()) {
				if (p.getRead().getIntegerAttribute(tag) != null) {
					if (p.getRead().getIntegerAttribute(tag) == 2)

						if (p.getRead().getIntegerAttribute(tag) == 1)
							coverageMarked++;
				}

			}

			for (PileupElement p : downsampledPileup) {

				if (p.getRead().getIntegerAttribute(tag) != null) {

				}
				if (coverageMarked >= covergaeLimit)
					break;
				if (p.getRead().getIntegerAttribute(tag) == null) {
					samWriter.addAlignment(p.getRead());
					p.getRead().setAttribute(tag, 1);
					coverageMarked++;
				}

			}
			for (PileupElement p : rawContext.getBasePileup()) {
				if (p.getRead().getIntegerAttribute(tag) == null)
					p.getRead().setAttribute(tag, 2);
			}

		}
		return downsampledCoverage;
	}

	public void downsamplingBamFileLikeGATK(AlignmentContext rawContext, GenomeAnalysisEngine toolkit) {
		if (rawContext.hasReads()) {
			String tag = "Xi";
			Integer coverageMarked = 0;
			int covergaeLimit = toolkit.getArguments().downsampleCoverage;
			ReadBackedPileup downsampledPileup = rawContext.getBasePileup().getDownsampledPileup(covergaeLimit);
			for (PileupElement p : rawContext.getBasePileup()) {
				if (p.getRead().getIntegerAttribute(tag) != null) {
					if (p.getRead().getIntegerAttribute(tag) == 2)
						System.out.println("loc: " + rawContext.getLocation().getStart() + " tag: " + p.getRead().getIntegerAttribute(tag));
					if (p.getRead().getIntegerAttribute(tag) == 1)
						coverageMarked++;
				}

			}

			for (PileupElement p : downsampledPileup) {

				if (coverageMarked >= covergaeLimit)
					break;
				if (p.getRead().getIntegerAttribute(tag) == null) {
					samWriter.addAlignment(p.getRead());
					p.getRead().setAttribute(tag, 1);
					coverageMarked++;
				}

			}
			for (PileupElement p : rawContext.getBasePileup()) {
				if (p.getRead().getIntegerAttribute(tag) == null)
					p.getRead().setAttribute(tag, 2);
			}

		}

	}

	public ReadBackedPileup getDownsampledPileup(ReadBackedPileup pileup, int desiredCov) {
		if (pileup.depthOfCoverage() <= desiredCov)
			return pileup;

		// randomly choose numbers corresponding to positions in the reads list
		Random generator = new Random();
		TreeSet<Integer> positions = new TreeSet<Integer>();
		for (int i = 0; i < desiredCov; /* no update */) {
			if (positions.add(generator.nextInt(pileup.depthOfCoverage())))
				i++;
		}
		GenomeLoc loc = pileup.getLocation();
		List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
		;
		List<Integer> elementOffsets = new ArrayList<Integer>();
		int i = 0;
		for (PileupElement p : pileup) {
			if (positions.contains(i)) {
				int elementOffset = p.getOffset();
				if (elementOffset < 0 || elementOffset > p.getRead().getReadLength() - 1)
					continue;
				elementOffsets.add(elementOffset);
				reads.add(p.getRead());
			}
			i++;
		}
		ReadBackedPileup downsampledPileup = new ReadBackedPileupImpl(loc, reads, elementOffsets);

		return downsampledPileup;

	}

}
