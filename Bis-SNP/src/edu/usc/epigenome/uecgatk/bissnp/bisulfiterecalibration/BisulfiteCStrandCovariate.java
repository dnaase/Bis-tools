/**
 * measure possible strand bias of bisulfite-seq, put read into bisulfiteC strand as "C" by C or T's
 * base number in the reads is bigger than expectation. Expectation=readLen*0.25
 */
package edu.usc.epigenome.uecgatk.bissnp.bisulfiterecalibration;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Apr 22, 2012 8:54:13 PM
 * 
 */
public class BisulfiteCStrandCovariate implements StandardCovariate {

	@Override
	public Comparable getValue(String str) {
		// TODO Auto-generated method stub
		return str;
	}

	@Override
	public void getValues(GATKSAMRecord read, Comparable[] comparable) {
		// TODO Auto-generated method stub

	}

	public void getValues(final GATKSAMRecord read, final Comparable[] comparable, ReferenceContext ref) {
		// TODO Auto-generated method stub
		final String BisulfiteStrandId;
		final boolean negativeStrand = read.getReadNegativeStrandFlag();
		final int readLength = read.getReadLength();
		byte[] bases = read.getReadBases();
		byte[] refBases = new byte[readLength];
		boolean secondPairFlag = false;
		if (read.getReadPairedFlag()) {
			secondPairFlag = read.getSecondOfPairFlag();
		}
		int start = read.getAlignmentStart();

		int end = read.getUnclippedEnd();
		GenomeLoc locBoundary = ref.getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(), start, end);

		for (int i = 0, cor = start; i < readLength; i++) {
			GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(), cor);
			if (!ref.getWindow().overlapsP(locBoundary))
				return;

			ReferenceContext tmpRef = new ReferenceContext(ref.getGenomeLocParser(), loc, ref.getWindow(), ref.getBases());

			refBases[i] = tmpRef.getBase();

			cor++;

		}

		if (negativeStrand) {
			bases = BaseUtils.simpleReverseComplement(bases);
			refBases = BaseUtils.simpleReverseComplement(refBases);

		}
		int numC = 0;

		for (int i = 0; i < readLength; i++) {
			if (secondPairFlag) {
				if (negativeStrand) {
					if (refBases[i] == BaseUtils.G && bases[i] == BaseUtils.G) {
						numC++;
					}
				} else {
					if (refBases[i] == BaseUtils.C && bases[i] == BaseUtils.C) {
						numC++;
					}
				}

			} else {
				if (negativeStrand) {
					if (refBases[i] == BaseUtils.G) {
						if (bases[i] == BaseUtils.G || bases[i] == BaseUtils.A) {
							numC++;
						}
					}
				} else {
					if (refBases[i] == BaseUtils.C) {
						if (bases[i] == BaseUtils.C || bases[i] == BaseUtils.T) {
							numC++;
						}
					}
				}

			}

		}
		if (numC > read.getReadLength() * 0.25) {
			BisulfiteStrandId = "C";
		} else {
			BisulfiteStrandId = "G";
		}
		for (int i = 0; i < read.getReadLength(); i++) {
			comparable[i] = BisulfiteStrandId;
		}
	}

	@Override
	public void initialize(RecalibrationArgumentCollection RAC) {
		// TODO Auto-generated method stub

	}

}
