package edu.usc.epigenome.uecgatk.bissnp.bisulfiterecalibration;

import org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatumOptimized;
import org.broadinstitute.sting.utils.BaseUtils;

import edu.usc.epigenome.uecgatk.bissnp.BaseUtilsMore;
import edu.usc.epigenome.uecgatk.bissnp.BisSNPUtils;

public class BisulfiteRecalDatumOptimized extends RecalDatumOptimized {

	public BisulfiteRecalDatumOptimized() {
		// TODO Auto-generated constructor stub
	}

	public BisulfiteRecalDatumOptimized(long numObservations, long numMismatches) {
		super(numObservations, numMismatches);
		// TODO Auto-generated constructor stub
	}

	public BisulfiteRecalDatumOptimized(RecalDatumOptimized copy) {
		super(copy);
		// TODO Auto-generated constructor stub
	}

	public synchronized void incrementBaseCountsBisulfite(final byte curBase, final byte refBase, boolean negStrand, boolean secondEnd) {

		if (BaseUtils.simpleBaseToBaseIndex(curBase) != BaseUtils.simpleBaseToBaseIndex(refBase)) {
			if (BaseUtilsMore.isBisulfiteMismatch(refBase,curBase, negStrand,secondEnd)) {
				increment(1, 1);
			} else {
				increment(1, 0);
			}
		} else {

			increment(1, 0);
		}

	//	increment(1, (BaseUtilsMore.isBisulfiteMismatch(refBase,curBase, negStrand,secondEnd)) ? 1 : 0);
	}
}
