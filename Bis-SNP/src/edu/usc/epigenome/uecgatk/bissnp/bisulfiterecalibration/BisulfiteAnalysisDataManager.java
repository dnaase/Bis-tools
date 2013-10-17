/**
 * 
 */
package edu.usc.epigenome.uecgatk.bissnp.bisulfiterecalibration;

import java.util.ArrayList;

import org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum;
import org.broadinstitute.sting.utils.collections.NestedHashMap;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Apr 13, 2012 8:01:29 PM
 * 
 */
public class BisulfiteAnalysisDataManager {
	private ArrayList<NestedHashMap> dataCollapsedByCovariate;
	private NestedHashMap dataCollapsedReadGroup;

	BisulfiteAnalysisDataManager() {
	}

	BisulfiteAnalysisDataManager(final int numCovariates) {
		dataCollapsedReadGroup = new NestedHashMap();
		dataCollapsedByCovariate = new ArrayList<NestedHashMap>();
		for (int iii = 0; iii < numCovariates - 1; iii++) {
			dataCollapsedByCovariate.add(new NestedHashMap());
		}
	}

	/**
	 * Add the given mapping to all of the collapsed hash tables
	 * 
	 * @param key
	 *            The list of comparables that is the key for this mapping
	 * @param fullDatum
	 *            The RecalDatum which is the data for this mapping
	 * @param IGNORE_QSCORES_LESS_THAN
	 *            The threshold in report quality for adding to the aggregate collapsed table
	 */
	public final void addToAllTables(final Object[] key, final RecalDatum fullDatum, final int IGNORE_QSCORES_LESS_THAN) {

		int qscore = Integer.parseInt(key[1].toString());
		RecalDatum collapsedDatum;
		final Object[] readGroupCollapsedKey = new Object[1];
		final Object[] covariateCollapsedKey = new Object[2];

		if (!(qscore < IGNORE_QSCORES_LESS_THAN)) {
			// Create dataCollapsedReadGroup, the table where everything except
			// read group has been collapsed
			readGroupCollapsedKey[0] = key[0]; // Make a new key with just the
												// read group
			collapsedDatum = (RecalDatum) dataCollapsedReadGroup.get(readGroupCollapsedKey);
			if (collapsedDatum == null) {
				dataCollapsedReadGroup.put(new RecalDatum(fullDatum), readGroupCollapsedKey);
			} else {
				collapsedDatum.combine(fullDatum);
			}
		}

		// Create dataCollapsedByCovariate's, the tables where everything except
		// read group and given covariate has been collapsed
		for (int iii = 0; iii < dataCollapsedByCovariate.size(); iii++) {
			if (iii == 0 || !(qscore < IGNORE_QSCORES_LESS_THAN)) {
				covariateCollapsedKey[0] = key[0]; // Make a new key with the
													// read group ...
				Object theCovariateElement = key[iii + 1]; // and the given
															// covariate
				if (theCovariateElement != null) {
					covariateCollapsedKey[1] = theCovariateElement;
					collapsedDatum = (RecalDatum) dataCollapsedByCovariate.get(iii).get(covariateCollapsedKey);
					if (collapsedDatum == null) {
						dataCollapsedByCovariate.get(iii).put(new RecalDatum(fullDatum), covariateCollapsedKey);
					} else {
						collapsedDatum.combine(fullDatum);
					}
				}
			}
		}
	}

	/**
	 * Get the appropriate collapsed table out of the set of all the tables held by this Object
	 * 
	 * @param covariate
	 *            Which covariate indexes the desired collapsed HashMap
	 * @return The desired collapsed HashMap
	 */
	public final NestedHashMap getCollapsedTable(final int covariate) {
		if (covariate == 0) {
			return dataCollapsedReadGroup; // Table where everything except read
											// group has been collapsed
		} else {
			return dataCollapsedByCovariate.get(covariate - 1);
		}
	}

}
