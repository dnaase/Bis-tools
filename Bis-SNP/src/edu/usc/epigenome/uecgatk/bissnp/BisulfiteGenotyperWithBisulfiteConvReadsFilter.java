/**
 * 
 */
package edu.usc.epigenome.uecgatk.bissnp;

import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;

import edu.usc.epigenome.uecgatk.bissnp.filters.BisulfiteIncompleteConvReadsFilter;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Feb 8, 2013 4:45:49 PM
 * 
 */
@ReadFilters( {BisulfiteIncompleteConvReadsFilter.class} )
public class BisulfiteGenotyperWithBisulfiteConvReadsFilter extends BisulfiteGenotyper {

	/**
	 * 
	 */
	public BisulfiteGenotyperWithBisulfiteConvReadsFilter() {
		// TODO Auto-generated constructor stub
	}

}
