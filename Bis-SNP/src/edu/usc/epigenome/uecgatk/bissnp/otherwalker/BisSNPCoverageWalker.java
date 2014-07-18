/**
 * 
 */
package edu.usc.epigenome.uecgatk.bissnp.otherwalker;

import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.coverage.DepthOfCoverageWalker;

import edu.usc.epigenome.uecgatk.bissnp.filters.BisulfiteFivePrimeConvReadsFilter;
import edu.usc.epigenome.uecgatk.bissnp.filters.BisulfiteIncompleteConvReadsFilter;
import edu.usc.epigenome.uecgatk.bissnp.filters.BisulfiteMismatchReadsFilter;
import edu.usc.epigenome.uecgatk.bissnp.filters.InvertedDupsReadFilter;
import edu.usc.epigenome.uecgatk.bissnp.filters.MappingQualityFilter;
import edu.usc.epigenome.uecgatk.bissnp.filters.NotProperPairedReadFilter;

/**
 * Add filters into DepthOfCoverageWalker to control the estimation in different bad reads condition
 * 
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Oct 30, 2013 5:57:16 PM
 * 
 */
@ReadFilters({UnmappedReadFilter.class, DuplicateReadFilter.class, NotPrimaryAlignmentFilter.class, FailsVendorQualityCheckFilter.class,  MappingQualityZeroFilter.class, NotProperPairedReadFilter.class, InvertedDupsReadFilter.class, BisulfiteMismatchReadsFilter.class, BisulfiteIncompleteConvReadsFilter.class})
@By(DataSource.READS)
public class BisSNPCoverageWalker extends DepthOfCoverageWalker{

	/**
	 * 
	 */
	

}
