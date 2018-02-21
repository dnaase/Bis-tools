/**
 * 
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp.otherwalker;

import org.broadinstitute.gatk.engine.filters.DuplicateReadFilter;
import org.broadinstitute.gatk.engine.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.gatk.engine.filters.MappingQualityZeroFilter;
import org.broadinstitute.gatk.engine.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.gatk.engine.filters.UnmappedReadFilter;
import org.broadinstitute.gatk.engine.walkers.By;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.ReadFilters;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.tools.walkers.coverage.DepthOfCoverage;

import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.BisulfiteIncompleteConvReadsFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.BisulfiteMismatchReadsFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.InvertedDupsReadFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.NotProperPairedReadFilter;

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
public class BisSNPCoverageWalker extends DepthOfCoverage{

	/**
	 * 
	 */
	

}
