/**
 * get rid of Not proper paired reads
 */
package edu.usc.epigenome.uecgatk.bissnp.filters;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.ReadFilter;


/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Dec 11, 2012 11:14:32 AM
 * 
 */
public class NotProperPairedReadFilter extends ReadFilter {
	@Argument(fullName = "use_badly_mated_reads", shortName = "badMate", doc = "use badly mated reads for calling", required = false)
	public boolean USE_BADLY_MATED_READS = false;
	
	public boolean filterOut( final SAMRecord read ) {
        return read.getReadPairedFlag() && !USE_BADLY_MATED_READS && BadMateFilter.hasBadMate(read) && !read.getProperPairFlag();
    }


}
