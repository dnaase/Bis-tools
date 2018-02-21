/**
 * 
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp.filters;

import htsjdk.samtools.SAMRecord;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.engine.filters.ReadFilter;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Sep 15, 2013 10:54:38 PM
 * 
 */
public class MappingQualityFilter extends ReadFilter {

	@Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum mapping quality required to consider a base for calling", required = false)
	public int MIN_MAPPING_QUALTY_SCORE = 30;
	 
	/* (non-Javadoc)
	 * @see net.sf.picard.filter.SamRecordFilter#filterOut(htsjdk.samtools.SAMRecord)
	 */
	@Override
	public boolean filterOut(SAMRecord arg0) {
		// TODO Auto-generated method stub
		return false;
	}

}
