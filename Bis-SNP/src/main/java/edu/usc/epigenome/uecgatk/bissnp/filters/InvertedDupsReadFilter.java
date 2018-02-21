/**
 * 
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp.filters;

import htsjdk.samtools.SAMRecord;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.engine.filters.ReadFilter;

import main.java.edu.usc.epigenome.uecgatk.bissnp.BisulfiteEnums.INVERT_DUPS;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Dec 11, 2012 11:22:51 AM
 * 
 */
public class InvertedDupsReadFilter extends ReadFilter {

	@Argument(fullName = "invert_dups_usage", shortName = "invDups", doc = "how to use the end when met inverted dups reads[USE_ONLY_1ST_END, USE_BOTH_END, NOT_TO_USE]. (Default: USE_ONLY_1ST_END)", required = false)
    public INVERT_DUPS invDups = INVERT_DUPS.USE_ONLY_1ST_END;
	
	/* (non-Javadoc)
	 * @see net.sf.picard.filter.SamRecordFilter#filterOut(htsjdk.samtools.SAMRecord)
	 */
	@Override
	public boolean filterOut(SAMRecord samRecord) {
		if (samRecord.getReadPairedFlag() && (samRecord.getAlignmentStart() == samRecord.getMateAlignmentStart() && samRecord.getReadNegativeStrandFlag() == samRecord.getMateNegativeStrandFlag()))
		{
			if(invDups == INVERT_DUPS.USE_ONLY_1ST_END){
				if (samRecord.getSecondOfPairFlag()) 
					return true;
			}
			else if(invDups == INVERT_DUPS.NOT_TO_USE){
				return true;
			}
		}
		return false;
	}


}
