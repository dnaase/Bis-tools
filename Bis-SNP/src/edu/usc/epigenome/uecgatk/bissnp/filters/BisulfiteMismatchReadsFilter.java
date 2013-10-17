/**
 * 
 */
package edu.usc.epigenome.uecgatk.bissnp.filters;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.utils.BaseUtils;

import edu.usc.epigenome.uecgatk.bissnp.BaseUtilsMore;
import edu.usc.epigenome.uecgatk.bissnp.BisSNPUtils;
import edu.usc.epigenome.uecgatk.bissnp.BisulfiteSAMConstants;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Sep 15, 2013 10:49:14 PM
 * 
 */
public class BisulfiteMismatchReadsFilter extends ReadFilter {

	@Argument(fullName = "max_mismatches", shortName = "mm", doc = "Maximum percentage of non-bisulfite mismatches within a read for a read to be used for calling. Default: 0.05(5% of mismatches allowed)", required = false)
	public static double MAX_MISMATCHES = 0.05;
	
	/* (non-Javadoc)
	 * @see net.sf.picard.filter.SamRecordFilter#filterOut(net.sf.samtools.SAMRecord)
	 */
	@Override
	public boolean filterOut(SAMRecord read) {
		if(read.getStringAttribute(BisulfiteSAMConstants.MD_TAG) == null){
			return false;
		}
		try {
			return hasTooManyBisulfiteMismathces(read);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return false;
	}
	
	public static boolean hasTooManyBisulfiteMismathces(SAMRecord read) throws Exception{
		int readLength = read.getReadLength();
		boolean negativeStrand = read.getReadNegativeStrandFlag();
		boolean secondEnd = read.getReadPairedFlag() && read.getSecondOfPairFlag();

		byte[] bases = read.getReadBases();
		
		//TODO: need to check window length in boundary...
		byte[] refBases = BisSNPUtils.refStrFromMd(read);
		//System.err.println(new String(refBasesInWidnow) + "\t" + refBasesInWidnow.length);
		//System.err.println(new String(refBases) + "\t" + len);
		//System.err.println(new String(bases) + "\t" + readLength + "\t" + negativeStrand + "\t" + secondEnd);
		int numberOfMismatches = 0;
		for(int i = 0; i < readLength; i++){
			if( !BaseUtils.basesAreEqual(refBases[i], bases[i]) && BaseUtilsMore.isBisulfiteMismatch(refBases[i], bases[i],negativeStrand, secondEnd))
				numberOfMismatches++;
			if(numberOfMismatches > MAX_MISMATCHES * readLength)
				return true;
		}
		//System.err.println(numberOfMismatches);
		return false;
	}

}
