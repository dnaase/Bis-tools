/**
 * 
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp.filters;

import main.java.edu.usc.epigenome.uecgatk.bissnp.BaseUtilsMore;
import main.java.edu.usc.epigenome.uecgatk.bissnp.BisSNPUtils;
import main.java.edu.usc.epigenome.uecgatk.bissnp.BisulfiteSAMConstants;
import htsjdk.samtools.SAMRecord;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.engine.filters.ReadFilter;
import org.broadinstitute.gatk.utils.BaseUtils;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Sep 15, 2013 10:49:14 PM
 * 
 */
public class BisulfiteMismatchReadsFilter extends ReadFilter {

	@Argument(fullName = "max_mismatches", shortName = "mm", doc = "Maximum percentage of non-bisulfite mismatches within a read for a read to be used for calling. Default: 0.1(10% of mismatches allowed)", required = false)
	public static double MAX_MISMATCHES = 0.1;
	
	/* (non-Javadoc)
	 * @see net.sf.picard.filter.SamRecordFilter#filterOut(htsjdk.samtools.SAMRecord)
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
		
		boolean negativeStrand = read.getReadNegativeStrandFlag();
		boolean secondEnd = read.getReadPairedFlag() && read.getSecondOfPairFlag();

		byte[] refBases = BaseUtilsMore.toUpperCase(BisSNPUtils.modifyRefSeqByCigar(BisSNPUtils.refStrFromMd(read), read.getCigarString()));
		
		
		byte[] bases = BaseUtilsMore.toUpperCase(BisSNPUtils.getClippedReadsBase(read));
		if (negativeStrand) {
			bases = BaseUtils.simpleReverseComplement(bases);
			refBases = BaseUtils.simpleReverseComplement(refBases);

		}
		
		int numberOfMismatches = 0;
		for(int i = 0; i < bases.length; i++){
			if( !BaseUtils.basesAreEqual(refBases[i], bases[i]) && BaseUtilsMore.isBisulfiteMismatch(refBases[i], bases[i],negativeStrand, secondEnd))
				numberOfMismatches++;
			if(numberOfMismatches > MAX_MISMATCHES * bases.length){
				//System.err.println(read.getReadString() + "\t" + read.getCigarString());
			//	System.err.println(new String(refBases) + "\t" + refBases.length);
			//	System.err.println(new String(bases) + "\t" + bases.length + "\t" + negativeStrand + "\t" + secondEnd + "\t" + numberOfMismatches);

				return true;
			}
				
		}
		//System.err.println(numberOfMismatches);
		return false;
	}

	@Override
	public boolean filterOut(SAMRecord arg0, SAMRecord arg1) {
		// TODO Auto-generated method stub
		return false;
	}

}
