/**
 * 
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp.filters;

import main.java.edu.usc.epigenome.uecgatk.bissnp.BaseUtilsMore;
import main.java.edu.usc.epigenome.uecgatk.bissnp.BisSNPUtils;
import main.java.edu.usc.epigenome.uecgatk.bissnp.BisulfiteSAMConstants;
import htsjdk.samtools.SAMRecord;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.filters.ReadFilter;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;


/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Feb 8, 2013 4:48:00 PM
 * 
 */
public class BisulfiteIncompleteConvReadsFilter extends ReadFilter {

	@Argument(fullName = "minmum_pattern_converted", shortName = "minPatConv", doc = "minimum frequency of C pattern (e.g.HCH) allowed to get methylated inside the reads (methylated C/ reference C in the reads). Default: 0.4 ", required = false)
	public double minPatConv = 0.4;
	
	@Argument(fullName = "conversion_pattern_to_check", shortName = "patConv", doc = "define the methylation pattern to check for bisulfite_conversion. Default: HCH", required = false)
    public String patConv = "CH";
	
	public BisulfiteIncompleteConvReadsFilter(){
		
	}
	
	@Override
	public boolean filterOut(SAMRecord read) {
		if(read.getStringAttribute(BisulfiteSAMConstants.MD_TAG) == null){
			return false;
		}
		try {
			return incompleteConversion(read);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return false;
	}
	
	private boolean incompleteConversion(SAMRecord read) throws Exception{
		
		boolean negativeStrand = read.getReadNegativeStrandFlag();
		
		byte[] refBases = BaseUtilsMore.toUpperCase(BisSNPUtils.modifyRefSeqByCigar(BisSNPUtils.refStrFromMd(read), read.getCigarString()));
		
		
		byte[] bases = BaseUtilsMore.toUpperCase(BisSNPUtils.getClippedReadsBase(read));
		//if(negativeStrand){
		//	bases = BisSNPUtils.complementArray(bases);
		//}
		
		if(read.getReadPairedFlag() && read.getSecondOfPairFlag())
			negativeStrand = !negativeStrand;
		
		String pattern = patConv;
		
		int numberOfPatternInRef = read.getReadUnmappedFlag() ? -1 : 0;

		
		int readLength = bases.length;

		if(negativeStrand){
			bases = BaseUtils.simpleReverseComplement(bases);

			if(numberOfPatternInRef != -1){
				refBases = BaseUtils.simpleReverseComplement(refBases);
			}
			
		}
		
		if(readLength < pattern.length())
			return true;
		int numberOfPattern = 0;
		
		byte[] patterns = pattern.getBytes();
		
		
		for(int i = 0; i <= readLength-patterns.length; i++){
			short numMatchesInRef=0;
			short numMatchesInReads=0;
			for(int j = i, index = 0; index < patterns.length; index++, j++){
				
				if(numberOfPatternInRef != -1){ //mean reads is mapped
					
					if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(patterns[index], refBases[j])){
						numMatchesInRef++;
						
					}
					
				}
				if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(patterns[index], bases[j])){
					numMatchesInReads++;
				}
				
			}
			if(numMatchesInReads == patterns.length){
				numberOfPattern++;
			}
			if(numberOfPatternInRef != -1 && numMatchesInRef == patterns.length){
				numberOfPatternInRef++;
			}
		}
		if(((numberOfPatternInRef == -1 || numberOfPatternInRef == 0) && numberOfPattern == 0) || (numberOfPatternInRef != -1 && (double)numberOfPattern/(double)numberOfPatternInRef < minPatConv)){
			return false;
		}else{
			//System.err.println(new String(refBasesInWidnow) + "\t" + refBasesInWidnow.length);
			//		System.err.println(new String(refBases) + "\t" + len);
			//		System.err.println(new String(bases) + "\t" + readLength + "\t" + negativeStrand + "\t" + read.getSecondOfPairFlag() + "\t" + numberOfPattern + "\t" + numberOfPatternInRef);
			//		System.err.println(new String(baseQs) + "\t" + read.getReadName());
			return true;
		}
		
	}

}
