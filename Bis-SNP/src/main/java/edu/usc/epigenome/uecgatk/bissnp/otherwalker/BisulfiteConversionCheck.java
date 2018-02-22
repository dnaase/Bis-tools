/**
 * This walker is used for check percentage of CH/HCH get methylated in genome wide to estimate the incomplete bisulfite conversion reads frequency
 * It will output the histogram of # of methylated CH/HCH in a reads. and also the filtered BAM file. 
 * This may be later:  Also a plot of methylation level along cycle CH/HCH
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp.otherwalker;

import java.io.File;
import java.io.PrintStream;


import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterFactory;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.walkers.By;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.ReadFilters;
import org.broadinstitute.gatk.engine.walkers.ReadWalker;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.Requires;
import org.broadinstitute.gatk.engine.walkers.WalkerName;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.engine.filters.*;

import main.java.edu.usc.epigenome.uecgatk.bissnp.BaseUtilsMore;


/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Feb 7, 2013 4:53:14 PM
 * 
 */
@WalkerName("BisulfiteConversionCheck")
@Reference(window = @Window(start = -500, stop = 500))
@Requires({DataSource.READS, DataSource.REFERENCE})
@By(DataSource.READS)
@ReadFilters( {FailsVendorQualityCheckFilter.class} )
public class BisulfiteConversionCheck extends ReadWalker<SAMRecord, Long> {

	@Output(fullName = "histogram_methylated_pattern_in_reads", shortName = "patternHist", doc = "output histogram of number of methylated pattern(HCH) in the reads ", required = false)
	public PrintStream patternHistWriter = null;
	
	@Output(fullName = "cleanBam", shortName = "cleanBam", doc = "output of BAM file with filtered incomplete bisulfite conversion reads(also filter out Non-unique aligned reads, PCR duplicate reads)", required = false)
	public String cleanBam = null;
	
	@Output(fullName = "dirtyBam", shortName = "dirtyBam", doc = "output of BAM file contain reads that been filtered out(not include Non-unique aligned reads, PCR duplicate reads, and some -filterInvertedDupsReads filters)", required = false)
	public String dirtyBam = null;
	
	@Output(fullName = "methylated_pattern_along_cycle", shortName = "patternCycle", doc = "output pattern(HCH) methylation pattern along cycle in the reads, -100 to -1 is for 2nd end in PE ", required = false)
	public PrintStream patternMethyAlongCycleWriter = null;
	
	@Argument(fullName = "methylation_pattern_to_check", shortName = "pattern", doc = "define the methylation pattern to check for bisulfite_conversion. Default: HCH", required = false)
    public String pattern = "HCH";
	
	@Argument(fullName = "minimum_base_quality", shortName = "minBaseQ", doc = "minimum base quality score required to check. Default: 5", required = false)
    public int minBaseQ = 5;
	
	@Argument(fullName = "minimum_number_C_pattern_methylated", shortName = "minNumMethyC", doc = "minimum number of C pattern allowed to get methylated inside the reads (probably use percentage of methylated C in the reads??). Default: 2", required = false)
    public int minNumMethyC = 2;
	
	@Argument(fullName = "only_mapped_reads", shortName = "onlyMappedReads", doc = "filter out unmapped reads before doing stats. Default: false", required = false)
    public boolean onlyMappedReads = false;
	
	@Argument(fullName = "only_proper_paired_reads", shortName = "onlyProperPairedReads", doc = "filter out non-proper paired reads before doing stats. Default: false", required = false)
    public boolean onlyProperPairedReads = false;
	
	@Argument(fullName = "filter_invert_dups_reads", shortName = "filterInvertedDupsReads", doc = "filter out inverted dups reads before doing stats. Default: false", required = false)
    public boolean filterInvertedDupsReads = false;
	
	@Argument(fullName = "filter_invert_dups_reads_complete", shortName = "filterInvertedDupsReadsComplete", doc = "filter out inverted dups reads completedly (even 1st end is not left)before doing stats. Default: false", required = false)
    public boolean filterInvertedDupsReadsComplete = false;
	
//	@Argument(fullName = "filter_overlapped_fragment", shortName = "filterOverlappedFragment", doc = "filter out filter Overlapped Fragment when two ends overlapped. Default: false", required = false)
  //  public boolean filterOverlappedFragment = false;
	
	@Argument(fullName = "filter_pcr_duplicate_reads", shortName = "filterPcrDupsReads", doc = "filter out PCR dups reads before doing stats. Default: true", required = false)
    public boolean filterPcrDupsReads = true;
	
	@Argument(fullName = "filter_no_primary_aligned_reads", shortName = "filterNoPrimaryAlignedReads", doc = "filter out no primary aligned reads before doing stats. Default: true", required = false)
    public boolean filterNoPrimaryAlignedReads = true;
	
	@Argument(fullName = "minimum_freq_C_pattern_methylated", shortName = "minFreqMethyC", doc = "minimum frequency of C pattern allowed to get methylated inside the reads (methylated C/ reference C in the reads). Default: 0.4", required = false)
    public double minFreqMethyC = 0.4;
	
	
	//private HashMap<Integer, Pair<Long,Long>> patternMethyCycleMap = null;
	private SAMFileWriter bamWriter = null;
	private SAMFileWriter bamDirtyWriter = null;
	
	public void initialize(){
		//Pair<Long,Long> initialNum = new Pair<Long,Long>(0L,0L);
		//patternMethyCycleMap = new HashMap<Integer, Pair<Long,Long>>();
		if(cleanBam != null){
			bamWriter = new SAMFileWriterFactory().makeBAMWriter(getToolkit().getSAMFileHeader(),false,new File(cleanBam));
		}
		if(dirtyBam != null){
			bamDirtyWriter = new SAMFileWriterFactory().makeBAMWriter(getToolkit().getSAMFileHeader(),false,new File(dirtyBam));
		}
		
	}
	
	/* (non-Javadoc)
	 * @see org.broadinstitute.gatk.engine.walkers.ReadWalker#map(org.broadinstitute.gatk.engine.contexts.ReferenceContext, org.broadinstitute.gatk.utils.sam.GATKSAMRecord, org.broadinstitute.gatk.engine.refdata.ReadMetaDataTracker)
	 */
	@Override
	public SAMRecord map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
		int readLength = read.getReadLength();
		boolean negativeStrand = read.getReadNegativeStrandFlag();
		boolean secondPair = false;
		int numberOfPatternInRef = 0;
		
		if (read.getReadPairedFlag()) {
			secondPair = read.getSecondOfPairFlag();
			if(secondPair)
				negativeStrand = !negativeStrand;
		}
		
		if(onlyMappedReads && read.getReadUnmappedFlag())
			return null;
		if(filterPcrDupsReads && read.getDuplicateReadFlag())
			return null;
		if(filterNoPrimaryAlignedReads && read.getNotPrimaryAlignmentFlag())
			return null;
		if(read.getReadPairedFlag()){
			if(onlyProperPairedReads && !read.getProperPairFlag())
				return null;
			if(filterInvertedDupsReads && (read.getAlignmentStart() == read.getMateAlignmentStart() && read.getReadNegativeStrandFlag() == read.getMateNegativeStrandFlag())){
				if(filterInvertedDupsReadsComplete){
					return null;
				}
				else{
					if(secondPair){
						return null;
					}
				}
			}
				
		}
		
		
		
		byte[] bases = read.getReadBases();
		byte[] baseQs = read.getBaseQualities();
		byte[] refBases = new byte[readLength];
		
		int start = read.getAlignmentStart();
		int end = read.getAlignmentEnd();
		if(ref == null || ref.getLocus() == null || ref.getLocus().getContig() == null || start <= 1 || read.getReadUnmappedFlag()){
			numberOfPatternInRef = -1;
		}
		else{
			//GenomeLoc locBoundary = ref.getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(), start-1, end+1);

			for (int i = 0, cor = start; i < readLength; i++) {
				GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(), cor);
				//if( !ref.getWindow().containsP(loc) )
				//	continue;
				//if (!ref.getWindow().overlapsP(locBoundary)){
				if( !ref.getWindow().containsP(loc) ){
					numberOfPatternInRef = -1;
					break;
				}
					

				ReferenceContext tmpRef = new ReferenceContext(ref.getGenomeLocParser(), loc, ref.getWindow(), ref.getBases());
				if(tmpRef.getLocus().getStart() - tmpRef.getWindow().getStart() >= 0 && tmpRef.getLocus().getStart() - tmpRef.getWindow().getStart() < readLength){
					//System.err.println(read.getReadName() + "\t" + i + "\t" + readLength + "\t" + tmpRef.getLocus().getStart() + "\t" + tmpRef.getWindow().getStart() + "\t" + tmpRef.getBases().length + "\t" + new String(tmpRef.getBases()));
					refBases[i] = tmpRef.getBase();
				}
				else{
					numberOfPatternInRef = -1;
					break;
				}
				

				cor++;

			}
		}
		
		
		
		if(negativeStrand){
			bases = BaseUtils.simpleReverseComplement(bases);
			baseQs = BaseUtilsMore.simpleReverse(baseQs);
			if(numberOfPatternInRef != -1){
				refBases = BaseUtils.simpleReverseComplement(refBases);
			}
			
		}
		
		if(readLength < pattern.length())
			return null;
		int numberOfPattern = 0;
		
		byte[] patterns = pattern.getBytes();
		
		
		for(int i = 0; i <= readLength-patterns.length; i++){
			boolean flag = true;
			boolean flagRef = true;
			for(int j = i, index = 0; index < patterns.length; index++, j++){
				if(baseQs[j] < minBaseQ){
					flag = false;
					flagRef = false;
					break;
				}
				if(numberOfPatternInRef != -1){ //mean reads is mapped
					
					if(!BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(patterns[index], refBases[j])){
						flagRef = false;
						
					}
					else{
					//	System.err.println(patterns[index] + "\t" + refBases[j]);
					}
				}
				if(!BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(patterns[index], bases[j])){
					flag = false;
					break;
				}
				
			}
			if(flag){
				numberOfPattern++;
			}
			if(numberOfPatternInRef != -1 && flagRef){
				numberOfPatternInRef++;
			}
		}
		if(patternHistWriter != null){
			patternHistWriter.println(numberOfPattern + "\t" + numberOfPatternInRef + "\t" + readLength);
		}
		
		//System.err.println(read.getReadName() + "\t" + numberOfPattern + "\t" + numberOfPatternInRef);
		//if(numberOfPatternInRef == -1 || numberOfPatternInRef == 0)
		//	System.err.println(read.getReadName() + "\t" + numberOfPattern + "\t" + numberOfPatternInRef);
		
		
		//if(((numberOfPatternInRef == -1 || numberOfPatternInRef == 0) && numberOfPattern == 0) || (numberOfPattern <= minNumMethyC &&  (double)numberOfPattern/(double)numberOfPatternInRef <= minFreqMethyC)){
		if(((numberOfPatternInRef == -1 || numberOfPatternInRef == 0) && numberOfPattern == 0) || (numberOfPatternInRef != -1 && (double)numberOfPattern/(double)numberOfPatternInRef <= minFreqMethyC)){
			if(cleanBam != null)
				bamWriter.addAlignment(read);
		}
		else{
			if(dirtyBam != null)
				bamDirtyWriter.addAlignment(read);
		}
		
		return null;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.gatk.engine.walkers.Walker#reduceInit()
	 */
	@Override
	public Long reduceInit() {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.gatk.engine.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Long reduce(SAMRecord value, Long sum) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public void onTraversalDone(Long result) {
		//	output.close();
			if(cleanBam != null)
				bamWriter.close();
			if(dirtyBam != null)
				bamDirtyWriter.close();
			logger.info("Finshed!~");
	}
	

}
