/**
 * NOTES: Only use good reads!!
 * quickly check the cytosine context methylation level provided, e.g. use the genome wide CG_methylation_sum/ CG_number, CG is determined by reference genome
 * also, check output the methylation value along the sequence cycle: 2nd end in paired end is -, while 1st end or single end is always + coordinate
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp.otherwalker;

import java.io.PrintStream;
import java.util.TreeMap;

import main.java.edu.usc.epigenome.uecgatk.bissnp.BaseUtilsMore;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.InvertedDupsReadFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.NotProperPairedReadFilter;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.filters.BadMateFilter;
import org.broadinstitute.gatk.engine.filters.DuplicateReadFilter;
import org.broadinstitute.gatk.engine.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.gatk.engine.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.gatk.engine.filters.UnmappedReadFilter;
import org.broadinstitute.gatk.engine.filters.MappingQualityFilter;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
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
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.BisulfiteMismatchReadsFilter;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jul 27, 2012 7:51:31 PM
 * 
 */
@WalkerName("QuickMethylationLevel")
@Reference(window=@Window(start=-500,stop=500))
@By(DataSource.READS)
@Requires({DataSource.REFERENCE, DataSource.READS})
@ReadFilters( {UnmappedReadFilter.class, MappingQualityFilter.class, BadMateFilter.class, NotPrimaryAlignmentFilter.class, DuplicateReadFilter.class, FailsVendorQualityCheckFilter.class, InvertedDupsReadFilter.class, NotProperPairedReadFilter.class, BisulfiteMismatchReadsFilter.class} )

public class QuickMethylationLevelWalker extends ReadWalker<QuickMethylationLevelWalker.Datum, QuickMethylationLevelWalker.Datum> {

	@Input(fullName = "pattern_to_search", shortName = "pattern", doc = "cytosines pattern to search in provided BAM", required = false)
	public String pattern = "ACT";
	
	@Output(fullName = "histogram_methylated_pattern_along_cycle", shortName = "patternHist", doc = "output histogram of number of methylated pattern(ACA) along the cycle", required = true)
	public PrintStream patternHistWriter = null;
	
	@Argument(fullName = "methylation_count_position", shortName = "methyPos", doc = "methylation count position in the pattern provided, 1-based. Default: -methyPos 2 for ACT pattern", required = false)
    public int methyPos = 2;
	
	@Argument(fullName = "minimum_base_quality", shortName = "minBaseQ", doc = "minimum base quality score required to check. Default: 5", required = false)
    public int minBaseQ = 5;
	
	

	

	/* (non-Javadoc)
	 * @see org.broadinstitute.gatk.engine.walkers.ReadWalker#map(org.broadinstitute.gatk.engine.contexts.ReferenceContext, org.broadinstitute.gatk.utils.sam.GATKSAMRecord, org.broadinstitute.gatk.engine.refdata.ReadMetaDataTracker)
	 */
	@Override
	public Datum map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
		int readLength = read.getReadLength();
		boolean negativeStrand = read.getReadNegativeStrandFlag();
		boolean secondPair = false;

		
		if (read.getReadPairedFlag()) {
			secondPair = read.getSecondOfPairFlag();
			if(secondPair)
				negativeStrand = !negativeStrand;
		}
	
		
		byte[] bases = read.getReadBases();
		byte[] baseQs = read.getBaseQualities();
		byte[] refBases = new byte[readLength];
		
		int start = read.getAlignmentStart();
		int end = read.getAlignmentEnd();
		if(ref == null || ref.getLocus() == null || ref.getLocus().getContig() == null || start <= 1){
			return null;
		}
		else{
			for (int i = 0, cor = start; i < readLength; i++) {
				GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(), cor);
				//if( !ref.getWindow().containsP(loc) )
				//	continue;
				//if (!ref.getWindow().overlapsP(locBoundary)){
				if( !ref.getWindow().containsP(loc) ){
					return null;
				}
					

				ReferenceContext tmpRef = new ReferenceContext(ref.getGenomeLocParser(), loc, ref.getWindow(), ref.getBases());
				if(tmpRef.getLocus().getStart() - tmpRef.getWindow().getStart() >= 0 && tmpRef.getLocus().getStart() - tmpRef.getWindow().getStart() < readLength){
					//System.err.println(read.getReadName() + "\t" + i + "\t" + readLength + "\t" + tmpRef.getLocus().getStart() + "\t" + tmpRef.getWindow().getStart() + "\t" + tmpRef.getBases().length + "\t" + new String(tmpRef.getBases()));
					refBases[i] = tmpRef.getBase();
				}
				else{
					return null;
				}
				

				cor++;

			}
		}
		
		
		
		if(negativeStrand){
			bases = BaseUtils.simpleReverseComplement(bases);
			baseQs = BaseUtilsMore.simpleReverse(baseQs);
			
			refBases = BaseUtils.simpleReverseComplement(refBases);

			
		}
		
		if(readLength < pattern.length())
			return null;
		long numberOfPattern = 0;
		long numberOfPatternMethy = 0;
		
		byte[] patterns = pattern.getBytes();
		Datum value = new Datum(pattern);
		
		for(int i = 0; i <= readLength-patterns.length; i++){
			boolean flag = true;
			long flagMethy = -1;
			for(int j = i, index = 0; index < patterns.length; index++, j++){
				if(baseQs[j] < minBaseQ){
					break;
				}
	
				if(!BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(patterns[index], refBases[j])){
					flag = false;
					break;	
				}
				else{
					if(index == (methyPos-1)){
						if(BaseUtils.basesAreEqual(bases[j], BaseUtils.Base.C.base)){
							flagMethy=1;
						}
						else if(BaseUtils.basesAreEqual(bases[j], BaseUtils.Base.T.base)){
							flagMethy=0;
						}
						else{
							flagMethy=-1;
						}
						
					}
				}

				
			}
			if(flag && flagMethy != -1){
				numberOfPattern++;
				numberOfPatternMethy += flagMethy;
				int pos = i + methyPos;
				if(secondPair)
					pos = 0-pos;
				value.motifStatAlongCycle.put(pos, new Pair<Long, Long>(flagMethy, 1L));
			}
			
		}
		value.motifStat = new Pair<Long, Long>(numberOfPatternMethy, numberOfPattern);
		return value;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.gatk.engine.walkers.Walker#reduceInit()
	 */
	@Override
	public Datum reduceInit() {
		Datum tmp = new Datum(pattern);
		tmp.motifStat = new Pair<Long,Long>(0L,0L);
		return tmp;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.gatk.engine.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Datum reduce(Datum value, Datum sum) {
		if(value != null){
			sum.motifStat.first +=  value.motifStat.first;
			sum.motifStat.second +=  value.motifStat.second;
			if(!value.motifStatAlongCycle.isEmpty()){
				for(Integer key : value.motifStatAlongCycle.keySet()){
					if(sum.motifStatAlongCycle.containsKey(key)){
						long methy = sum.motifStatAlongCycle.get(key).first + value.motifStatAlongCycle.get(key).first;
						long patternNum = sum.motifStatAlongCycle.get(key).second + value.motifStatAlongCycle.get(key).second;
						sum.motifStatAlongCycle.put(key, new Pair<Long, Long>(methy,patternNum));
					}
					else{
						long methy = value.motifStatAlongCycle.get(key).first;
						long patternNum = value.motifStatAlongCycle.get(key).second;
						sum.motifStatAlongCycle.put(key, new Pair<Long, Long>(methy,patternNum));
					}
				}
			}
			
		}
		return sum;
	}
	
	
	public void onTraversalDone(Datum result) {
		if(result != null){
			logger.info("Total methylated pattern(" + pattern + ") number is: " + result.motifStat.first);
			logger.info("Total pattern(" + pattern + ") number in the reference genome is: " + result.motifStat.second);
			logger.info("% of methylation level of pattern(" + pattern + ") number in the reference genome is: " + String.format("%.2f", 100 * (double)result.motifStat.first / (double)result.motifStat.second));
			if(!result.motifStatAlongCycle.isEmpty()){
				for(Integer key : result.motifStatAlongCycle.keySet()){
					patternHistWriter.println(key + "\t" + result.motifStatAlongCycle.get(key).first + "\t" + result.motifStatAlongCycle.get(key).second + "\t" + String.format("%.2f", 100 * (double)result.motifStatAlongCycle.get(key).first / (double)result.motifStatAlongCycle.get(key).second));
				}
			}
		}
		
		logger.info("Finshed!~");
	}
	
	public class Datum{

		public Pair<Long,Long> motifStat = null;
		public TreeMap<Integer,Pair<Long,Long>> motifStatAlongCycle = null; //paired methylated pattern number, total pattern number 
		public Datum(String motifs){

			motifStatAlongCycle =  new TreeMap<Integer,Pair<Long,Long>>();

		}
		
		
	}

}
