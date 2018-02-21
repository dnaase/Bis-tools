/**
 * 
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp.otherwalker;

import java.io.File;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;

import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.walkers.Downsample;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.filters.DuplicateReadFilter;
import org.broadinstitute.gatk.engine.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.gatk.engine.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.gatk.engine.filters.UnmappedReadFilter;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.By;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.Downsample;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.engine.walkers.ReadFilters;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.utils.pileup.PileupElement;

import main.java.edu.usc.epigenome.uecgatk.bissnp.BisSNPUtils;
import main.java.edu.usc.epigenome.uecgatk.bissnp.BisulfiteArgumentCollection;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.BisulfiteFivePrimeConvReadsFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.BisulfiteIncompleteConvReadsFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.BisulfiteMismatchReadsFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.InvertedDupsReadFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.MappingQualityFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.NotProperPairedReadFilter;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 18, 2012 3:15:17 PM
 * 
 */
@Reference(window = @Window(start = -500, stop = 500))
@By(DataSource.READS)
@Downsample(by = DownsampleType.NONE)
@ReadFilters({UnmappedReadFilter.class, DuplicateReadFilter.class, NotPrimaryAlignmentFilter.class, FailsVendorQualityCheckFilter.class, MappingQualityFilter.class, NotProperPairedReadFilter.class, InvertedDupsReadFilter.class, BisulfiteMismatchReadsFilter.class, BisulfiteIncompleteConvReadsFilter.class, BisulfiteFivePrimeConvReadsFilter.class})
public class BamFilterWalker extends LocusWalker<BamFilterWalker.Condition, BamFilterWalker.Condition> implements TreeReducible<BamFilterWalker.Condition> {

	@ArgumentCollection
	private static BisulfiteArgumentCollection BAC = new BisulfiteArgumentCollection();
	
	@Output(fullName = "file_name_output_reads_after_downsampling", shortName = "reads", doc = "output Bam file's name that after filtering", required = false)
	public String reads = null;
	
	
	private SAMFileWriter samWriter = null;

	
	public void initialize() {
		File outputBamFile = new File(reads);
		SAMFileWriterFactory samFileWriterFactory = new SAMFileWriterFactory();
		samFileWriterFactory.setCreateIndex(true);
		samWriter = samFileWriterFactory.makeBAMWriter(getToolkit().getSAMFileHeader(), false, outputBamFile);
		
	}
	
	/* (non-Javadoc)
	 * @see org.broadinstitute.gatk.engine.walkers.LocusWalker#map(org.broadinstitute.gatk.engine.refdata.RefMetaDataTracker, org.broadinstitute.gatk.engine.contexts.ReferenceContext, org.broadinstitute.gatk.engine.contexts.AlignmentContext)
	 */
	@Override
	public Condition map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
		Condition value = new Condition();
		String tag = "Xi";
		int covergaeLimit = BAC.orcad;
		int coverageMarked = 0;
		for (PileupElement p : rawContext.getBasePileup().getOverlappingFragmentFilteredPileup()) {
			if (p.getRead().getIntegerAttribute(tag) != null) {
				continue;
			}
			value.readsTotal++;
			if (!BisSNPUtils.goodBaseInPileupElement(p, BAC, refContext)) {
				continue;
			}
			
			if (coverageMarked >= covergaeLimit)
				break;
			if (p.getRead().getIntegerAttribute(tag) == null) {
				samWriter.addAlignment(p.getRead());
				p.getRead().setAttribute(tag, 1);
				coverageMarked++;
				value.readsAfterFilter ++;
			}
			
		}
		
		return value;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.gatk.engine.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Condition reduce(Condition value, Condition sum) {
		sum.readsTotal += value.readsTotal;
		sum.readsAfterFilter += value.readsAfterFilter;
		return sum;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.gatk.engine.walkers.Walker#reduceInit()
	 */
	@Override
	public Condition reduceInit() {
		// TODO Auto-generated method stub
		
		return new Condition();
	}
	
	/* (non-Javadoc)
	 * @see org.broadinstitute.gatk.engine.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Condition treeReduce(Condition lhs, Condition rhs) {
		lhs.readsTotal += rhs.readsTotal;
		lhs.readsAfterFilter += rhs.readsAfterFilter;
		return lhs;

	}
	
	@Override
	public void onTraversalDone(Condition sum) {
		samWriter.close();

		logger.info("Reads after filtering: " + sum.readsAfterFilter);
	}
	
	public class Condition{
		public long readsTotal = 0;
		public long readsAfterFilter = 0;
		
	}

}
