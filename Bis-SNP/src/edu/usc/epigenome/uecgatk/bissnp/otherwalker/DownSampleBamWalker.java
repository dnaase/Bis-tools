/**
 * 
 */
package edu.usc.epigenome.uecgatk.bissnp.otherwalker;

import java.io.File;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;

import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;

import edu.usc.epigenome.uecgatk.bissnp.BisulfiteArgumentCollection;
import edu.usc.epigenome.uecgatk.bissnp.DownsamplingBAM;
import edu.usc.epigenome.uecgatk.bissnp.BisulfiteGenotyper.ContextCondition;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 18, 2012 3:15:17 PM
 * 
 */
@Reference(window = @Window(start = -500, stop = 500))
@By(DataSource.REFERENCE)
@Downsample(by = DownsampleType.NONE)
@ReadFilters()
public class DownSampleBamWalker extends LocusWalker<DownSampleBamWalker.Condition, DownSampleBamWalker.Condition> implements TreeReducible<DownSampleBamWalker.Condition> {

	@ArgumentCollection
	private static BisulfiteArgumentCollection BAC = new BisulfiteArgumentCollection();
	
	@Output(fullName = "file_name_output_reads_after_downsampling", shortName = "reads", doc = "output Bam file's name that after downsapling", required = false)
	public String reads = null;
	
	@Output(fullName = "coverage_of_sample_before_downsampling", shortName = "coverage", doc = "coverage_of_sample_before_downsampling", required = false)
	public int coverage = 13;
	
	private DownsamplingBAM downsampleBam = null;
	private SAMFileWriter samWriter = null;

	
	public void initialize() {
		File outputBamFile = new File(reads);
		SAMFileWriterFactory samFileWriterFactory = new SAMFileWriterFactory();
		samFileWriterFactory.setCreateIndex(true);
		samWriter = samFileWriterFactory.makeBAMWriter(getToolkit().getSAMFileHeader(), false, outputBamFile);
		downsampleBam = new DownsamplingBAM(BAC, samWriter, coverage);
	}
	
	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.LocusWalker#map(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */
	@Override
	public Condition map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
		int coverageAfterDownsampling = downsampleBam.downsamplingBamFile(rawContext);
		Condition value = new Condition();
		value.coverageSumPre = rawContext.getBasePileup().depthOfCoverage();
		value.coverageSumAfter = coverageAfterDownsampling;
		value.Loci = 1;
		return value;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Condition reduce(Condition value, Condition sum) {
		sum.Loci += value.Loci;
		sum.coverageSumPre += value.coverageSumPre;
		sum.coverageSumAfter += value.coverageSumAfter;
		return sum;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduceInit()
	 */
	@Override
	public Condition reduceInit() {
		// TODO Auto-generated method stub
		return new Condition();
	}
	
	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Condition treeReduce(Condition lhs, Condition rhs) {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public void onTraversalDone(Condition sum) {
		samWriter.close();
		logger.info("Loci traversaled: " + sum.Loci);
		logger.info("Covergae before downsampling: " + (double)sum.coverageSumPre/(double)sum.Loci);
		logger.info("Covergae after downsampling: " + (double)sum.coverageSumAfter/(double)sum.Loci);
	}
	
	public class Condition{
		public long coverageSumPre = 0;
		public long coverageSumAfter = 0;
		public long Loci = 0;
		
	}

}
