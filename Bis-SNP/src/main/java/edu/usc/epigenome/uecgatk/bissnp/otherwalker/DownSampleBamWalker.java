/**
 * 
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp.otherwalker;

import java.io.File;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;

import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.By;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.Downsample;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.engine.walkers.ReadFilters;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.engine.walkers.Window;

import main.java.edu.usc.epigenome.uecgatk.bissnp.BisulfiteArgumentCollection;
import main.java.edu.usc.epigenome.uecgatk.bissnp.DownsamplingBAM;

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
	 * @see org.broadinstitute.gatk.engine.walkers.LocusWalker#map(org.broadinstitute.gatk.engine.refdata.RefMetaDataTracker, org.broadinstitute.gatk.engine.contexts.ReferenceContext, org.broadinstitute.gatk.engine.contexts.AlignmentContext)
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
	 * @see org.broadinstitute.gatk.engine.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Condition reduce(Condition value, Condition sum) {
		sum.Loci += value.Loci;
		sum.coverageSumPre += value.coverageSumPre;
		sum.coverageSumAfter += value.coverageSumAfter;
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
