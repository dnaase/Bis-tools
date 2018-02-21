/**
 * This is to implement Robert Shoemaker 2010 Genome Research paper's algorithm's comparison with
 * BisSNP discard reads in multiple location, bad mates were discarded. T/C ratio bigger than A/G
 * ratio is identified as bisulfite-conversion reads, then non-bisulfite-conversion read were
 * converted to its reverse complement reads. Cs in reads are converted to Ts, this is called in
 * silico demethylation.
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp.otherwalker;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.MissingResourceException;
import java.util.ResourceBundle;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;

import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.filters.BadMateFilter;
import org.broadinstitute.gatk.engine.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.gatk.engine.filters.UnmappedReadFilter;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.ReadFilters;
import org.broadinstitute.gatk.engine.walkers.ReadWalker;
import org.broadinstitute.gatk.engine.walkers.Requires;
import org.broadinstitute.gatk.engine.walkers.WalkerName;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.text.TextFormattingUtils;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Apr 29, 2012 2:48:50 PM
 * 
 */
@WalkerName("ShoemakerReadsConversion")
@Requires({ DataSource.READS })
@ReadFilters({ UnmappedReadFilter.class, BadMateFilter.class, NotPrimaryAlignmentFilter.class })
public class ShoemakerReadsConversionWalker extends ReadWalker<SAMRecord, SAMFileWriter> {

	public static final String PROGRAM_RECORD_NAME = "ShoemakerReadsConversionWalker";
	/**
	 * 
	 */

	@Output(fullName = "output", shortName = "o", doc = "output Bam file's name", required = true)
	private String OUTPUT_BAM = null;

	private SAMFileWriter samWriter = null;

	@Override
	public void initialize() {
		final SAMFileHeader header = getToolkit().getSAMFileHeader().clone();

		final SAMProgramRecord programRecord = new SAMProgramRecord(PROGRAM_RECORD_NAME);
		final ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText",null);
		try {
			final String version = headerInfo.getString("org.broadinstitute.gatk.engine.version");
			programRecord.setProgramVersion(version);
		} catch (MissingResourceException e) {
		}

		StringBuffer sb = new StringBuffer();
		sb.append(getToolkit().createApproximateCommandLineArgumentString(getToolkit(), this));

		programRecord.setCommandLine(sb.toString());

		List<SAMProgramRecord> oldRecords = header.getProgramRecords();
		List<SAMProgramRecord> newRecords = new ArrayList<SAMProgramRecord>(oldRecords.size() + 1);
		for (SAMProgramRecord record : oldRecords) {
			if (!record.getId().startsWith(PROGRAM_RECORD_NAME))
				newRecords.add(record);
		}
		newRecords.add(programRecord);
		header.setProgramRecords(newRecords);

		// Write out the new header
		samWriter = new SAMFileWriterFactory().makeBAMWriter(header,false,new File(OUTPUT_BAM));
	}

	@Override
	public SAMRecord map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
		boolean isBsConvertedStrand = isBisulfiteConversionStrand(read);
		if (isBsConvertedStrand) {
			read.setReadNegativeStrandFlag(!read.getReadNegativeStrandFlag());
		}
		GATKSAMRecord newRead = demethylatedReads(read);

		return newRead;
	}

	@Override
	public void onTraversalDone(SAMFileWriter output) {
		// output.close();

	}

	@Override
	public SAMFileWriter reduce(SAMRecord read, SAMFileWriter output) {
		if (output != null) {
			output.addAlignment(read);
		}
		return output;
	}

	@Override
	public SAMFileWriter reduceInit() {
		// TODO Auto-generated method stub
		return samWriter;
	}

	private double basesRatioInReads(GATKSAMRecord read, byte a, byte b) {
		byte[] bases = read.getReadNegativeStrandFlag() ? BaseUtils.simpleReverseComplement(read.getReadBases()) : read.getReadBases();
		int countBaseA = 0;
		int countBaseB = 0;
		for (byte base : bases) {
			if (BaseUtils.basesAreEqual(a, base)) {
				countBaseA++;
			}
			if (BaseUtils.basesAreEqual(b, base)) {
				countBaseB++;
			}
		}
		return (double) countBaseA / (double) countBaseB;
	}

	private GATKSAMRecord demethylatedReads(GATKSAMRecord read) {
		byte[] bases = read.getReadBases();
		if (read.getReadNegativeStrandFlag()) {
			for (int i = 0; i < bases.length; i++) {
				// in negative strand, gatk return 'c' as 'g'...
				if (BaseUtils.basesAreEqual(bases[i], BaseUtils.Base.G.base)) {
					bases[i] = BaseUtils.Base.A.base;
				}
			}
		} else {
			for (int i = 0; i < bases.length; i++) {
				if (BaseUtils.basesAreEqual(bases[i], BaseUtils.Base.C.base)) {
					bases[i] = BaseUtils.Base.T.base;
				}
			}
		}

		read.setReadBases(bases);
		return read;
	}

	private boolean isBisulfiteConversionStrand(GATKSAMRecord read) {
		double tcRatio = basesRatioInReads(read, BaseUtils.Base.T.base, BaseUtils.Base.C.base);
		double agRatio = basesRatioInReads(read, BaseUtils.Base.A.base, BaseUtils.Base.G.base);
		return tcRatio > agRatio ? true : false;

	}
}
