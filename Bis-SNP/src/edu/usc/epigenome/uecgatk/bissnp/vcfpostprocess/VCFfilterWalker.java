/**
 * 
 */
package edu.usc.epigenome.uecgatk.bissnp.vcfpostprocess;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time May 6, 2012 10:24:45 PM
 * 
 */
// @By(DataSource.REFERENCE_ORDERED_DATA)
@By(DataSource.READS)
@Requires({ DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_ORDERED_DATA, DataSource.REFERENCE_BASES })
@PartitionBy(PartitionType.LOCUS)
public class VCFfilterWalker extends LocusWalker<Integer, Integer> implements TreeReducible<Integer> {

	@Output(fullName = "filteredMinorAllele", shortName = "filteredMinorAllele", doc = "output filtered Minor Allele", required = false)
	public PrintStream filteredMinorAllele = null;

	@Output(fullName = "new_vcf", shortName = "newVcf", doc = "filtered vcf file", required = true)
	public VCFWriter newVcf = null;

	/**
	 * 
	 */
	@Input(fullName = "old_vcfs", shortName = "oldVcfs", doc = "input vcf file", required = true)
	public List<RodBinding<VariantContext>> oldVcfs;

	@Input(fullName = "snp_vcfs", shortName = "snpVcf", doc = "input SNP vcf file, used to filter out adjacent SNP files", required = true)
	// used to filter SNP clustered
	public RodBinding<VariantContext> snpVcf;

	@Argument(shortName = "filterOneAlleleRead", doc = "filter out loci that have only one allele's reads(for heterozygous SNPs statistics only)", required = false)
	protected boolean filterOneAlleleRead = false;

	@Argument(shortName = "filterSB", doc = "filter out loci that have high Strand bias", required = false)
	protected boolean filterSB = false;

	@Argument(shortName = "maxCov", doc = "maximum covergae required for the position in VCF file, default:120", required = false)
	protected int maxCov = 120;
	@Argument(shortName = "minBaseQ", doc = "Minimum base quality in the covered loci when counting coverage, default:0", required = false)
	protected int minBaseQ = 0;
	@Argument(shortName = "minCov", doc = "minimum covergae required for the position in VCF file, default:1", required = false)
	protected int minCov = 1;
	@Argument(shortName = "minMapQ", doc = "Minimum mapping quality in the covered loci when counting coverage, default:30", required = false)
	protected int minMapQ = 30;
	@Argument(shortName = "minNumOneStrand", doc = "Minimum number of one strand only when less than 10 coverages, default:2", required = false)
	protected int minNumOneStrand = 2;
	@Argument(shortName = "minPercentOneStrand", doc = "Minimum percentage of one strand only when total sequence coverage > 10X, default:0.2", required = false)
	protected double minPercentOneStrand = 0.2;
	@Argument(shortName = "minReadsMinorAllele", doc = "minimum number of reads in minor allele, default:1", required = false)
	protected int minReadsMinorAllele = 1;
	@Argument(shortName = "minSNPinWind", doc = "minimum number of SNPs in the window, default:2", required = false)
	protected int minSNPinWind = 2;

	@Argument(shortName = "windSizeForSNPfilter", doc = "window size for detect SNP cluster, default:10, means +/- 10bp distance, no second SNP there", required = false)
	protected int windSizeForSNPfilter = 10;

	@Override
	public void initialize() {
		Map<String, VCFHeader> headers = VCFUtils.getVCFHeadersFromRods(getToolkit(), oldVcfs);
		for (VCFHeader header : headers.values()) {
			newVcf.writeHeader(header);
		}
	}

	@Override
	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
		// TODO Auto-generated method stub
		if (tracker == null) {
			return null;
		}
		// First, filter the VariantContext to represent only the samples for
		// evaluation
		for (RodBinding<VariantContext> oldVcf : oldVcfs) {
			for (VariantContext vc_input : tracker.getValues(oldVcf, ref.getLocus())) {
				if (vc_input != null) {
					// oldRecord++;
					if (passFilter(context, ref, vc_input)) {
						// newRecord++;

						newVcf.add(vc_input);
					}
				}
			}
		}

		return null;
	}

	@Override
	public void onTraversalDone(Integer value) {

		logger.info("Finished!");
	}

	@Override
	public Integer reduce(Integer value, Integer sum) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Integer reduceInit() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Integer treeReduce(Integer lhs, Integer rhs) {
		// TODO Auto-generated method stub
		return null;
	}

	private boolean passFilter(AlignmentContext context, ReferenceContext ref, VariantContext vc) {
		HashMap<Byte, Integer> baseContainer = new HashMap<Byte, Integer>();
		if (context.hasBasePileup()) {
			int coverage = 0;
			int negStarndCovFirstEnd = 0;
			int posStarndCovFirstEnd = 0;
			int negStarndCovSecondEnd = 0;
			int posStarndCovSecondEnd = 0;
			int negStarndCov = 0;
			boolean paired = false;

			for (PileupElement p : context.getBasePileup()) {

				SAMRecord samRecord = p.getRead();
				paired = samRecord.getReadPairedFlag();
				if (samRecord.getDuplicateReadFlag() || samRecord.getNotPrimaryAlignmentFlag() || samRecord.getReadFailsVendorQualityCheckFlag() || samRecord.getReadUnmappedFlag()
						|| samRecord.getMappingQuality() < minMapQ || p.getQual() < minBaseQ)
					continue;
				if (samRecord.getReadPairedFlag()) {
					if (samRecord.getAlignmentStart() == samRecord.getMateAlignmentStart() && samRecord.getReadNegativeStrandFlag() == samRecord.getMateNegativeStrandFlag()) {
						if (samRecord.getSecondOfPairFlag())
							continue;
					}
					if (!samRecord.getProperPairFlag())
						continue;
					if (samRecord.getReadNegativeStrandFlag()) {
						if (samRecord.getSecondOfPairFlag()) {
							negStarndCovSecondEnd++;
							if (BaseUtils.basesAreEqual(ref.getBase(), BaseUtils.C) && BaseUtils.basesAreEqual(p.getBase(), BaseUtils.T)) {
								Integer num = 0;
								if (baseContainer.containsKey(BaseUtils.C))
									num = baseContainer.get(BaseUtils.C);
								num++;
								baseContainer.put(BaseUtils.C, num);

							} else {
								Integer num = 0;
								if (baseContainer.containsKey(p.getBase()))
									num = baseContainer.get(p.getBase());
								num++;
								baseContainer.put(p.getBase(), num);
							}
						} else {
							negStarndCovFirstEnd++;
							if (BaseUtils.basesAreEqual(ref.getBase(), BaseUtils.G) && BaseUtils.basesAreEqual(p.getBase(), BaseUtils.A)) {
								Integer num = 0;
								if (baseContainer.containsKey(BaseUtils.G))
									num = baseContainer.get(BaseUtils.G);
								num++;
								baseContainer.put(BaseUtils.G, num);
							} else {
								Integer num = 0;
								if (baseContainer.containsKey(p.getBase()))
									num = baseContainer.get(p.getBase());
								num++;
								baseContainer.put(p.getBase(), num);
							}
						}
					} else {
						if (samRecord.getSecondOfPairFlag()) {
							posStarndCovSecondEnd++;
							if (BaseUtils.basesAreEqual(ref.getBase(), BaseUtils.G) && BaseUtils.basesAreEqual(p.getBase(), BaseUtils.A)) {
								Integer num = 0;
								if (baseContainer.containsKey(BaseUtils.G))
									num = baseContainer.get(BaseUtils.G);
								num++;
								baseContainer.put(BaseUtils.G, num);
							} else {
								Integer num = 0;
								if (baseContainer.containsKey(p.getBase()))
									num = baseContainer.get(p.getBase());
								num++;
								baseContainer.put(p.getBase(), num);
							}
						} else {
							posStarndCovFirstEnd++;
							if (BaseUtils.basesAreEqual(ref.getBase(), BaseUtils.C) && BaseUtils.basesAreEqual(p.getBase(), BaseUtils.T)) {
								Integer num = 0;
								if (baseContainer.containsKey(BaseUtils.C))
									num = baseContainer.get(BaseUtils.C);
								num++;
								baseContainer.put(BaseUtils.C, num);
							} else {
								Integer num = 0;
								if (baseContainer.containsKey(p.getBase()))
									num = baseContainer.get(p.getBase());
								num++;
								baseContainer.put(p.getBase(), num);
							}
						}
					}
				} else {
					if (samRecord.getReadNegativeStrandFlag()) {
						negStarndCov++;
					}
				}
				coverage++;

			}
			if (coverage >= minCov && coverage <= maxCov) {
				if (filterOneAlleleRead) {
					if (baseContainer.keySet().size() < 2) {

						return false;
					} else {
						int tmpNum = Integer.MAX_VALUE;
						byte minorAllele = BaseUtils.NO_CALL_INDEX;
						int totalCoverage = 0;
						for (Byte base : baseContainer.keySet()) {
							if (baseContainer.get(base) <= tmpNum) {
								minorAllele = base;
								tmpNum = baseContainer.get(base);
							}
							totalCoverage += baseContainer.get(base);

						}

						if (tmpNum <= minReadsMinorAllele) {
							if (filteredMinorAllele != null) {
								filteredMinorAllele.println(ref.getLocus().getContig() + "\t" + ref.getLocus().getStart() + "\t" + (char) ref.getBase() + "\t" + (char) minorAllele + "\t" + tmpNum
										+ "\t" + totalCoverage + "\t" + vc.toString());
							}
							return false;
						}

						// System.err.println(ref.getLocus().getStart());
					}
				}
				if (filterSB) {
					if (paired) {
						if ((coverage > 10 && (double) (negStarndCovFirstEnd + posStarndCovSecondEnd) / (double) coverage > minPercentOneStrand && (double) (negStarndCovSecondEnd + posStarndCovFirstEnd)
								/ (double) coverage > minPercentOneStrand)
								|| (coverage <= 10 && (negStarndCovFirstEnd + posStarndCovSecondEnd) > minNumOneStrand && (negStarndCovSecondEnd + posStarndCovFirstEnd) > minNumOneStrand)) {
							return true;
						}
					} else {
						if ((coverage > 10 && (double) negStarndCov / (double) coverage > minPercentOneStrand && (double) (coverage - negStarndCov) / (double) coverage > minPercentOneStrand)
								|| (coverage <= 10 && negStarndCov > minNumOneStrand && (coverage - negStarndCov) > minNumOneStrand)) {
							return true;
						}
					}

				} else {
					return true;
				}

			}
		}
		return false;
	}

}
