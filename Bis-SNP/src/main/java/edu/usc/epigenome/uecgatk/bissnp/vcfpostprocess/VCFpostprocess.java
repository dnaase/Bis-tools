/**
 * After BisSNP raw genotyping, this walker was used to generate filtered SNPs by apply some hard
 * filters, then get into VCFstatisticsWalker to get summary statistics for different cytosines
 * (maybe it is better to integrated here..).
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp.vcfpostprocess;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.SAMSequenceDictionary;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.gatk.utils.refdata.utils.RODRecordList;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

import main.java.edu.usc.epigenome.uecgatk.bissnp.BaseUtilsMore;
import main.java.edu.usc.epigenome.uecgatk.bissnp.BisSNPUtils;
import main.java.edu.usc.epigenome.uecgatk.bissnp.BisulfiteVCFConstants;
import main.java.edu.usc.epigenome.uecgatk.bissnp.writer.SortingTcgaVCFWriter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.writer.TcgaVCFWriter;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time May 15, 2012 5:57:06 PM
 * 
 */
@Reference(window = @Window(start = -200, stop = 200))
public class VCFpostprocess extends RodWalker<VCFpostprocess.VCFplusRef, VCFpostprocess.VCFCondition> implements TreeReducible<VCFpostprocess.VCFCondition> {

	@Argument(fullName = "cytosine_contexts_checked", shortName = "C", doc = "Specify the cytosine contexts to check (e.g. -C CG,1 -C CH,1... You could specify '-C' multiple times for different cytosine pattern). default: CG,1, CH,1", required = false)
	public List<String> cytosineContextsChecked = new ArrayList<String>();

	@Argument(fullName = "max_coverage", shortName = "maxCov", doc = "maximum coverage filter for heterozygous SNP, default: 250", required = false)
	public int maxCov = 250;

	@Argument(fullName = "min_bq", shortName = "minBQ", doc = "minimum base quality for both of strand, default: 10, not use this option yet", required = false)
	public int minBQ = 10;

	@Argument(fullName = "min_ct_coverage", shortName = "minCT", doc = "minimum number of CT reads for count methylation level, default: 0", required = false)
	// but when statitics in methylation, numC+numT==0 sites are not included in
	// statitics
	public int minCT = 0;

	@Argument(fullName = "mapping_quality_zero", shortName = "mq0", doc = "fraction of mapping_quality_zero filter for heterozygous SNP, default: 0.1", required = false)
	public double mq0 = 0.1;

	@Output(fullName = "new_vcf", shortName = "newVcf", doc = "filtered vcf file, without Strand bias, no clustered SNPs, above quality threshold", required = true)
	public String newVcf = null;

	/**
	 * 
	 */
	@Input(fullName = "old_vcf", shortName = "oldVcf", doc = "input vcf file", required = true)
	public RodBinding<VariantContext> oldVcf;

	@Output(doc = "Output summary statistics", required = true)
	public String outFile;

	@Argument(fullName = "quality_by_depth", shortName = "qd", doc = "quality by depth filter for heterozygous SNP, default: 1.0", required = false)
	public double qd = 1.0;

	@Argument(fullName = "genotype_qual", shortName = "qual", doc = "genotype quality score filter for heterozygous SNP, default: 20", required = false)
	public double qual = 20;

	@Argument(fullName = "strand_bias", shortName = "sb", doc = "strand bias filter for heterozygous SNP, default: -0.02", required = false)
	public double sb = -0.02;

	@Input(fullName = "snp_vcf", shortName = "snpVcf", doc = "input raw SNP vcf file(not filtered by SNP cluster or SB yet), used to filter out adjacent SNPs", required = true)
	// used to filter SNP clustered
	public RodBinding<VariantContext> snpVcf;

	@Argument(shortName = "minSNPinWind", doc = "minimum number of SNPs in the window, default:2", required = false)
	protected int minSNPinWind = 2;

	protected SortingTcgaVCFWriter multiThreadWriter = null;

	@Argument(shortName = "windSizeForSNPfilter", doc = "window size for detect SNP cluster, default:10, means +/- 10bp distance, no second SNP there", required = false)
	protected int windSizeForSNPfilter = 10;

	protected TcgaVCFWriter writer = null;

	//private int cytosinePosInPattern = 0; // for CpG or CpH, it is 0; if GCH,
											// then it should be 1;

	private int MAXIMUM_CACHE_FOR_OUTPUT_VCF = 10000000;

	private ReferenceOrderedDataSource rodIt = null;

	private PrintStream statWriter = null;

	@Override
	public void initialize() {

		rodIt = getToolkit().getRodDataSources().get(1);
		SAMSequenceDictionary refDict = getToolkit().getMasterSequenceDictionary();

		File newVcfFile = new File(newVcf);
		try{
			writer = new TcgaVCFWriter(newVcfFile, new FileOutputStream(newVcfFile), refDict,
					true, false, true,false);
		}catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		if (getToolkit().getArguments().numberOfCPUThreadsPerDataThread > 1) {
			multiThreadWriter = new SortingTcgaVCFWriter(writer, MAXIMUM_CACHE_FOR_OUTPUT_VCF);

		}
		List<RodBinding<VariantContext>> oldVcfs = new ArrayList<RodBinding<VariantContext>>();
		oldVcfs.add(oldVcf);

		Map<String, VCFHeader> headers = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), oldVcfs);
		for (VCFHeader header : headers.values()) {
			if (getToolkit().getArguments().numberOfCPUThreadsPerDataThread > 1) {
				multiThreadWriter.writeHeader(header);
			} else {
				writer.writeHeader(header);
			}

		}

		try {
			statWriter = new PrintStream(new File(outFile));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if (cytosineContextsChecked.isEmpty()) {
			cytosineContextsChecked.add("CG,1");
			cytosineContextsChecked.add("CH,1");
		}
	}

	@Override
	public VCFplusRef map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
		// TODO Auto-generated method stub
		if (tracker == null)
			return null;
		List<VariantContext> oldVcf_bindings = tracker.getValues(oldVcf);

		if (!oldVcf_bindings.isEmpty()) {

			VariantContext vc = oldVcf_bindings.get(0);

			if (passFilter(vc)) {
				// one genotype passed, then it pass all of genotypes in this location
				for (Genotype genotype : vc.getGenotypes()) {

					if (!genotype.isHomRef()) {
						if (passSNPclusterFilter(ref)) {
							if (getToolkit().getArguments().numberOfCPUThreadsPerDataThread > 1) {
								multiThreadWriter.add(vc);
							} else {
								writer.add(vc);
							}
							return new VCFplusRef(ref, vc, true);
						}
					} else if (genotype.getAttributeAsString(BisulfiteVCFConstants.C_STRAND_KEY, ".").equalsIgnoreCase("+")) {
						for (String content : cytosineContextsChecked) {
							String[] contentsArray = content.split(",");
							String cytosine = contentsArray[0].toUpperCase();
							int cytosinePosInPattern = Integer.parseInt(contentsArray[1]) - 1;
							
							if (BaseUtilsMore.iupacCodeOverlap(genotype.getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, "."), cytosine)
									&& !BisSNPUtils.isRefCytosinePattern(ref, cytosine, false)) {
								if (passSNPclusterFilter(ref)) {
									for (int i = 0 - cytosinePosInPattern; i < cytosine.length() - cytosinePosInPattern; i++) {
										if (!passSNPfilter(ref, i)) {
											return new VCFplusRef(ref, vc, false);
										}

									}
									if (getToolkit().getArguments().numberOfCPUThreadsPerDataThread > 1) {
										multiThreadWriter.add(vc);
									} else {
										writer.add(vc);
									}
									return new VCFplusRef(ref, vc, true);
								}
							}
						}
					} else if (genotype.getAttributeAsString(BisulfiteVCFConstants.C_STRAND_KEY, ".").equalsIgnoreCase("-")) {
						for (String content : cytosineContextsChecked) {
							String[] contentsArray = content.split(",");
							String cytosine = contentsArray[0].toUpperCase();
							int cytosinePosInPattern = Integer.parseInt(contentsArray[1]) - 1;
							if (BaseUtilsMore.iupacCodeOverlap(genotype.getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, "."), cytosine)
									&& !BisSNPUtils.isRefCytosinePattern(ref, cytosine, true)) {
								if (passSNPclusterFilter(ref)) {
									for (int i = 0 - cytosinePosInPattern; i < cytosine.length() - cytosinePosInPattern; i++) {
										if (!passSNPfilter(ref, 0 - i)) {
											return new VCFplusRef(ref, vc, false);
										}

									}
									if (getToolkit().getArguments().numberOfCPUThreadsPerDataThread > 1) {
										multiThreadWriter.add(vc);
									} else {
										writer.add(vc);
									}
									return new VCFplusRef(ref, vc, true);
								}
							}
						}
					} else {
						if (getToolkit().getArguments().numberOfCPUThreadsPerDataThread > 1) {
							multiThreadWriter.add(vc);
						} else {
							writer.add(vc);
						}
						return new VCFplusRef(ref, vc, true);
					}
				}
			}
		}

		return null;
	}

	@Override
	public void onTraversalDone(VCFCondition result) {

		statWriter.println("nCpgTotal\t" + result.nCpgTotal + "\tsumCpgTotal\t" + 100 * result.sumCpgTotal / result.nCpgTotal);
		statWriter.println("nCphTotal\t" + result.nCphTotal + "\tsumCphTotal\t" + 100 * result.sumCphTotal / result.nCphTotal);

		// hom cytosine:
		statWriter.println("nCpgHom\t" + result.nCpgHom + "\t" + "sumCpgHom\t" + 100 * result.sumCpgHom / result.nCpgHom);
		statWriter.println("nCpgHomRefCpg\t" + result.nCpgHomRefCpg + "\t" + "sumCpgHomRefCpg\t" + 100 * result.sumCpgHomRefCpg / result.nCpgHomRefCpg);
		statWriter.println("nCpgHomRefCph\t" + result.nCpgHomRefCph + "\t" + "sumCpgHomRefCph\t" + 100 * result.sumCpgHomRefCph / result.nCpgHomRefCph);
		statWriter.println("nCpgHomRefNotC\t" + result.nCpgHomRefNotC + "\t" + "sumCpgHomRefNotC\t" + 100 * result.sumCpgHomRefNotC / result.nCpgHomRefNotC);

		statWriter.println("nCphHom\t" + result.nCphHom + "\t" + "sumCphHom\t" + 100 * result.sumCphHom / result.nCphHom);
		statWriter.println("nCphHomRefCpg\t" + result.nCphHomRefCpg + "\t" + "sumCphHomRefCpg\t" + 100 * result.sumCphHomRefCpg / result.nCphHomRefCpg);
		statWriter.println("nCphHomRefCph\t" + result.nCphHomRefCph + "\t" + "sumCphHomRefCph\t" + 100 * result.sumCphHomRefCph / result.nCphHomRefCph);
		statWriter.println("nCphHomRefNotC\t" + result.nCphHomRefNotC + "\t" + "sumCphHomRefNotC\t" + 100 * result.sumCphHomRefNotC / result.nCphHomRefNotC);

		// het cytosines:
		statWriter.println("nCpgHet\t" + result.nCpgHet + "\t" + "sumCpgHet\t" + 100 * result.sumCpgHet / result.nCpgHet);
		statWriter.println("nCpgHetRefCpg\t" + result.nCpgHetRefCpg + "\t" + "sumCpgHetRefCpg\t" + 100 * result.sumCpgHetRefCpg / result.nCpgHetRefCpg);
		statWriter.println("nCpgHetRefCph\t" + result.nCpgHetRefCph + "\t" + "sumCpgHetRefCph\t" + 100 * result.sumCpgHetRefCph / result.nCpgHetRefCph);
		statWriter.println("nCpgHetAtC_RefCpg\t" + result.nCpgHetAtC_RefCpg + "\t" + "sumCpgHetAtC_RefCpg\t" + 100 * result.sumCpgHetAtC_RefCpg / result.nCpgHetAtC_RefCpg);
		statWriter.println("nCpgHetAtC_RefCph\t" + result.nCpgHetAtC_RefCph + "\t" + "sumCpgHetAtC_RefCph\t" + 100 * result.sumCpgHetAtC_RefCph / result.nCpgHetAtC_RefCph);
		statWriter.println("nCpgHetAtC_RefNotC\t" + result.nCpgHetAtC_RefNotC + "\t" + "sumCpgHetAtC_RefNotC\t" + 100 * result.sumCpgHetAtC_RefNotC / result.nCpgHetAtC_RefNotC);
		statWriter.println("nCpgHetAtG_RefCpg\t" + result.nCpgHetAtG_RefCpg + "\t" + "sumCpgHetAtG_RefCpg\t" + 100 * result.sumCpgHetAtG_RefCpg / result.nCpgHetAtG_RefCpg);
		statWriter.println("nCpgHetAtG_RefCph\t" + result.nCpgHetAtG_RefCph + "\t" + "sumCpgHetAtG_RefCph\t" + 100 * result.sumCpgHetAtG_RefCph / result.nCpgHetAtG_RefCph);
		statWriter.println("nCpgHetAtG_RefNotC\t" + result.nCpgHetAtG_RefNotC + "\t" + "sumCpgHetAtG_RefNotC\t" + 100 * result.sumCpgHetAtG_RefNotC / result.nCpgHetAtG_RefNotC);

		statWriter.println("nCphHet\t" + result.nCphHet + "\t" + "sumCphHet\t" + 100 * result.sumCphHet / result.nCphHet);
		statWriter.println("nCphHetRefCpg\t" + result.nCphHetRefCpg + "\t" + "sumCphHetRefCpg\t" + 100 * result.sumCphHetRefCpg / result.nCphHetRefCpg);
		statWriter.println("nCphHetRefCph\t" + result.nCphHetRefCph + "\t" + "sumCphHetRefCph\t" + 100 * result.sumCphHetRefCph / result.nCphHetRefCph);
		statWriter.println("nCphHetAtC_RefCpg\t" + result.nCphHetAtC_RefCpg + "\t" + "sumCphHetAtC_RefCpg\t" + 100 * result.sumCphHetAtC_RefCpg / result.nCphHetAtC_RefCpg);
		statWriter.println("nCphHetAtC_RefCph\t" + result.nCphHetAtC_RefCph + "\t" + "sumCphHetAtC_RefCph\t" + 100 * result.sumCphHetAtC_RefCph / result.nCphHetAtC_RefCph);
		statWriter.println("nCphHetAtC_RefNotC\t" + result.nCphHetAtC_RefNotC + "\t" + "sumCphHetAtC_RefNotC\t" + 100 * result.sumCphHetAtC_RefNotC / result.nCphHetAtC_RefNotC);

		statWriter.close();
		if (getToolkit().getArguments().numberOfCPUThreadsPerDataThread > 1) {
			multiThreadWriter.close();
		} else {
			writer.close();
		}

		logger.info("Finished!");
		logger.info(String.format("Visited Loci: %d", result.nLociVisited));

	}

	@Override
	public VCFCondition reduce(VCFplusRef value, VCFCondition sum) {
		if (value == null || !value.passFilter) {
			return sum;
		}
		sum.nLociVisited++;

		if (value.vc.hasAttribute(BisulfiteVCFConstants.CYTOSINE_TYPE) && !value.vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".").equalsIgnoreCase(".")) {

			int numC = value.vc.getGenotype(0).getAttributeAsInt(BisulfiteVCFConstants.NUMBER_OF_C_KEY, -1);
			int numT = value.vc.getGenotype(0).getAttributeAsInt(BisulfiteVCFConstants.NUMBER_OF_T_KEY, -1);
			if (numC == -1 || numT == -1 || (numC == 0 && numT == 0) || (numC + numT < minCT))
				return sum;
			double methyValue = (double) numC / (double) (numC + numT);

			String strand = value.vc.getAttributeAsString(BisulfiteVCFConstants.C_STRAND_KEY, ".");
			if (BisSNPUtils.isHomoC(value.vc)) {
				if (BisSNPUtils.isRefCpg(value.refContext)) { // ref cpg

					String pattern = value.vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".");

					if (pattern.equalsIgnoreCase("CG")) {
						sum.nCpgTotal++;
						sum.sumCpgTotal += methyValue;
						sum.nCpgHom++;
						sum.sumCpgHom += methyValue;
						sum.nCpgHomRefCpg++;
						sum.sumCpgHomRefCpg += methyValue;

					} else if (pattern.equalsIgnoreCase("CH")) {
						sum.nCphTotal++;
						sum.sumCphTotal += methyValue;
						sum.nCphHom++;
						sum.sumCphHom += methyValue;
						sum.nCphHomRefCpg++;
						sum.sumCphHomRefCpg += methyValue;
					} else if (pattern.equalsIgnoreCase(".") || pattern.equalsIgnoreCase("C")) {

					} else {
						byte[] bases = pattern.getBytes();
						if (BaseUtils.basesAreEqual(bases[0], BaseUtilsMore.C) && BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[1], BaseUtilsMore.G)) {
							// if(BisSNPUtils.isHetCpg_at_G(value.vc)){

							sum.nCpgTotal++;
							sum.sumCpgTotal += methyValue;
							sum.nCpgHet++;
							sum.sumCpgHet += methyValue;
							sum.nCpgHetRefCpg++;
							sum.sumCpgHetRefCpg += methyValue;
							sum.nCpgHetAtG_RefCpg++;
							sum.sumCpgHetAtG_RefCpg += methyValue;

						}
					}

				} else if (BisSNPUtils.isRefCph(value.refContext)) { // ref cph
					String pattern = value.vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".");
					if (pattern.equalsIgnoreCase("CG")) {

						sum.nCpgTotal++;
						sum.sumCpgTotal += methyValue;
						sum.nCpgHom++;
						sum.sumCpgHom += methyValue;
						sum.nCpgHomRefCph++;
						sum.sumCpgHomRefCph += methyValue;

					} else if (pattern.equalsIgnoreCase("CH")) {
						sum.nCphTotal++;
						sum.sumCphTotal += methyValue;
						sum.nCphHom++;
						sum.sumCphHom += methyValue;
						sum.nCphHomRefCph++;
						sum.sumCphHomRefCph += methyValue;
					} else if (pattern.equalsIgnoreCase(".") || pattern.equalsIgnoreCase("C")) {

					} else {
						byte[] bases = pattern.getBytes();
						if (BaseUtils.basesAreEqual(bases[0], BaseUtilsMore.C) && BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[1], BaseUtilsMore.G)) {

							sum.nCpgTotal++;
							sum.sumCpgTotal += methyValue;
							sum.nCpgHet++;
							sum.sumCpgHet += methyValue;
							sum.nCpgHetRefCph++;
							sum.sumCpgHetRefCph += methyValue;
							sum.nCpgHetAtG_RefCph++;
							sum.sumCpgHetAtG_RefCph += methyValue;

						}
					}
				} else { // ref not C
					String pattern = value.vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".");
					if (pattern.equalsIgnoreCase("CG")) {
						sum.nCpgTotal++;
						sum.sumCpgTotal += methyValue;
						sum.nCpgHom++;
						sum.sumCpgHom += methyValue;
						sum.nCpgHomRefNotC++;
						sum.sumCpgHomRefNotC += methyValue;
					} else if (pattern.equalsIgnoreCase("CH")) {
						sum.nCphTotal++;
						sum.sumCphTotal += methyValue;
						sum.nCphHom++;
						sum.sumCphHom += methyValue;
						sum.nCphHomRefNotC++;
						sum.sumCphHomRefNotC += methyValue;
					} else if (pattern.equalsIgnoreCase(".") || pattern.equalsIgnoreCase("C")) {

					} else {
						byte[] bases = pattern.getBytes();
						if (BaseUtils.basesAreEqual(bases[0], BaseUtilsMore.C) && BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[1], BaseUtilsMore.G)) {
							sum.nCpgTotal++;
							sum.sumCpgTotal += methyValue;
							sum.nCpgHet++;
							sum.sumCpgHet += methyValue;
							sum.nCpgHetRefNotC++;
							sum.sumCpgHetRefNotC += methyValue;
							sum.nCpgHetAtG_RefNotC++;
							sum.sumCpgHetAtG_RefNotC += methyValue;
						}
					}
				}
			} else {// exclude C/T SNP at cytosine sites!!
				if (BisSNPUtils.isRefCpg(value.refContext)) {
					String pattern = value.vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".");
					if (pattern.contains("G") && !pattern.contains("Y")) {
						sum.nCpgTotal++;
						sum.sumCpgTotal += methyValue;
						sum.nCpgHet++;
						sum.sumCpgHet += methyValue;
						sum.nCpgHetRefCpg++;
						sum.sumCpgHetRefCpg += methyValue;
						sum.nCpgHetAtC_RefCpg++;
						sum.sumCpgHetAtC_RefCpg += methyValue;
					} else if (pattern.contains("H") && !pattern.contains("Y")) {
						sum.nCphTotal++;
						sum.sumCphTotal += methyValue;
						sum.nCphHet++;
						sum.sumCphHet += methyValue;
						sum.nCphHetRefCpg++;
						sum.sumCphHetRefCpg += methyValue;
						sum.nCphHetAtC_RefCpg++;
						sum.sumCphHetAtC_RefCpg += methyValue;
					}

				} else if (BisSNPUtils.isRefCph(value.refContext)) {
					String pattern = value.vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".");
					if (pattern.contains("G") && !pattern.contains("Y")) {
						sum.nCpgTotal++;
						sum.sumCpgTotal += methyValue;
						sum.nCpgHet++;
						sum.sumCpgHet += methyValue;
						sum.nCpgHetRefCph++;
						sum.sumCpgHetRefCph += methyValue;
						sum.nCpgHetAtC_RefCph++;
						sum.sumCpgHetAtC_RefCph += methyValue;
					} else if (pattern.contains("H") && !pattern.contains("Y")) {
						sum.nCphTotal++;
						sum.sumCphTotal += methyValue;
						sum.nCphHet++;
						sum.sumCphHet += methyValue;
						sum.nCphHetRefCph++;
						sum.sumCphHetRefCph += methyValue;
						sum.nCphHetAtC_RefCph++;
						sum.sumCphHetAtC_RefCph += methyValue;
					}
				} else {
					String pattern = value.vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".");
					if (pattern.contains("G") && !pattern.contains("Y")) {
						sum.nCpgTotal++;
						sum.sumCpgTotal += methyValue;
						sum.nCpgHet++;
						sum.sumCpgHet += methyValue;
						sum.nCpgHetRefNotC++;
						sum.sumCpgHetRefNotC += methyValue;
						sum.nCpgHetAtC_RefNotC++;
						sum.sumCpgHetAtC_RefNotC += methyValue;
					} else if (pattern.contains("H") && !pattern.contains("Y")) {
						sum.nCphTotal++;
						sum.sumCphTotal += methyValue;
						sum.nCphHet++;
						sum.sumCphHet += methyValue;
						sum.nCphHetRefNotC++;
						sum.sumCphHetRefNotC += methyValue;
						sum.nCphHetAtC_RefNotC++;
						sum.sumCphHetAtC_RefNotC += methyValue;
					}
				}
			}
		}

		return sum;
	}

	@Override
	public VCFCondition reduceInit() {
		VCFCondition vcfCondition = new VCFCondition();
		return vcfCondition;
	}

	@Override
	public VCFCondition treeReduce(VCFCondition lhs, VCFCondition rhs) {
		lhs.nLociVisited += rhs.nLociVisited;

		lhs.nCpgTotal += rhs.nCpgTotal;
		lhs.nCpgHom += rhs.nCpgHom;
		lhs.nCpgHomRefCpg += rhs.nCpgHomRefCpg;
		lhs.nCpgHomRefCph += rhs.nCpgHomRefCph;
		lhs.nCpgHomRefNotC += rhs.nCpgHomRefNotC;
		lhs.nCpgHet += rhs.nCpgHet;
		lhs.nCpgHetRefCpg += rhs.nCpgHetRefCpg;
		lhs.nCpgHetRefCph += rhs.nCpgHetRefCph;

		lhs.nCpgHetAtC_RefCpg += rhs.nCpgHetAtC_RefCpg;
		lhs.nCpgHetAtG_RefCpg += rhs.nCpgHetAtG_RefCpg;
		lhs.nCpgHetAtC_RefCph += rhs.nCpgHetAtC_RefCph;
		lhs.nCpgHetAtG_RefCph += rhs.nCpgHetAtG_RefCph;
		lhs.nCpgHetAtC_RefNotC += rhs.nCpgHetAtC_RefNotC;
		lhs.nCpgHetAtG_RefNotC += rhs.nCpgHetAtG_RefNotC;

		lhs.nCphTotal += rhs.nCphTotal;
		lhs.nCphHom += rhs.nCphHom;
		lhs.nCphHomRefCpg += rhs.nCphHomRefCpg;
		lhs.nCphHomRefCph += rhs.nCphHomRefCph;
		lhs.nCphHomRefNotC += rhs.nCphHomRefNotC;
		lhs.nCphHet += rhs.nCphHet;
		lhs.nCphHetRefCpg += rhs.nCphHetRefCpg;
		lhs.nCphHetRefCph += rhs.nCphHetRefCph;

		lhs.nCphHetAtC_RefCpg += rhs.nCphHetAtC_RefCpg;
		lhs.nCphHetAtC_RefCph += rhs.nCphHetAtC_RefCph;
		lhs.nCphHetAtC_RefNotC += rhs.nCphHetAtC_RefNotC;
		lhs.nCphHetAtC_RefNotC += rhs.nCphHetAtC_RefNotC;

		// lhs.sumCytosineTotal += rhs.sumCytosineTotal;

		lhs.sumCpgTotal += rhs.sumCpgTotal;
		lhs.sumCpgHom += rhs.sumCpgHom;
		lhs.sumCpgHomRefCpg += rhs.sumCpgHomRefCpg;
		lhs.sumCpgHomRefCph += rhs.sumCpgHomRefCph;
		lhs.sumCpgHomRefNotC += rhs.sumCpgHomRefNotC;
		lhs.sumCpgHet += rhs.sumCpgHet;
		lhs.sumCpgHetRefCpg += rhs.sumCpgHetRefCpg;
		lhs.sumCpgHetRefCph += rhs.sumCpgHetRefCph;

		lhs.sumCpgHetAtC_RefCpg += rhs.sumCpgHetAtC_RefCpg;
		lhs.sumCpgHetAtG_RefCpg += rhs.sumCpgHetAtG_RefCpg;
		lhs.sumCpgHetAtC_RefCph += rhs.sumCpgHetAtC_RefCph;
		lhs.sumCpgHetAtG_RefCph += rhs.sumCpgHetAtG_RefCph;
		lhs.sumCpgHetAtC_RefNotC += rhs.sumCpgHetAtC_RefNotC;
		lhs.sumCpgHetAtG_RefNotC += rhs.sumCpgHetAtG_RefNotC;

		lhs.sumCphTotal += rhs.sumCphTotal;
		lhs.sumCphHom += rhs.sumCphHom;
		lhs.sumCphHomRefCpg += rhs.sumCphHomRefCpg;
		lhs.sumCphHomRefCph += rhs.sumCphHomRefCph;
		lhs.sumCphHomRefNotC += rhs.sumCphHomRefNotC;
		lhs.sumCphHet += rhs.sumCphHet;
		lhs.sumCphHetRefCpg += rhs.sumCphHetRefCpg;
		lhs.sumCphHetRefCph += rhs.sumCphHetRefCph;

		lhs.sumCphHetAtC_RefCpg += rhs.sumCphHetAtC_RefCpg;
		lhs.sumCphHetAtC_RefCph += rhs.sumCphHetAtC_RefCph;
		lhs.sumCphHetAtC_RefNotC += rhs.sumCphHetAtC_RefNotC;
		lhs.sumCphHetAtC_RefNotC += rhs.sumCphHetAtC_RefNotC;

		return lhs;
	}

	private boolean passFilter(VariantContext vc) {
		if (vc.getPhredScaledQual() <= qual)
			return false;
		if (vc.hasAttribute(VCFConstants.DEPTH_KEY) && vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0) > maxCov) {
			return false;
		}
		if (vc.hasAttribute(VCFConstants.STRAND_BIAS_KEY) && vc.getAttributeAsDouble(VCFConstants.STRAND_BIAS_KEY, 0.0) > sb) {
			return false;
		}
		if (vc.hasAttribute(VCFConstants.MAPPING_QUALITY_ZERO_KEY)) {
			if (vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0) > 40 && vc.getAttributeAsDouble(VCFConstants.MAPPING_QUALITY_ZERO_KEY, 0) / vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0) > mq0)
				return false;
		}
		if (vc.hasAttribute(BisulfiteVCFConstants.QUAL_BY_DEPTH) && vc.getAttributeAsDouble(BisulfiteVCFConstants.QUAL_BY_DEPTH, 10000.0) < qd) {
			return false;
		}
		if (vc.getGenotype(0).hasAnyAttribute(BisulfiteVCFConstants.C_STATUS)) {
			String cytosineStatus = vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.C_STATUS, ".");
			if (!cytosineStatus.equalsIgnoreCase(".")) {
				String[] tmp = cytosineStatus.split(",");
				int numC = Integer.parseInt(tmp[0]);
				int numT = Integer.parseInt(tmp[1]);

				if (numC + numT < minCT) {
					return false;
				}
			}

		}
		return true;
	}

	private boolean passSNPclusterFilter(ReferenceContext ref) {
		GenomeLoc searchLoc = getToolkit().getGenomeLocParser().createGenomeLoc(ref.getLocus().getLocation().getContig(), ref.getLocus().getLocation().getStart() - windSizeForSNPfilter,
				ref.getLocus().getLocation().getStart() + windSizeForSNPfilter);
		LocationAwareSeekableRODIterator locRodIt = rodIt.seek(searchLoc);
		if (locRodIt.hasNext()) {
			RODRecordList rodList = locRodIt.seekForward(searchLoc);
			if (rodList != null && rodList.size() >= minSNPinWind) {
				rodIt.close(locRodIt);
				return false;
			}
		}
		rodIt.close(locRodIt);
		return true;
	}

	// to test position next to cytosine, but inside cytosine pattern.
	private boolean passSNPfilter(ReferenceContext ref, int pos) {
		GenomeLoc searchLoc = getToolkit().getGenomeLocParser().createGenomeLoc(ref.getLocus().getLocation().getContig(), ref.getLocus().getLocation().getStart() + pos);
		LocationAwareSeekableRODIterator locRodIt = rodIt.seek(searchLoc);
		if (locRodIt.hasNext()) {
			RODRecordList rodList = locRodIt.seekForward(searchLoc);
			if (rodList == null) {
				rodIt.close(locRodIt);
				return false;
			} else {
				VariantContext snp = (VariantContext) rodList.get(0).getUnderlyingObject();
				if (passFilter(snp)) {
					rodIt.close(locRodIt);
					return true;
				}
			}
		}
		rodIt.close(locRodIt);
		return false;
	}

	public static class VCFCondition {
		long nCpgHet = 0;

		long nCpgHetAtC_RefCpg = 0;

		long nCpgHetAtC_RefCph = 0;

		long nCpgHetAtC_RefNotC = 0;

		long nCpgHetAtG_RefCpg = 0;

		long nCpgHetAtG_RefCph = 0;

		long nCpgHetAtG_RefNotC = 0;

		long nCpgHetRefCpg = 0;

		long nCpgHetRefCph = 0;

		long nCpgHetRefNotC = 0;

		long nCpgHom = 0;

		long nCpgHomRefCpg = 0;

		long nCpgHomRefCph = 0;

		long nCpgHomRefNotC = 0;

		// double sumCytosineTotal = 0;
		long nCpgTotal = 0;

		long nCphHet = 0;

		long nCphHetAtC_RefCpg = 0;

		long nCphHetAtC_RefCph = 0;

		long nCphHetAtC_RefNotC = 0;

		long nCphHetRefCpg = 0;

		long nCphHetRefCph = 0;

		long nCphHetRefNotC = 0;

		long nCphHom = 0;

		long nCphHomRefCpg = 0;

		long nCphHomRefCph = 0;

		long nCphHomRefNotC = 0;

		long nCphTotal = 0;

		long nLociVisited = 0;

		double sumCpgHet = 0;

		double sumCpgHetAtC_RefCpg = 0;

		double sumCpgHetAtC_RefCph = 0;

		double sumCpgHetAtC_RefNotC = 0;

		double sumCpgHetAtG_RefCpg = 0;

		double sumCpgHetAtG_RefCph = 0;

		double sumCpgHetAtG_RefNotC = 0;

		double sumCpgHetRefCpg = 0;

		double sumCpgHetRefCph = 0;

		double sumCpgHetRefNotC = 0;

		double sumCpgHom = 0;

		double sumCpgHomRefCpg = 0;

		double sumCpgHomRefCph = 0;

		double sumCpgHomRefNotC = 0;

		double sumCpgTotal = 0;

		double sumCphHet = 0;

		double sumCphHetAtC_RefCpg = 0;

		double sumCphHetAtC_RefCph = 0;

		double sumCphHetAtC_RefNotC = 0;

		double sumCphHetRefCpg = 0;

		double sumCphHetRefCph = 0;

		double sumCphHetRefNotC = 0;

		double sumCphHom = 0;

		double sumCphHomRefCpg = 0;

		double sumCphHomRefCph = 0;

		double sumCphHomRefNotC = 0;

		double sumCphTotal = 0;

	}

	public class VCFplusRef {
		public boolean passFilter;
		public ReferenceContext refContext;
		public VariantContext vc;

		public VCFplusRef(ReferenceContext ref, VariantContext vc_eval, boolean pass) {
			refContext = ref;
			vc = vc_eval;
			passFilter = pass;
		}

	}

}
