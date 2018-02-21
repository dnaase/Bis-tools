package main.java.edu.usc.epigenome.uecgatk.bissnp;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.FileOutputStream;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.HashMap;
import java.io.OutputStream;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.filters.DuplicateReadFilter;
import org.broadinstitute.gatk.engine.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.gatk.engine.filters.UnmappedReadFilter;
import org.broadinstitute.gatk.engine.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.gatk.engine.filters.MappingQualityZeroFilter;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.By;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.Downsample;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.engine.walkers.ReadFilters;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.samtools.util.BlockCompressedOutputStream;

import main.java.edu.usc.epigenome.uecgatk.bissnp.writer.FormatWriterBase;
import main.java.edu.usc.epigenome.uecgatk.bissnp.writer.SortingTcgaVCFWriter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.writer.TcgaVCFWriter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.writer.cpgReads;
import main.java.edu.usc.epigenome.uecgatk.bissnp.writer.cpgReadsWriterImp;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.BisulfiteIncompleteConvReadsFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.NotProperPairedReadFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.InvertedDupsReadFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.MappingQualityFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.BisulfiteFivePrimeConvReadsFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.BisulfiteMismatchReadsFilter;


/*
 * Bis-SNP/BisSNP: It is a genotyping and methylation calling in bisulfite treated massively
 * parallel sequencing (Bisulfite-seq and NOMe-seq) on Illumina platform Copyright (C) <2011>
 * <Yaping Liu: lyping1986@gmail.com>
 * 
 * This program is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program. If
 * not, see <http://www.gnu.org/licenses/>.
 */

//TO DO: add confident interval into VCF outpout for methylation level: ProportionCI(long x, long n, double confidenceCoeff) function in jsc.onesample.ProportionCI could provide Clopper-Pearson "exact" confidence interval 

/**
 * A Bisulfite genotyper. Works for single-sample data right now.
 */

//@BisBAQMode(QualityMode = BisBAQ.QualityMode.ADD_TAG, ApplicationTime = BisBAQ.ApplicationTime.ON_INPUT)
@Reference(window = @Window(start = -500, stop = 500))
@By(DataSource.READS)
@ReadFilters({UnmappedReadFilter.class, DuplicateReadFilter.class, NotPrimaryAlignmentFilter.class, FailsVendorQualityCheckFilter.class, MappingQualityFilter.class, NotProperPairedReadFilter.class, InvertedDupsReadFilter.class, BisulfiteMismatchReadsFilter.class, BisulfiteIncompleteConvReadsFilter.class, BisulfiteFivePrimeConvReadsFilter.class})
@Downsample(by = DownsampleType.NONE)
public class BisulfiteGenotyper extends LocusWalker<BisulfiteVariantCallContext, BisulfiteGenotyper.ContextCondition> implements
		TreeReducible<BisulfiteGenotyper.ContextCondition> {

	@ArgumentCollection
	private static BisulfiteArgumentCollection BAC = new BisulfiteArgumentCollection();
	

	private static int MAXIMUM_CACHE_FOR_OUTPUT_VCF;

	private static Set<String> samples = null;
	@Output(fullName = "file_name_output_cpg_reads_detail", shortName = "cpgreads", doc = "output Haplotype CpG reads bed file that contain each CpG's position, methylation and reads name info ", required = false)
	public String cpgreads = null;
	@Output(fullName = "file_name_output_gch_reads_detail", shortName = "gchreads", doc = "output Haplotype GCH reads bed file that contain each GCH position, methylation and reads name info ", required = false)
	public String gchreads = null;

//	@Output(fullName = "file_name_output_reads_after_downsampling", shortName = "reads", doc = "output Bam file's name that after downsapling, for performance test only", required = false)
//	public String reads = null;

	@Output(fullName = "vcf_file_name_1", shortName = "vfn1", doc = "output Vcf file, when used for [DEFAULT_FOR_TCGA] output mode, it is used to store all CpG sites. While the original vcf file is to store all CpG sites", required = true)
	public String vfn1 = null;

	@Output(fullName = "vcf_file_name_2", shortName = "vfn2", doc = "output Vcf file 2, only used for [DEFAULT_FOR_TCGA] output mode, it is used to store all SNP sites. While the original vcf file is to store all CpG sites", required = false)
	public String vfn2 = null;
	
	@Output(fullName = "coverage_metric", shortName = "coverage", doc = "output coverage metric files, 1st column is coverage, 2nd column is '1' when give -cgi option and inside this option's file region, otherwise they are all '0' ??should we seperate it to another walker?", required = false)
	public String coverage = null;
	//To Do!!
	@Input(fullName = "cgi_file", shortName = "cgi", doc = "Give the CGI bed file for the coverage estimation", required = false)
	public RodBinding<BEDFeature> cgi = null;

	private TcgaVCFWriter additionalWriterForDefaultTcgaMode = null;

	private String argCommandline = "";

	private BisulfiteGenotyperEngine BG_engine = null;

	

	private SortingTcgaVCFWriter multiAdditionalWriterForDefaultTcgaMode = null;


	private SortingTcgaVCFWriter multiThreadWriter = null;

	// cpgReads output
	private cpgReadsWriterImp readsWriter = null;
	private cpgReadsWriterImp gchReadsWriter = null;

	// only works with single core right now
	

	private TcgaVCFWriter verboseWriter = null;

	private TcgaVCFWriter writer = null;
	
	private PrintWriter covWriter = null;

	public static BisulfiteArgumentCollection getBAC() {
		return BAC;
	}

	/**
	 * Initialize the samples, output, and genotype calculation model
	 * 
	 **/
	@Override
	public void initialize() {
		samples = new TreeSet<String>();
		samples = ReadUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
		logger.info("samples provided: " + samples.toString());
		if (samples.isEmpty()) {	
				throw new UserException("No sample name provided, program will automately provide the bam file header: " + samples.toString());

		}
		MAXIMUM_CACHE_FOR_OUTPUT_VCF = BAC.vcfCache;
		// initiate BisulfiteGenotyperEngine

		SAMSequenceDictionary refDict = getToolkit().getMasterSequenceDictionary();
		// initialize the header
		argCommandline = BisSNP.getBisSNPArgumantsInput();
		
		initiateVCFInDifferentOutmode(refDict);
		// in the first iteration, initiate CytosineTypeStatus
		
		BAC.makeCytosine();
		
		if(coverage != null){
			try {
				covWriter = new PrintWriter(new File(coverage));
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		 
	}

	/**
	 * Compute at a given locus.
	 * 
	 * @param tracker
	 *            the meta data tracker
	 * @param refContext
	 *            the reference base
	 * @param rawContext
	 *            contextual information around the locus
	 * @return the VariantCallContext object
	 */
	@Override
	public BisulfiteVariantCallContext map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
		
		if (rawContext.size() > BAC.toCoverage && !rawContext.getContig().equalsIgnoreCase(BisulfiteSAMConstants.CONT_CHRM_UCSC) && !rawContext.getContig().equalsIgnoreCase(BisulfiteSAMConstants.CONT_CHRM_NCBI)){ // get rid of region with // abnormal read coverage
			if(BAC.includeAllInRef){
				return new BisulfiteVariantCallContext(null, null, rawContext, refContext);
			}
			else{
				return null;
			}
			
		}
		BG_engine = new BisulfiteGenotyperEngine(tracker, refContext, rawContext, BAC, getToolkit());

	//	if (BAC.orad) {
	//		downsampleBam.downsamplingBamFile(rawContext);
	//	}

		// calculation LikelihoodsAndGenotypes for this loci
		return BG_engine.getBisulfiteVariantCallContext();
	}

	@Override
	public void onTraversalDone(ContextCondition sum) {
		String fn = vfn1 + ".MethySummarizeList.txt";
		PrintWriter outWriter = null;
		String outLine = null;
		try {
			outWriter = new PrintWriter(new File(fn));

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		logger.info(String.format("Visited bases                                %d", sum.nBasesVisited));
		logger.info(String.format("Callable bases                               %d", sum.nBasesCallable));
		logger.info(String.format("Confidently called bases                     %d", sum.nBasesCalledConfidently));
		logger.info(String.format("Visited bases in Cpg island                     %d", sum.nBasesInCgi));
		logger.info(String.format("Callable bases in Cpg island                     %d", sum.nBasesCallableInCgi));
		logger.info(String.format("Confidently called bases in Cpg island                     %d", sum.nBasesConfidentlyInCgi));
		logger.info(String.format("%% callable bases of all loci                 %3.3f", sum.percentCallableOfAll()));
		logger.info(String.format("%% confidently called bases of all loci       %3.3f", sum.percentCalledOfAll()));
		logger.info(String.format("%% confidently called bases of callable loci  %3.3f", sum.percentCalledOfCallable()));
		logger.info(String.format("Actual calls made                            %d", sum.nCallsMade));
		logger.info(String.format("Average good reads coverage in all visited loci                            %3.1f", sum.meanCoverageOfAll()));
		logger.info(String.format("Average good reads coverage in callable position                            %3.1f", sum.meanCoverageOfCallable()));
		logger.info(String.format("Average good reads coverage in all position of Cpg island                            %3.1f", sum.meanCoverageOfCpgIsland()));
		logger.info(String.format("Average good reads coverage in callable position of Cpg island                            %3.1f", sum.meanCoverageOfCallableInCpgIsland()));
		outLine = String.format("BisSNP version                                %s", BisSNP.getBisSNPVersionNumber()) + "\n";
		outLine = outLine +  String.format("Visited bases                                %d", sum.nBasesVisited) + "\n";
		outLine = outLine + String.format("Callable bases                               %d", sum.nBasesCallable) + "\n";
		outLine = outLine + String.format("Confidently called bases                     %d", sum.nBasesCalledConfidently) + "\n";
		outLine = outLine + String.format("Visited bases in Cpg island                               %d", sum.nBasesInCgi) + "\n";
		outLine = outLine + String.format("Callable bases in Cpg island                               %d", sum.nBasesCallableInCgi) + "\n";
		outLine = outLine + String.format("Confidently called bases in Cpg island                     %d", sum.nBasesConfidentlyInCgi) + "\n";
		outLine = outLine + String.format("%% callable bases of all loci                 %3.3f", sum.percentCallableOfAll()) + "\n";
		outLine = outLine + String.format("%% confidently called bases of all loci       %3.3f", sum.percentCalledOfAll()) + "\n";
		outLine = outLine + String.format("%% confidently called bases of callable loci  %3.3f", sum.percentCalledOfCallable()) + "\n";
		outLine = outLine + String.format("Actual calls made                            %d", sum.nCallsMade) + "\n";
		outLine = outLine
				+ String.format("Average good reads coverage in all visited loci                            %3.1f", sum.meanCoverageOfAll()) + "\n";
		outLine = outLine
				+ String.format("Average good reads coverage in callable position                            %3.1f", sum.meanCoverageOfCallable())
				+ "\n";
		outLine = outLine
				+ String.format("Average good reads coverage in all position of Cpg island                           %3.1f", sum.meanCoverageOfCpgIsland())
				+ "\n";
		outLine = outLine
				+ String.format("Average good reads coverage in callable position of Cpg island                           %3.1f", sum.meanCoverageOfCallableInCpgIsland())
				+ "\n";
		if (!sum.cytosineMethySummary.isEmpty()) {
			// logger.info(String.format("C/T heterozygous at cytosine loci is not included into summary statistics, they are not in cytosine.vcf but in snp.vcf file"));
			logger.info(String.format("-------------Methylation summary(confidently called) in total "));
			outLine = outLine + String.format("##Methylation summary in total:") + "\n";
			for (String key : sum.cytosineMethySummary.keySet()) {
				Pair<Integer, Double> methyValue = sum.cytosineMethySummary.get(key);
				logger.info(String.format("%% average methylation level of confidently called %s loci in total:       %3.3f", key, 100 * methyValue.getSecond()
						/ methyValue.getFirst()));
				logger.info(String.format(" Number of confidently called %s loci in total:       %d", key, methyValue.getFirst() / samples.size()));
				outLine = outLine + key + ":\t" + methyValue.getFirst() / samples.size() + "\t"
						+ String.format("%3.3f", 100 * methyValue.getSecond() / methyValue.getFirst()) + "%\n";
			}
		}

		if (!sum.cytosineMethySummaryByReadsGroup.isEmpty()) {

			for (String sampleKey : sum.cytosineMethySummaryByReadsGroup.keySet()) {
			
				logger.info(String.format("---------- Methylation summary(confidently called) in Read Group:        %s", sampleKey));
				outLine = outLine + String.format("##Methylation summary in Read Group:\t%s", sampleKey) + "\n";
				TreeMap<String, Pair<Integer, Double>> formattedMap = sorHashMapByKeyValue(sum.cytosineMethySummaryByReadsGroup.get(sampleKey));
				for (String key : formattedMap.keySet()) {
					Pair<Integer, Double> methyValue = formattedMap.get(key);
					logger.info(String.format("%% methylation level of confidently called %s loci in total:       %3.3f", key,
							100 * methyValue.getSecond() / methyValue.getFirst()));
					logger.info(String.format(" Number of confidently called %s loci in total:       %d", key, methyValue.getFirst()));
					outLine = outLine + key + ":\t" + methyValue.getFirst() + "\t"
							+ String.format("%3.3f", 100 * methyValue.getSecond() / methyValue.getFirst()) + "%\n";
				}
			}

		}

	//	if (BAC.orad) {
	//		samWriter.close();
	//	}

		if (getToolkit().getArguments().numberOfCPUThreadsPerDataThread > 1) {
			multiThreadWriter.writerFlush();
			multiThreadWriter.close();

			if (BAC.OutputMode == BisulfiteEnums.OUTPUT_MODE.DEFAULT_FOR_TCGA || BAC.OutputMode == BisulfiteEnums.OUTPUT_MODE.EMIT_VARIANT_AND_CYTOSINES ||  BAC.OutputMode == BisulfiteEnums.OUTPUT_MODE.NOMESEQ_MODE) {
				multiAdditionalWriterForDefaultTcgaMode.writerFlush();
				multiAdditionalWriterForDefaultTcgaMode.close();
			}
		}
		if (BAC.OutputMode == BisulfiteEnums.OUTPUT_MODE.DEFAULT_FOR_TCGA || BAC.OutputMode == BisulfiteEnums.OUTPUT_MODE.EMIT_VARIANT_AND_CYTOSINES ||  BAC.OutputMode == BisulfiteEnums.OUTPUT_MODE.NOMESEQ_MODE) {
			additionalWriterForDefaultTcgaMode.writeFlush();
			additionalWriterForDefaultTcgaMode.close();
		}

		writer.close();

		if (cpgreads != null) {
				readsWriter.close();
				if (BAC.sequencingMode == BisulfiteEnums.MethylSNPModel.GM && gchreads != null) {
					gchReadsWriter.close();
				} 
				
		}
		outWriter.println(outLine);
		outWriter.close();
		if(coverage != null){
				covWriter.close();		
		}

	}

	/**
	 * calculate statistics number in each reduce steps.
	 */
	@Override
	public ContextCondition reduce(BisulfiteVariantCallContext value, ContextCondition sum) {
		// the base vistited
		sum.nBasesVisited++;

		// not call the locus when there is no coverage
		if (value == null || value.getVariantContext() == null){
			
				return sum;
			
		}
		else{
			sum.nBasesInCgi += value.isCgi ? 1 : 0;
		}
			

		sum.nBasesCallable++;
		sum.nBasesCallableInCgi += value.isCgi ? 1 : 0;

		// the base was confidently callable
		sum.nBasesCalledConfidently += value.confidentlyCalled ? 1 : 0;
		sum.nBasesConfidentlyInCgi += value.confidentlyCalled ? (value.isCgi ? 1 : 0) : 0;

		// the base was confidently callable
		sum.nReadsCallable += value.getSummaryAcrossRG().getGoodReadsCoverage();
		sum.nReadsCgi += value.isCgi ? value.getSummaryAcrossRG().getGoodReadsCoverage() : 0;

		// can't make a confident variant call here
		if (!value.shouldEmit)
			return sum;

		// cytosines provided summary in each readGroup

		//for (String sampleKey : value.getBisulfiteContextsGenotypeLikelihoods().keySet()) { // seperated
		for(Entry<String, BisulfiteContextsGenotypeLikelihoods> valueEntrySet : value.getBisulfiteContextsGenotypeLikelihoods().entrySet()){																					// by
			String sampleKey = valueEntrySet.getKey();																			// ReadGroup
			BisulfiteContextsGenotypeLikelihoods bcglTmp = valueEntrySet.getValue();
			HashMap<String, CytosineParameters> cytosineParameters = bcglTmp.getCytosineParameters();
			HashMap<String, Pair<Integer, Double>> cytosineMethySummaryInEachReadGroup = sum.cytosineMethySummaryByReadsGroup.get(sampleKey);
			//for (String cytosineKey : cytosineParameters.keySet()) {
			for (Entry<String, CytosineParameters> cytosineParametersEntrySet : cytosineParameters.entrySet()) {
				String cytosineKey = cytosineParametersEntrySet.getKey();
				if (cytosineParameters.get(cytosineKey).isCytosinePattern || cytosineParameters.get(cytosineKey).isReferenceCytosinePattern) {
					if (cytosineParameters.get(cytosineKey).isCytosinePattern) {
						Integer cytosineNum = cytosineMethySummaryInEachReadGroup.get(cytosineKey).getFirst() + 1;
						Double methyValue = cytosineMethySummaryInEachReadGroup.get(cytosineKey).getSecond();
						if (!Double.isNaN(bcglTmp.getMethylationLevel())) {
							methyValue += bcglTmp.getMethylationLevel();
						}
						else{
							cytosineNum--; // identified CpGs without methylation information are not take into statistics account..
						}

						Pair<Integer, Double> sumValue = cytosineMethySummaryInEachReadGroup.get(cytosineKey);
						sumValue.set(cytosineNum, methyValue);
						cytosineMethySummaryInEachReadGroup.put(cytosineKey, sumValue);
					}

				}
			}
			sum.cytosineMethySummaryByReadsGroup.put(sampleKey, cytosineMethySummaryInEachReadGroup);
		}

		// cytosines provided summary in total

		for (Entry<String, CytosineParameters> cytosineEntrySet : BAC.cytosineDefined.getContextDefined().entrySet()) {
			String cytosineKey = cytosineEntrySet.getKey();
			Integer cytosineNumTotal = 0;
			Double methyValueTotal = 0.0;
			for (HashMap<String, Pair<Integer, Double>> cytosineMethyEachRG : sum.cytosineMethySummaryByReadsGroup.values()) {
				Pair<Integer, Double> tmp = cytosineMethyEachRG.get(cytosineKey);
				if(tmp != null){
					cytosineNumTotal += tmp.getFirst();
					methyValueTotal += tmp.getSecond();
				}
				/*
				if (cytosineMethyEachRG.containsKey(cytosineKey)) {
					cytosineNumTotal += cytosineMethyEachRG.get(cytosineKey).getFirst();
					methyValueTotal += cytosineMethyEachRG.get(cytosineKey).getSecond();
				}
				*/

			}
			Pair<Integer, Double> sumValueInTotal = sum.cytosineMethySummary.get(cytosineKey);
			sumValueInTotal.set(cytosineNumTotal, methyValueTotal);

			sum.cytosineMethySummary.put(cytosineKey, sumValueInTotal);

		}

		try {
			// actually making a call
			sum.nCallsMade++;
			
			readsReport(value);

			outputVCFInDifferentOutmode(value);

		} catch (IllegalArgumentException e) {
			throw new IllegalArgumentException(e.getMessage()
					+ "; this is often caused by using the --assume_single_sample_reads argument with the wrong sample name");
		}

		return sum;
	}

	/**
	 * Initiate statistics object.
	 */
	@Override
	public ContextCondition reduceInit() {
		ContextCondition initiated = new ContextCondition();
		initiated.makeCytosineMap();
		return initiated;

	}

	@Override
	public ContextCondition treeReduce(ContextCondition lhs, ContextCondition rhs) {
		lhs.nBasesCallable += rhs.nBasesCallable;
		lhs.nBasesCalledConfidently += rhs.nBasesCalledConfidently;
		lhs.nBasesVisited += rhs.nBasesVisited;
		lhs.nCallsMade += rhs.nCallsMade;
		lhs.nReadsCallable += rhs.nReadsCallable;
		lhs.nBasesInCgi += rhs.nBasesInCgi;
		lhs.nBasesCallableInCgi += rhs.nBasesCallableInCgi;
		lhs.nBasesConfidentlyInCgi += rhs.nBasesConfidentlyInCgi;
		lhs.nReadsCgi += rhs.nReadsCgi;
		if (!rhs.cytosineMethySummary.isEmpty()) {
			// System.err.println( "rudcue\t");
			//for (String key : rhs.cytosineMethySummary.keySet()) {
			for (Entry<String, Pair<Integer, Double>>  cytosineMethySummaryEntrySet: rhs.cytosineMethySummary.entrySet()) {
				String key = cytosineMethySummaryEntrySet.getKey();
				Pair<Integer, Double> rhsValue = cytosineMethySummaryEntrySet.getValue();
				// System.err.println(key + "rudcue\t" + rhsValue.toString());
				Pair<Integer, Double> lhsValue = lhs.cytosineMethySummary.get(key);
				if(lhsValue != null){
					lhsValue.first += rhsValue.getFirst();
					lhsValue.second += rhsValue.getSecond();
					lhs.cytosineMethySummary.put(key, lhsValue);
				}
				/*
				if (lhs.cytosineMethySummary.containsKey(key)) {
					lhsValue = lhs.cytosineMethySummary.get(key);
					lhsValue.first += rhsValue.getFirst();
					lhsValue.second += rhsValue.getSecond();
					lhs.cytosineMethySummary.put(key, lhsValue);

				}
				*/

			}
		}

		if (!rhs.cytosineMethySummaryByReadsGroup.isEmpty()) {
			//for (String sampleKey : rhs.cytosineMethySummaryByReadsGroup.keySet()) {
			for (Entry<String, HashMap<String, Pair<Integer, Double>>> cytosineMethySummaryByReadsGroupEntrySet: rhs.cytosineMethySummaryByReadsGroup.entrySet()) {
				String sampleKey = cytosineMethySummaryByReadsGroupEntrySet.getKey();
				//if (lhs.cytosineMethySummaryByReadsGroup.containsKey(sampleKey)) {
					HashMap<String, Pair<Integer, Double>> cytosineMethySummaryLhs = lhs.cytosineMethySummaryByReadsGroup.get(sampleKey);
					if(cytosineMethySummaryLhs != null){
						//for (String cytosineKey : rhs.cytosineMethySummaryByReadsGroup.get(sampleKey).keySet()) {
						for ( Entry<String, Pair<Integer, Double>>  cytosineMethySummaryEntrySet:  rhs.cytosineMethySummaryByReadsGroup.get(sampleKey).entrySet()) {
							String cytosineKey = cytosineMethySummaryEntrySet.getKey();
							Pair<Integer, Double> rhsValue = cytosineMethySummaryEntrySet.getValue();
							Pair<Integer, Double> lhsValue = cytosineMethySummaryLhs.get(cytosineKey);
							if(lhsValue != null){
								lhsValue.first += rhsValue.getFirst();
								lhsValue.second += rhsValue.getSecond();
								cytosineMethySummaryLhs.put(cytosineKey, lhsValue);
							}
							/*
							if (cytosineMethySummaryLhs.containsKey(cytosineKey)) {
								lhsValue = cytosineMethySummaryLhs.get(cytosineKey);
								lhsValue.first += rhsValue.getFirst();
								lhsValue.second += rhsValue.getSecond();
								cytosineMethySummaryLhs.put(cytosineKey, lhsValue);
							}
							*/
						}
						lhs.cytosineMethySummaryByReadsGroup.put(sampleKey, cytosineMethySummaryLhs);
					}
					
				//}
			}
		}

		return lhs;
	}

	/**
	 * get VCF header for the output VCF file
	 * 
	 **/
	private Set<VCFHeaderLine> getHeaderInfo() {
		Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();

		if (!BAC.NO_SLOD)
			headerInfo.add(new VCFInfoHeaderLine(VCFConstants.STRAND_BIAS_KEY, 1, VCFHeaderLineType.Float, "Strand Bias"));

		// Bisulfite-seq own INFO column
		headerInfo.add(new VCFInfoHeaderLine(BisulfiteVCFConstants.CYTOSINE_TYPE, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,
				"Cytosine Context. e.g. 'CG' means homozygous CG sites, while heterozygous CG sites will be output as IUPAC code"));

		headerInfo.add(new VCFInfoHeaderLine(BisulfiteVCFConstants.C_STRAND_KEY, 1, VCFHeaderLineType.Character,
				"Strand of cytosine relative to reference genome. Does not apply to non-cytosine positions"));

		headerInfo.add(new VCFInfoHeaderLine(BisulfiteVCFConstants.NUM_OF_SAMPLES, 1, VCFHeaderLineType.Integer, "Number of Samples With Data"));
		headerInfo.add(new VCFInfoHeaderLine(BisulfiteVCFConstants.REF_C_PATTERN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,
				"Reference Cytosine Context. if all cytosine context found are in reference genome, REF=0; otherwise, e.g. it will be REF=CG for position not confidently called as CG or real context is CH."));

		// General VCF INFO column
		headerInfo.add(new VCFInfoHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer,
				"Total Depth, not filtered by mapping quality and base quality score criteria yet"));
		headerInfo.add(new VCFInfoHeaderLine(BisulfiteVCFConstants.QUAL_BY_DEPTH, 1, VCFHeaderLineType.Float, "Variant Confidence/Quality by Depth"));
		headerInfo
				.add(new VCFInfoHeaderLine(VCFConstants.MAPPING_QUALITY_ZERO_KEY, 1, VCFHeaderLineType.Integer, "Total Mapping Quality Zero Reads"));
		headerInfo.add(new VCFInfoHeaderLine(VCFConstants.HAPLOTYPE_QUALITY_KEY, 1, VCFHeaderLineType.Float,
				"Consistency of the site with at most two segregating haplotypes"));
		

		// check in dbSNP or not
		if (BAC.dbsnp.isBound())
			headerInfo.add(new VCFInfoHeaderLine(VCFConstants.DBSNP_KEY, 0, VCFHeaderLineType.Flag, "dbSNP Membership"));

		// FORMAT fields

		headerInfo.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));
		headerInfo.add(new VCFFormatHeaderLine(BisulfiteVCFConstants.NUMBER_OF_C_KEY, 1, VCFHeaderLineType.Integer,
				"Number of Unconverted Cytosines in this position(filtered by minConv)"));
		headerInfo.add(new VCFFormatHeaderLine(BisulfiteVCFConstants.NUMBER_OF_T_KEY, 1, VCFHeaderLineType.Integer,
				"Number of Converted Cytosines in this position(filtered by minConv)"));
		headerInfo
				.add(new VCFFormatHeaderLine(BisulfiteVCFConstants.READS_SUPPORT_ALT, 4, VCFHeaderLineType.Integer,
						"Reads supporting ALT, only keep good base. Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles"));
		headerInfo.add(new VCFFormatHeaderLine(BisulfiteVCFConstants.AVE_BASE_QUALITY_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float,
				"Average base quality for reads supporting alleles. For each allele, in the same order as listed"));
		headerInfo.add(new VCFFormatHeaderLine(BisulfiteVCFConstants.SOMATIC_STAT_VAR, 1, VCFHeaderLineType.Integer,
				"Somatic status of the variant. 1) wildtype; 2) germline, somatic; 3) LOH; 4) post-transcriptional modification; 5) unknown"));
		headerInfo.add(new VCFFormatHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer,
				"Approximate read depth, not filtered by mapping quality and base quality score criteria yet"));

		headerInfo.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY, 1, VCFHeaderLineType.Float, "Genotype Quality"));

		headerInfo.add(new VCFFormatHeaderLine(BisulfiteVCFConstants.BEST_C_PATTERN, 1, VCFHeaderLineType.String, "Best Cytosine Context"));
		headerInfo.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_POSTERIORS_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Integer,
				"Normalized, Phred-scaled posteriors for genotypes as defined in the VCF specification"));
		headerInfo
				.add(new VCFFormatHeaderLine(
						BisulfiteVCFConstants.C_STATUS,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.Integer,
						"Bisulfite read counts(not filtered by minConv): 1) number of C in cytosine strand, 2) number of T in cytosine strand, 3) number of A/G/N in cytosine strand,"
								+ " 4) number of G in guanine strand, 5) number of A in guanine strand, 6) number of C/T/N in guanine strand."));

		// FILTER fields
		if (BAC.STANDARD_CONFIDENCE_FOR_EMITTING < BAC.STANDARD_CONFIDENCE_FOR_CALLING)
			headerInfo.add(new VCFFilterHeaderLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME,
					"Low quality. When QUAL is less than Confidance score criteria user provided"));

		// Program commandLine fields
		headerInfo.add(new VCFHeaderLine(BisulfiteVCFConstants.PROGRAM_ARGS, argCommandline));

		return headerInfo;
	}

	private void initiateVCFInDifferentOutmode(SAMSequenceDictionary refDict) {
		File outputVcfFile = new File(vfn1);
		if(vfn1.endsWith(".gz")){
			writer = new TcgaVCFWriter(outputVcfFile, new BlockCompressedOutputStream(outputVcfFile), refDict,
					false, false, true,false);

		}else{
			try{
				writer = new TcgaVCFWriter(outputVcfFile, new FileOutputStream(outputVcfFile), refDict,
						false, false, true,false);
			}catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}



		writer.writeHeader(new VCFHeader(getHeaderInfo(), samples));

		if (getToolkit().getArguments().numberOfCPUThreadsPerDataThread > 1) {
			multiThreadWriter = new SortingTcgaVCFWriter(writer, MAXIMUM_CACHE_FOR_OUTPUT_VCF);
			multiThreadWriter.enableDiscreteLoci(BAC.lnc);
			if (BAC.ovd) {
				File outputVerboseFile = new File(BAC.fnovd);
				//verboseWriter = new TcgaVCFWriter(outputVerboseFile, refDict, false);
				if(BAC.fnovd.endsWith(".gz")){
					verboseWriter = new TcgaVCFWriter(outputVerboseFile, new BlockCompressedOutputStream(outputVerboseFile), refDict,
							false, false, true,false);
				}else{
					try{
						verboseWriter = new TcgaVCFWriter(outputVerboseFile, new FileOutputStream(outputVerboseFile), refDict,
								false, false, true,false);
					}catch (FileNotFoundException e) {
						e.printStackTrace();
					}
				}

				verboseWriter.writeHeader(new VCFHeader(getHeaderInfo(), samples));
			}

		}

		if (BAC.OutputMode == BisulfiteEnums.OUTPUT_MODE.DEFAULT_FOR_TCGA || BAC.OutputMode == BisulfiteEnums.OUTPUT_MODE.EMIT_VARIANT_AND_CYTOSINES ||  BAC.OutputMode == BisulfiteEnums.OUTPUT_MODE.NOMESEQ_MODE) {
			File outputAdditionalVcfFile = new File(vfn2);
			//additionalWriterForDefaultTcgaMode = new TcgaVCFWriter(outputAdditionalVcfFile, refDict, false);
			if(vfn2.endsWith(".gz")){
				try{
				additionalWriterForDefaultTcgaMode = new TcgaVCFWriter(outputAdditionalVcfFile, new BlockCompressedOutputStream(new FileOutputStream(outputAdditionalVcfFile), outputAdditionalVcfFile), refDict,
						false, false, true,false);
				}catch (FileNotFoundException e) {
					e.printStackTrace();
				}
			}else{
				try{
					additionalWriterForDefaultTcgaMode = new TcgaVCFWriter(outputAdditionalVcfFile, new FileOutputStream(outputAdditionalVcfFile), refDict,
							false, false, true,false);
				}catch (FileNotFoundException e) {
					e.printStackTrace();
				}
			}

			additionalWriterForDefaultTcgaMode.writeHeader(new VCFHeader(getHeaderInfo(), samples));
			if (getToolkit().getArguments().numberOfCPUThreadsPerDataThread > 1) {
				multiAdditionalWriterForDefaultTcgaMode = new SortingTcgaVCFWriter(additionalWriterForDefaultTcgaMode, MAXIMUM_CACHE_FOR_OUTPUT_VCF);
				multiAdditionalWriterForDefaultTcgaMode.enableDiscreteLoci(BAC.lnc);
			}

		}


		if (cpgreads != null) {
			if (BAC.sequencingMode == BisulfiteEnums.MethylSNPModel.GM && gchreads != null) {
					gchReadsWriter = new cpgReadsWriterImp(new File(gchreads));
					gchReadsWriter.addRefStrand(BAC.withRef);
					gchReadsWriter.addHeader(true);
					
					gchReadsWriter.setEncrypt(BAC.notEncrypt);
					
			} 
			readsWriter = new cpgReadsWriterImp(new File(cpgreads));

			if (getToolkit().getArguments().numberOfCPUThreadsPerDataThread > 1) {
				throw new UserException("cpgreads mode do not support multiple thread yet. please use -nt 1!!");
				// multiThreadCpgReadsWriter.enableDiscreteLoci(BAC.lnc);
			}
			readsWriter.addRefStrand(BAC.withRef);
			readsWriter.addHeader(true);
			
			readsWriter.setEncrypt(BAC.notEncrypt);
			
		}
	}

	private void outputVCFInDifferentOutmode(BisulfiteVariantCallContext value) {
		if (!value.shouldEmit)
			return;
		boolean firstStreamFlag = false;
		boolean secondStreamFlag = false;
		switch(BAC.OutputMode){
			case NOMESEQ_MODE:{
				
				if (value.isNomeseqC()){
					firstStreamFlag = true;
					//System.err.println(value.getVariantContext());
				}
					
				if (value.isVariant())
					secondStreamFlag = true;
				break;
			}
			case EMIT_ALL_CYTOSINES:{
				if (value.isC || value.isPatternInRef("C"))
					firstStreamFlag = true;
				break;
			}
			case EMIT_ALL_CPG:{
				if (value.isCpg || value.isPatternInRef("CG"))
					firstStreamFlag = true;
				break;
			}
			case EMIT_VARIANTS_ONLY:{
				if (value.isVariant())
					firstStreamFlag = true;
				break;
			}
			case EMIT_HET_SNPS_ONLY:{
				if (value.isHetSnp())
					firstStreamFlag = true;
				break;
			}
			case EMIT_VARIANT_AND_CYTOSINES:{
				if (value.isC || value.isPatternInRef("C"))
					firstStreamFlag = true;
				if (value.isVariant())
					secondStreamFlag = true;
				break;
			}
			case DEFAULT_FOR_TCGA:{
				if (value.isCpg || value.isPatternInRef("CG"))
					firstStreamFlag = true;
				if (value.isVariant())
					secondStreamFlag = true;
				break;
			}
			default:{
				firstStreamFlag = true;
			}
		}	
		
		if(firstStreamFlag){
			if (getToolkit().getArguments().numberOfCPUThreadsPerDataThread > 1) {
				multiThreadWriter.add(value.getVariantContext());
				//System.err.println(value.getVariantContext());
				if (value.ref.getLocus().getStart() % 1000000 == 0) {
					multiThreadWriter.writerFlush();
					System.gc();
				}
			} else {
				writer.add(value.getVariantContext());
			}
		}
		if(secondStreamFlag){
			if (getToolkit().getArguments().numberOfCPUThreadsPerDataThread > 1) {
				multiAdditionalWriterForDefaultTcgaMode.add(value.getVariantContext());
				if (value.ref.getLocus().getStart() % 1000000 == 0) {
					multiAdditionalWriterForDefaultTcgaMode.writerFlush();
				}
			} else {
				//System.err.println(value.getVariantContext());
				additionalWriterForDefaultTcgaMode.add(value.getVariantContext());
			}
		}
		
		
	}

	private void readsReport(BisulfiteVariantCallContext value) {
		//System.err.println(value.getSummaryAcrossRG().cytosinePatternConfirmedSet);
		if (BAC.sequencingMode == BisulfiteEnums.MethylSNPModel.GM) {
			
			if (value.getSummaryAcrossRG().cytosinePatternConfirmedSet.contains("GCH")) {
				if (gchreads != null) {
					if (value.getSummaryAcrossRG().cytosinePatternStrand == '+') {
						readsDetailReport(value, true, gchReadsWriter);
					} else {
						readsDetailReport(value, false, gchReadsWriter);
					}

				}
			}
			else if (value.getSummaryAcrossRG().cytosinePatternConfirmedSet.contains("HCG") || value.getSummaryAcrossRG().cytosinePatternConfirmedSet.contains("WCG")) {
				
				if (cpgreads != null) {
					if (value.getSummaryAcrossRG().cytosinePatternStrand == '+') {
						readsDetailReport(value, true, readsWriter);
					} else {
						readsDetailReport(value, false, readsWriter);
					}

				}
			}
		} else {			
			 if (value.isCpg && !value.isHetCpg) {
				if (cpgreads != null) {
					if (value.getSummaryAcrossRG().cytosinePatternStrand == '+') {
						readsDetailReport(value, true, readsWriter);
					} else {
						readsDetailReport(value, false, readsWriter);
					}

				}
			}
			
		}
		
		//if user want SNP inside the cpgreads file.
		if(value.isHetSnp() && value.confidentlyCalled && BAC.includeSnp){
			if (cpgreads != null) {
				if( BAC.onlyDbsnp && !value.getVariantContext().hasID()){
					
				}
				else{
					readsDetailReport(value, readsWriter, true);
					if (BAC.sequencingMode == BisulfiteEnums.MethylSNPModel.GM && gchreads != null) {
						readsDetailReport(value, gchReadsWriter, true);
					}
				}

			}
		}
		
	}


	public void readsDetailReport(BisulfiteVariantCallContext value, boolean posStrand, FormatWriterBase privateWriter) {
		AlignmentContext rawContext = value.rawContext;

		if (rawContext.hasReads()) {
			if (value.ref.getLocus().getStart() % 1000000 == 0) {
				privateWriter.writerFlush();
				System.gc();
			}
			String refStrand = ".";
			if(BaseUtils.basesAreEqual(value.ref.getBase(), BaseUtils.Base.C.base)){
				refStrand = "+";
			}else if(BaseUtils.basesAreEqual(value.ref.getBase(), BaseUtils.Base.G.base)){
				refStrand = "-";
			}else{
				if(value.getSummaryAcrossRG().cytosinePatternStrand == '+'){
					refStrand = refStrand + "+";
				}else{
					refStrand = refStrand + "-";
				}
			}
			for (PileupElement p : rawContext.getBasePileup().getOverlappingFragmentFilteredPileup()) {
				if (!BisSNPUtils.goodBaseInPileupElement(p, BAC, value.ref)) {
					continue;
				}
					boolean readNegStrand = p.getRead().getReadNegativeStrandFlag();
					char strand = readNegStrand ? '-' : '+';
				
				
					if (p.getRead().getReadPairedFlag() && p.getRead().getSecondOfPairFlag() && !BAC.nonDirectional) {
						readNegStrand = !readNegStrand;

					}
					
					
					if (readNegStrand) {
						//strand = '-';
						if (!posStrand) {
							char methyStatus;

							if (p.getBase() == BaseUtilsMore.G) {
								methyStatus = 'm';
								cpgReads cr;
								if(BAC.withRef){
										cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand,p.getRead().getReadName(), refStrand);
									
								}
								else{
									cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand,
											p.getRead().getReadName());
								} 
								//cpgReads cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand,
								//p.getRead().getReadName());	
								privateWriter.add(cr);

							} else if (p.getBase() == BaseUtilsMore.A) {
								methyStatus = 'u';
								cpgReads cr;
								if(BAC.withRef){
									cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand,p.getRead().getReadName(), refStrand);
									
								}
								else{
									cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand,
											p.getRead().getReadName());
								} 
								//cpgReads cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand,
								//p.getRead().getReadName());
								privateWriter.add(cr);

							}

						}
					} else {
						//strand = '+';
						if (posStrand) {
							char methyStatus;
							if (p.getBase() == BaseUtilsMore.C) {
								methyStatus = 'm';
								cpgReads cr;
								if(BAC.withRef){
									cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand,p.getRead().getReadName(), refStrand);
									
								}
								else{
									cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand,
											p.getRead().getReadName());
								} 
								//cpgReads cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand,
								//p.getRead().getReadName());
								privateWriter.add(cr);

							} else if (p.getBase() == BaseUtilsMore.T) {
								methyStatus = 'u';
								cpgReads cr;
								if(BAC.withRef){
									cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand,p.getRead().getReadName(), refStrand);
									
								}
								else{
									cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand,
											p.getRead().getReadName());
								} 
								//cpgReads cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), methyStatus, p.getQual(), strand,
								//p.getRead().getReadName());
								privateWriter.add(cr);
							}

						}
					}

			}

		}

	}
	
	public void readsDetailReport(BisulfiteVariantCallContext value, FormatWriterBase privateWriter, boolean hetSNP) {
		AlignmentContext rawContext = value.rawContext;
		
		if (rawContext.hasReads()) {
			if (value.ref.getLocus().getStart() % 1000000 == 0) {
				privateWriter.writerFlush();
				System.gc();
			}
			
			String refStrand = ".";
			if(BaseUtils.basesAreEqual(value.ref.getBase(), BaseUtils.Base.C.base)){
				refStrand = "+";
			}else if(BaseUtils.basesAreEqual(value.ref.getBase(), BaseUtils.Base.G.base)){
				refStrand = "-";
			}else{
				if(value.getSummaryAcrossRG().cytosinePatternStrand == '+'){
					refStrand = refStrand + "+";
				}else{
					refStrand = refStrand + "-";
				}
			}
			Allele altAllele = value.getVariantContext().getAltAlleleWithHighestAlleleCount();
			int i =0;
			for (PileupElement p : rawContext.getBasePileup().getOverlappingFragmentFilteredPileup()) {

				if (!BisSNPUtils.goodBaseInPileupElement(p, BAC, value.ref)) {
					continue;
				}
				char strand = p.getRead().getReadNegativeStrandFlag() ? '-' : '+';
				
				
				
				
				
				boolean negStrand = p.getRead().getReadNegativeStrandFlag();
				boolean secondStrand = p.getRead().getProperPairFlag() && p.getRead().getSecondOfPairFlag();
				byte base = negStrand ? BaseUtils.simpleComplement(p.getBase()) : p.getBase();
				//if(value.ref.getLocus().getStart() == 10784192 && p.getRead().getReadName().equalsIgnoreCase("HWI-ST550_0181:5:2206:16865:54452")){
				//	System.err.println((char)base + "\t" + (char)value.ref.getBase() + "\t" + (char)altAllele.getBases()[0] + "\t" + i + "\t" + privateWriter.toString());
				//}
				byte tmpBase = negStrand ? BaseUtils.simpleComplement(base) : base;
				if(secondStrand){
					tmpBase = BaseUtilsMore.iupacCodeComplement(tmpBase);
				}
				if(BaseUtils.basesAreEqual(value.ref.getBase(), tmpBase) || (altAllele.isCalled() && BaseUtils.basesAreEqual(altAllele.getBases()[0], tmpBase))){
					base = tmpBase;
				}else if(!BaseUtilsMore.isBisulfiteMismatch(value.ref.getBase(), base, negStrand, secondStrand)){
					base = value.ref.getBase();
				}else if(altAllele.isCalled() && !BaseUtilsMore.isBisulfiteMismatch(altAllele.getBases()[0], base, negStrand, secondStrand)){
					base = altAllele.getBases()[0];
				}else{
					base = '.';
				}
				
				cpgReads cr;
				if(BAC.withRef){
					cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), (char)base, p.getQual(), strand,p.getRead().getReadName(), refStrand);
					
				}
				else{
					cr = new cpgReads(rawContext.getContig(), rawContext.getLocation().getStart(), (char)base, p.getQual(), strand,
							p.getRead().getReadName());
				} 
				
				privateWriter.add(cr);
				i++;
			}

		}

	}

	
	private TreeMap<String, Pair<Integer, Double>> sorHashMapByKeyValue(HashMap<String, Pair<Integer, Double>> input) {

		TreeMap<String, Pair<Integer, Double>> output = new TreeMap<String, Pair<Integer, Double>>();
		for (String key : input.keySet()) {
			Pair<Integer, Double> methyValue = input.get(key);
			output.put(key, methyValue);
		}

		return output;
	}
	

	/**
	 * Inner class for collecting output statistics
	 */
	public static class ContextCondition {
		// Integer: number of cytosine;
		// Double: sumMethyLevel;
		HashMap<String, Pair<Integer, Double>> cytosineMethySummary = new HashMap<String, Pair<Integer, Double>>();

		// record other kind of cytosine statistics status user provided
		HashMap<String, HashMap<String, Pair<Integer, Double>>> cytosineMethySummaryByReadsGroup = new HashMap<String, HashMap<String, Pair<Integer, Double>>>();

		/**
		 * The number of bases that were potentially callable -- i.e., those not at excessive
		 * coverage or masked with N
		 */
		long nBasesCallable = 0;

		/**
		 * The number of bases called confidently (according to user threshold), either ref or other
		 */
		long nBasesCalledConfidently = 0;

		/** The total number of passes examined -- i.e., the number of map calls */
		long nBasesVisited = 0;

		/** The number of bases for which calls were emitted */
		long nCallsMade = 0;

		/** The total number of good reads bases in callable loci */
		long nReadsCallable = 0;
		
		/** The total number of good reads bases in cpg island */
		long nReadsCgi = 0;
		
		/** The number of bases in CpG island */
		long nBasesInCgi = 0;
		
		/** The number of bases callable and in CpG island */
		long nBasesCallableInCgi = 0;
		
		/** The number of bases confidently called and in CpG island */
		long nBasesConfidentlyInCgi = 0;
		

		void makeCytosineMap() {
			//for (String cytosineKey : BAC.cytosineDefined.getContextDefined().keySet()) {
			for(Entry<String, CytosineParameters> valueEntrySet : BAC.cytosineDefined.getContextDefined().entrySet()){
				Pair<Integer, Double> methySummary = new Pair<Integer, Double>(0, 0.0);
				cytosineMethySummary.put(valueEntrySet.getKey(), methySummary);

			}
			for (String sample : samples) {
				HashMap<String, Pair<Integer, Double>> cytosineMethySummaryTmp = new HashMap<String, Pair<Integer, Double>>();
				//for (String cytosineKey : BAC.cytosineDefined.getContextDefined().keySet()) {
				for (Entry<String, CytosineParameters> valueEntrySet : BAC.cytosineDefined.getContextDefined().entrySet()) {
					Pair<Integer, Double> methySummary = new Pair<Integer, Double>(0, 0.0);
					cytosineMethySummaryTmp.put(valueEntrySet.getKey(), methySummary);
				}
				cytosineMethySummaryByReadsGroup.put(sample, cytosineMethySummaryTmp);
			}

		}

		double meanCoverageOfAll() {
			return (double) nReadsCallable / (double) nBasesVisited;
		}

		double meanCoverageOfCallable() {
			return (double) nReadsCallable / (double) nBasesCallable;
		}
		
		double meanCoverageOfCpgIsland() {
			return (double) nReadsCgi / (double) nBasesInCgi;
		}
		
		double meanCoverageOfCallableInCpgIsland() {
			return (double) nReadsCgi / (double) nBasesCallableInCgi;
		}

		double percentCallableOfAll() {
			return (100.0 * nBasesCallable) / (nBasesVisited);
		}

		double percentCalledOfAll() {
			return (100.0 * nBasesCalledConfidently) / (nBasesVisited);
		}

		double percentCalledOfCallable() {
			return (100.0 * nBasesCalledConfidently) / (nBasesCallable);
		}
		
	
	}

}
