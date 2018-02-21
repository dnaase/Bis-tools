package main.java.edu.usc.epigenome.uecgatk.bissnp;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.gatk.utils.commandline.Advanced;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
//import org.broadinstitute.gatk.tools.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.gatk.engine.arguments.GenotypeCalculationArgumentCollection;
import htsjdk.variant.variantcontext.VariantContext;


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

public class BisulfiteArgumentCollection extends GenotypeCalculationArgumentCollection {

	@Argument(fullName = "non_directional_protocol", shortName = "nonDirectional", doc = "false: Illumina protocol which is often used, only bisulfite conversion strand is kept (Lister protocol, sequence 2 forward strands only); true: Cokus protocol, sequence all 4 bisulfite converted strands", required = false)
	public boolean nonDirectional = false;


	@Argument(fullName = "bisulfite_conversion_rate", shortName = "bsRate", doc = "bisulfite conversion rate", required = false)
	public double bsRate = 0.9975;


	@Advanced
	@Argument(fullName = "cytosine_contexts_acquired", shortName = "C", doc = "Specify the cytosine contexts to check (e.g. -C CG,1,0.7  CG is the methylation pattern to check, 1 is the C's position in CG pattern, 0.7 is the initial cytosine pattern methylation level. You could specify '-C' multiple times for different cytosine pattern)", required = false)
	public List<String> cytosineContextsAcquired = new ArrayList<String>();

	public CytosinePatternsUserDefined cytosineDefined = null;

	@Input(fullName = "dbsnp", shortName = "D", doc = "dbSNP file", required = false)
	public RodBinding<VariantContext> dbsnp;

	@Argument(fullName = "file_name_output_verbose_detail", shortName = "fnovd", doc = "output file that contain verbose information, for test only", required = false)
	public String fnovd = null;

	@Argument(fullName = "locus_not_continuous", shortName = "lnc", doc = "locu to look at is not continuous, if the distance is too large, it will make some trouble in multithread VCF writer, just enable this option in performance test only", required = false)
	public boolean lnc = false;


	
	@Argument(fullName = "original_quality_score", shortName = "obq", doc = "Use original rather than recalibrated base quality score for genotyping and methylaiton calling", required = false)
	public boolean useOriginalQual = false;

	@Argument(fullName = "NO_SLOD", shortName = "NO_SLOD", doc = "not estimate the strand bias", required = false)
	public boolean NO_SLOD = false;

	@Argument(fullName = "novelDbsnpHet", shortName = "ndh", doc = "heterozygous SNP rate when the loci is discovered as SNP in dbSNP and but not validated, the default value is human genome", required = false)
	public double novelDbsnpHet = 0.02;

	@Argument(fullName = "output_reads_after_downsampling", shortName = "orad", doc = "output Bam file that after downsapling, for performance test only", required = false)
	public boolean orad = false;

	@Argument(fullName = "output_reads_coverage_after_downsampling", shortName = "orcad", doc = "define output Bam file's mean coverage that after downsapling, for performance test only", required = false)
	public int orcad = 1;

	@Argument(fullName = "output_modes", shortName = "out_modes", doc = "Output modes[EMIT_VARIANTS_ONLY,EMIT_ALL_CONFIDENT_SITES,EMIT_ALL_SITES,EMIT_ALL_CPG, EMIT_ALL_CYTOSINES,EMIT_HET_SNPS_ONLY, DEFAULT_FOR_TCGA, EMIT_VARIANT_AND_CYTOSINES, NOMESEQ_MODE]", required = false)
	public BisulfiteEnums.OUTPUT_MODE OutputMode = BisulfiteEnums.OUTPUT_MODE.DEFAULT_FOR_TCGA;

	@Argument(fullName = "output_verbose_detail", shortName = "ovd", doc = "output_verbose_detail, for performance test only", required = false)
	public boolean ovd = false;

	@Argument(fullName = "over_conversion_rate", shortName = "overRate", doc = "cytosine over conversion rate. it is often 0", required = false)
	public double overRate = 0;

	@Argument(fullName = "reference_genome_error", shortName = "rge", doc = "Reference genome error, the default value is human genome, in hg16 it is 99.99% accurate,  in hg17/hg18/hg19, it is less than 1e-4 (USCS genome browser described); We define it here default for human genome assembly(hg18,h19) to be 1e-6 as GATK did ", required = false)
	public double referenceGenomeErr = 1e-6;

	@Argument(fullName = "sequencing_mode", shortName = "sm", doc = "Bisulfite mode: BM, GNOMe-seq mode: GM", required = false)
	public BisulfiteEnums.MethylSNPModel sequencingMode = BisulfiteEnums.MethylSNPModel.BM;

	@Argument(fullName = "test_location", shortName = "loc", doc = "for debug only, output the detail information in the location", required = false)
	public long testLocus = -1;

	@Argument(fullName = "ti_vs_tv", shortName = "tvt", doc = "Transition rate vs. Transversion rate, in human genome, the default is 2", required = false)
	public int tiVsTv = 2;

	// @Output(fullName = "file_name_output_cpg_reads_detail", shortName =
	// "fnocrd", doc =
	// "output CpG reads bed file that contain each CpG's position in reads information, for test only",
	// required = false)
	// public SortingFormatWriterBase readsWriter = null;

	@Argument(fullName = "maximum_read_cov", shortName = "toCoverage", doc = "maximum read coverage allowed. Default is: 250", required = false)
	public long toCoverage = 250;

	@Argument(fullName = "trim_3_end_bp", shortName = "trim3", doc = "how many bases at 3'end of the reads are discarded", required = false)
	public int trim3 = 0;

	@Argument(fullName = "trim_5_end_bp", shortName = "trim5", doc = "how many bases at 5'end of the reads are discarded", required = false)
	public int trim5 = 0;
	
	@Argument(fullName = "trim_only_2nd_end_reads", shortName = "trim2nd", doc = "Only trimmed at 2nd end reads, should be enabled for tagmentation-based whole genome bisulfite sequencing", required = false)
	public boolean trim2nd = false;


	
	@Argument(fullName = "use_baq_for_calculation", shortName = "useBAQ", doc = "use BAQ for genotype calculation", required = false)
	public boolean useBAQ = false;

	@Argument(fullName = "validateDbsnphet", shortName = "vdh", doc = "heterozygous SNP rate when the loci is discovered as SNP in dbSNP and is validated, the default value is human genome", required = false)
	public double validateDbsnpHet = 0.1;

	@Argument(fullName = "maximum_cache_for_output_vcf", shortName = "vcfCache", doc = "maximum cached position for multithreads output of VCF. Default is: 10,000,000", required = false)
	public int vcfCache = 1000000;
	
	@Argument(fullName = "include_cytosine_no_coverage", shortName = "includeAllInRef", doc = "including all of cytosine pattern in reference genome, even it does not have any reads covered. Default: false", required = false)
	public boolean includeAllInRef = false;
	
	@Argument(fullName = "include_snp_in_cpg_reads_file", shortName = "includeSnp", doc = "only output heterozygous SNP in cpg reads file. Default: false", required = false)
	public boolean includeSnp = false;
	
	@Argument(fullName = "only_output_dbsnp_in_cpg_reads_file", shortName = "onlyDbsnp", doc = "only output heterozygous SNP at known dbSNP position in cpg reads file (need to enable -includeSnp at the sam time). Default: false", required = false)
	public boolean onlyDbsnp = false;
	
	@Argument(fullName = "not_encrypt", shortName = "notEncrypt", doc = "In cpg reads file, output original reads ID rather than the output of encrypt id. Default: false", required = false)
	public boolean notEncrypt = false;
	
	@Argument(fullName = "with_ref", shortName = "withRef", doc = "In cpg reads file, output cpg's strand in the reference genome. Default: false", required = false)
	public boolean withRef = false;

    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", required = false)
    public int MIN_BASE_QUALTY_SCORE = 5;

	public Double heterozygosity = 0.001D;

	@Override
	public BisulfiteArgumentCollection clone() {
		BisulfiteArgumentCollection bac = new BisulfiteArgumentCollection();

		//bac.GLmodel = GLmodel;
		//bac.PCR_error = PCR_error;
		//bac.GenotypingMode = GenotypingMode;
		bac.OutputMode = OutputMode;
		bac.NO_SLOD = NO_SLOD;
		bac.useBAQ = useBAQ;

		bac.STANDARD_CONFIDENCE_FOR_CALLING = STANDARD_CONFIDENCE_FOR_CALLING;
		bac.STANDARD_CONFIDENCE_FOR_EMITTING = STANDARD_CONFIDENCE_FOR_EMITTING;
		bac.MIN_BASE_QUALTY_SCORE = MIN_BASE_QUALTY_SCORE;
		//bac.MIN_MAPPING_QUALTY_SCORE = MIN_MAPPING_QUALTY_SCORE;
		//bac.MAX_MISMATCHES = MAX_MISMATCHES;
		//bac.USE_BADLY_MATED_READS = USE_BADLY_MATED_READS;
		//bac.MAX_DELETION_FRACTION = MAX_DELETION_FRACTION;
		//bac.MIN_INDEL_COUNT_FOR_GENOTYPING = MIN_INDEL_COUNT_FOR_GENOTYPING;
		//bac.INDEL_HETEROZYGOSITY = indelHeterozygosity;

		bac.dbsnp = dbsnp;
		bac.sequencingMode = sequencingMode;

		bac.nonDirectional = nonDirectional;
		bac.cytosineContextsAcquired = cytosineContextsAcquired;
		bac.cytosineDefined = cytosineDefined;

		bac.testLocus = testLocus;
		bac.useOriginalQual = useOriginalQual;
		//bac.minConv = minConv;
		//bac.patConv5 = patConv5;
		//bac.patConv = patConv;
		//bac.minPatConv = minPatConv;
		bac.trim5 = trim5;
		bac.trim3 = trim3;
		bac.trim2nd = trim2nd;
		bac.bsRate = bsRate;
		bac.overRate = overRate;
		bac.validateDbsnpHet = validateDbsnpHet;
		bac.novelDbsnpHet = novelDbsnpHet;
		bac.referenceGenomeErr = referenceGenomeErr;
		bac.snpHeterozygosity = snpHeterozygosity;
		bac.heterozygosity = snpHeterozygosity;
		bac.tiVsTv = tiVsTv;
		bac.toCoverage = toCoverage;

		bac.orcad = orcad;

		bac.fnovd = fnovd;
		bac.ovd = ovd;
		bac.lnc = lnc;
		bac.vcfCache = vcfCache;
		bac.onlyDbsnp = onlyDbsnp;
		bac.includeSnp = includeSnp;
		//bac.invDups = invDups;
		bac.includeAllInRef = includeAllInRef;
		bac.notEncrypt = notEncrypt;
		bac.withRef = withRef;

		return bac;
	}

	public void makeCytosine() {
		cytosineDefined = new CytosinePatternsUserDefined(cytosineContextsAcquired, sequencingMode);
		
	}

	@Override
	public String toString() {
		String cls = new String();

		for (Field f : this.getClass().getDeclaredFields()) {
			cls = cls.concat("\n");
			cls = cls.concat(f.toString());

		}
		return cls;

	}

}
