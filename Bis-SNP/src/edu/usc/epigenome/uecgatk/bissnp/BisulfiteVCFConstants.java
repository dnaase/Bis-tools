package edu.usc.epigenome.uecgatk.bissnp;

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

public class BisulfiteVCFConstants {

	public static final String AVE_BASE_QUALITY_KEY = "BQ";
	public static final String BEST_C_PATTERN = "CP";
	
	// to record number of C, T,other_reads in bisulfite_conversion strand, number of A,G,
	// other_reads in genotype strand
	public static final String C_STATUS = "BRC6";
	public static final String C_STRAND_KEY = "CS";
	public static final String CYTOSINE_METHY_VALUE = "MR"; // in 0-100% scale
	public static final String CYTOSINE_TYPE = "Context";
	public static final String GENOTYPE_LIKELIHOODS_KEY = "GL";


	// INFO/FORMAT field keys for Bis-SNP VCF
	public static final String GENOTYPE_TYPE = "HOM_REF,HET,HOM_VAR";
	public static final String REF_C_PATTERN = "REF";
	public static final String ID_KEY = "ID";
	public static final String LESS_THAN_HALF_SAMPLES_HAVE_DATE = "s50";
	public static final String NUM_OF_SAMPLES = "NS";
	public static final String NUMBER_OF_C_KEY = "CM";
	public static final String NUMBER_OF_T_KEY = "CU";
	public static final String PROGRAM_ARGS = BisSNP.getBisSNPVersionNumber() + " Program Args";
	public static final String QUAL_BY_DEPTH = "QD";
	public static final String QUALITY_BELOW_10 = "q10";
	public static final String READS_SUPPORT_ALT = "DP4";
	public static final String SOMATIC_STAT_VAR = "SS";
	public static final String VCF_HEADER_VERSION_ASSEMBLY = "assembly";
	public static final String VCF_HEADER_VERSION_CENTER = "center";
	public static final String VCF_HEADER_VERSION_DATE = "fileDate";
	public static final String VCF_HEADER_VERSION_FORMAT = "fileformat";
	public static final String VCF_HEADER_VERSION_GAF = "geneAnno";
	public static final String VCF_HEADER_VERSION_LOG = "vcfProcessLog";
	public static final String VCF_HEADER_VERSION_PHASE = "phasing";
	public static final String VCF_HEADER_VERSION_REF = "reference";
	public static final String VCF_HEADER_VERSION_CONTIG = "contig";
	public static final String VCF_HEADER_VERSION_TCGA_VERSION = "tcgaversion";
	public static final String GLOBAL_MINOR_ALLELE_FRE = "GMAF";
	public static final String VALIDATED_SNP = "VLD";
}
