/**
 * 
 */
package edu.usc.epigenome.uecgatk.bissnp;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 19, 2012 4:09:03 PM
 * 
 */
public class BisulfiteEnums {

	public enum cytosineContextParameters {
		cytosineMethylation, cytosinePosition, cytosineStrand, methylationAutoestimated
	}

	public enum MethylSNPModel {
		BM, GM, NM
	}

	public enum OUTPUT_MODE {
		// output one file contained Cpg info
		// the other file contained variant info
		DEFAULT_FOR_TCGA, EMIT_ALL_CONFIDENT_SITES,
		//output GCH, HCG, GCG info and variants info
		NOMESEQ_MODE,
		// only confident Cpg
		EMIT_ALL_CPG,
		// only confident cytosines
		EMIT_ALL_CYTOSINES, EMIT_ALL_SITES, EMIT_HET_SNPS_ONLY,
		// output one file contained Cytosine info
		// the other file contained variant info
		EMIT_VARIANT_AND_CYTOSINES,
		// only confident variants
		EMIT_VARIANTS_ONLY
	}

	public enum INVERT_DUPS {
		USE_ONLY_1ST_END,
		USE_BOTH_END,
		NOT_TO_USE
	}
	
}
