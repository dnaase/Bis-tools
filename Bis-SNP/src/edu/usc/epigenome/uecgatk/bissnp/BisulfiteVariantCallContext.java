/**
 * 
 */
package edu.usc.epigenome.uecgatk.bissnp;

import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;


import java.util.Iterator;
import java.util.Map.Entry;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;



/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 20, 2012 11:46:27 AM
 * 
 */
public class BisulfiteVariantCallContext {

	public boolean confidentlyCalled = false;
	
	public boolean isCgi = false;

	// is one of provided cytosine pattern
	public boolean isC = false;

	public boolean isCpg = false;
	
	public boolean isNomeseqC = false; // GCH, HCG/WCG or GCG category

	public boolean isCph = false;
	public boolean isCThetLoci = false;
	public boolean isHetCpg = false;
	public boolean isHetCph = false;
	public AlignmentContext rawContext = null;
	public ReferenceContext ref = null;
	// Should this site be emitted?
	public boolean shouldEmit = true;

	private THashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs;
	private summaryAcrossReadsGroup summary;
	private VariantContext vc;

	/**
	 * 
	 */
	public BisulfiteVariantCallContext(THashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs, VariantContext vc, AlignmentContext rawContext, ReferenceContext ref) {
		// TODO Auto-generated constructor stub
		this.rawContext = rawContext;
		this.ref = ref;
		this.BCGLs = BCGLs;
		this.vc = vc;
		if(BCGLs == null && vc != null){
			 setSummaryOnVariantContext();
		}
		else if(BCGLs != null){
			setSummaryAcrossReadsGroup();
		}
		

	}

	public THashMap<String, BisulfiteContextsGenotypeLikelihoods> getBisulfiteContextsGenotypeLikelihoods() {
		return BCGLs;
	}

	public summaryAcrossReadsGroup getSummaryAcrossRG() {
		return summary;
	}

	public VariantContext getVariantContext() {
		return vc;
	}
	
	public boolean isNomeseqC() {
		if (summary.cytosinePatternConfirmedSet.contains("GCH") || summary.cytosinePatternConfirmedSet.contains("HCG") || summary.cytosinePatternConfirmedSet.contains("GCG") || summary.cytosinePatternConfirmedSet.contains("WCG") 
				|| summary.cytosinePatternInRefSet.contains("GCH") || summary.cytosinePatternInRefSet.contains("HCG") || summary.cytosinePatternInRefSet.contains("GCG") || summary.cytosinePatternInRefSet.contains("WCG")) {
			return true;
		}
		return false;
	}
	
	public boolean isPatternInRef(String pattern) {
		if (summary.cytosinePatternInRefSet.contains(pattern)) {
			return true;
		}
		return false;
	}

	public boolean isHetCytosine() {
		if (this.vc.hasGenotypes()) {
			Iterator<Genotype> it = vc.getGenotypes().iterator();
			while (it.hasNext()) {
				Genotype tmp = it.next();
				if (tmp.isHet() || tmp.isHomVar())
					return true;
			}
		}
		return false;
	}

	public boolean isHetSnp() {
		if (this.vc.hasGenotypes()) {
			Iterator<Genotype> it = vc.getGenotypes().iterator();
			while (it.hasNext()) {
				if (it.next().isHet())
					return true;
			}
		}
		return false;
	}

	public boolean isVariant() {
		if (this.vc.hasGenotypes()) {
			Iterator<Genotype> it = vc.getGenotypes().iterator();
			while (it.hasNext()) {
				Genotype tmp = it.next();
				if (tmp.isHet() || tmp.isHomVar())
					return true;
			}
		}
		return false;
	}
	
	private void setSummaryOnVariantContext() {
		summary = new summaryAcrossReadsGroup();
		summary.cytosinePatternConfirmedSet = new THashSet<String>(10);
		summary.cytosinePatternInRefSet = new THashSet<String>(10);
		if( !vc.getGenotypes().isEmpty()){
			summary.cytosinePatternConfirmedSet.add(vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, "."));
			summary.cytosinePatternInRefSet.add(vc.getAttributeAsString(BisulfiteVCFConstants.REF_C_PATTERN, "."));
			summary.numC += vc.getGenotype(0).getAttributeAsInt(BisulfiteVCFConstants.NUMBER_OF_C_KEY, 0);
			summary.numT += vc.getGenotype(0).getAttributeAsInt(BisulfiteVCFConstants.NUMBER_OF_T_KEY, 0);
			summary.numA += 0;
			summary.numG += 0;
			summary.cytosinePatternConfirmedList = summary.cytosinePatternConfirmedSet.toString();
		}
		

	}
	
 
	private void setSummaryAcrossReadsGroup() {
		summary = new summaryAcrossReadsGroup();
		summary.cytosinePatternConfirmedSet = new THashSet<String>();
		summary.cytosinePatternInRefSet = new THashSet<String>();
		if (!BCGLs.isEmpty()) {
			//for (BisulfiteContextsGenotypeLikelihoods BCGL : BCGLs.values()) {
			for ( Entry<String, BisulfiteContextsGenotypeLikelihoods> BCGLEntrySet: BCGLs.entrySet()) {
				BisulfiteContextsGenotypeLikelihoods BCGL = BCGLEntrySet.getValue();
				summary.numC += BCGL.getNumOfCReadsInBisulfiteCStrand();
				summary.numT += BCGL.getNumOfTReadsInBisulfiteCStrand();
				summary.numA += BCGL.getNumOfAReadsInGenotypeGStrand();
				summary.numG += BCGL.getNumOfGReadsInGenotypeGStrand();
				//for (String cytosinePattern : BCGL.getCytosineParameters().keySet()) {
				for ( Entry<String, CytosineParameters> cytosinePatternEntrySet: BCGL.getCytosineParameters().entrySet()) {
					String cytosinePattern = cytosinePatternEntrySet.getKey();
					CytosineParameters cp = cytosinePatternEntrySet.getValue();
					if(cp.isReferenceCytosinePattern){
						summary.cytosinePatternInRefSet.add(cytosinePattern);
					}
					
					if (cp.isCytosinePattern) {
						summary.cytosinePatternStrand = BCGL.getCytosineParameters().get(cytosinePattern).cytosineStrand;

						summary.cytosinePatternConfirmedSet.add(cytosinePattern);
						isC = true;
						if (cp.isCTHeterozygousLoci)
							isCThetLoci = true;

						if (cytosinePattern.equalsIgnoreCase("CG")) {
							isCpg = true; // exclude C/T heterozygous CpG
							if (BCGL.getCytosineParameters().get("CG").isHeterozygousCytosinePattern)
								isHetCpg = true;
						}
						if (cytosinePattern.equalsIgnoreCase("CH")) {
							isCph = true; // exclude C/T heterozygous CpG
							if (BCGL.getCytosineParameters().get("CH").isHeterozygousCytosinePattern)
								isHetCph = true;
						}

					}

				}

			}
			summary.cytosinePatternConfirmedList = summary.cytosinePatternConfirmedSet.toString();
		}

	}
	

	public class summaryAcrossReadsGroup {
		public String cytosinePatternConfirmedList;
		public THashSet<String> cytosinePatternConfirmedSet;
		public THashSet<String> cytosinePatternInRefSet;
		public char cytosinePatternStrand;
		public int numA = 0;
		public int numC = 0;
		public int numG = 0;
		public int numT = 0;

		public int getGoodReadsCoverage() {
			return numA + numC + numG + numT;
		}
	}

}
