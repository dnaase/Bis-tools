package edu.usc.epigenome.uecgatk.bissnp;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Map.Entry;


import java.util.List;


import org.broad.tribble.bed.BEDFeature;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import edu.usc.epigenome.uecgatk.bissnp.BisulfiteEnums.OUTPUT_MODE;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import gnu.trove.set.hash.TLinkedHashSet;

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

public class BisulfiteGenotyperEngine {

	public static final String LOW_QUAL_FILTER_NAME = "LowQual";

	// the standard filter to use for calls below the confidence threshold but
	// above the emit threshold
	private static final THashSet<String> filter = new THashSet<String>(20);

	protected double MAX_PHRED = 1000000;

	private BisulfiteArgumentCollection BAC = null;

	// the model used for calculating genotypes
	private ThreadLocal<BisulfiteSNPGenotypeLikelihoodsCalculationModel> bglcms = new ThreadLocal<BisulfiteSNPGenotypeLikelihoodsCalculationModel>();

	private BisulfiteVariantCallContext bisulfiteVariantCallContext = null;
	private int SOMATIC_STATS = 5;

	private GenomeAnalysisEngine toolkit = null;

	public BisulfiteGenotyperEngine(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, BisulfiteArgumentCollection BAC, GenomeAnalysisEngine toolkit) {
		this.BAC = BAC.clone();

		this.toolkit = toolkit;

		filter.add(LOW_QUAL_FILTER_NAME);

		calculateLikelihoodsAndGenotypes(tracker, refContext, rawContext);

	}

	protected static BisulfiteSNPGenotypeLikelihoodsCalculationModel getGenotypeLikelihoodsCalculationObject(BisulfiteArgumentCollection BAC) {
		return new BisulfiteSNPGenotypeLikelihoodsCalculationModel(BAC, BAC.useBAQ);
	}

	public VariantContext calculateLikelihoods(RefMetaDataTracker tracker, ReferenceContext refContext, THashMap<String, AlignmentContext> stratifiedContexts, AlignmentContextUtils.ReadOrientation type,
			THashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs) {
		if (stratifiedContexts == null) {
			return null;
		}

		// initialize the data for this thread if that hasn't been done yet
		if (bglcms.get() == null) {
			bglcms.set(BisulfiteGenotyperEngine.getGenotypeLikelihoodsCalculationObject(BAC));

		}

		BisulfiteSNPGenotypeLikelihoodsCalculationModel bglcm = bglcms.get();

		bglcm.setBsLikelihoods(tracker, refContext, stratifiedContexts, type, BCGLs);

		if (!BCGLs.isEmpty()) {

			return createVariantContextFromLikelihoods(refContext, bglcm.getRefAllele(), BCGLs);
		} else {
			return null;
		}

	}

	/**
	 * Compute full BisulfiteVariantCallContext at a given locus.
	 * 
	 * @param tracker
	 *            the meta data tracker
	 * @param refContext
	 *            the reference base
	 * @param rawContext
	 *            contextual information around the locus
	 * @return the BisulfiteVariantCallContext object
	 */
	public void calculateLikelihoodsAndGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
		
		//rawContext = filterBadReadsInAlignmentContext(rawContext, refContext, BAC);
		THashMap<String, AlignmentContext> stratifiedContexts = BisSNPUtils.splitContextBySampleName(rawContext);
		THashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs = new THashMap<String, BisulfiteContextsGenotypeLikelihoods>();
		// get likelihood and methylation pattern information from
		// BisulfiteSNPGenotypeLikelihoodsCalculationModel
		VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.COMPLETE, BCGLs);

		if (vc == null){
			bisulfiteVariantCallContext = new BisulfiteVariantCallContext(BCGLs, vc, rawContext, refContext);
		}
		else{
			bisulfiteVariantCallContext = calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, BCGLs, vc);
		}
		for (RODRecordList rods : tracker.getBoundRodTracks()) {

			for (GATKFeature vc_input : rods) {
				if (vc_input != null ) {
					if(vc_input.getUnderlyingObject() instanceof BEDFeature){
						bisulfiteVariantCallContext.isCgi = true;
					}
					
				}
			}
		}
		// including all Reads group genotypes information into vcc, and provide
		// most probable genotype and alt-allele for all of ReadsGroup
		

	}

	public BisulfiteVariantCallContext getBisulfiteVariantCallContext() {

		return bisulfiteVariantCallContext;
	}

	protected boolean passesCallThreshold(double conf) {
		return conf >= BAC.STANDARD_CONFIDENCE_FOR_CALLING;
	}

	protected boolean passesEmitThreshold(double conf) {
		return conf >= BAC.STANDARD_CONFIDENCE_FOR_EMITTING;
	}

	protected boolean passesEmitThreshold(double conf, int bestAFguess) {
		return (BAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES || bestAFguess != 0) && conf >= Math.min(BAC.STANDARD_CONFIDENCE_FOR_CALLING, BAC.STANDARD_CONFIDENCE_FOR_EMITTING);
	}

	private void assignAFPosteriors(double[] likelihoods, double[] log10AFPosteriors) {
		for (int i = 0; i < likelihoods.length; i++) {
			log10AFPosteriors[i] = likelihoods[i];
		}

	}

	/**
	 * Can be overridden by concrete subclasses
	 * 
	 * @param vc
	 *            variant context with genotype likelihoods
	 * @param log10AlleleFrequencyPosteriors
	 *            allele frequency results
	 * @param AFofMaxLikelihood
	 *            allele frequency of max likelihood
	 * 
	 * @return calls
	 */

	private void assignGenotypes(VariantContext vc, THashMap<String, Integer> bestAFs, GenotypesContext calls, THashMap<String, Double> logRatios, THashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs,
			THashSet<String> sampleNames) {
		if (!vc.isVariant())
			throw new UserException("The VCF record passed in does not contain an ALT allele at " + vc.getChr() + ":" + vc.getStart());

		GenotypesContext GLs = vc.getGenotypes();
		// Set<String> sampleNames = GLs.getSampleNames();
		for (String sampleName : sampleNames) {
			Genotype g = GLs.get(sampleName);

			// if ( !g.hasLikelihoods() )
			// continue;
			Double logRatio = -logRatios.get(sampleName) / 10.0;
			Integer bestAF = bestAFs.get(sampleName);
			Allele alleleA = BCGLs.get(sampleName).getAlleleA();
			Allele alleleB = BCGLs.get(sampleName).getAlleleB();
			ArrayList<Allele> myAlleles = new ArrayList<Allele>();
			//if(vc.getStart() == 7993150){
			//	System.err.println(bestAF + "\t" + alleleA + "\t" + alleleB);
			//}
			if (bestAF == 0) {
				if (alleleA.isReference()) {
					myAlleles.add(alleleA);
					myAlleles.add(alleleA);
				} else if (alleleB.isReference()) {
					myAlleles.add(alleleB);
					myAlleles.add(alleleB);
				}
				else {
					myAlleles.add(alleleA);
					myAlleles.add(alleleA);
				}

			} else if (bestAF == 1) {
				myAlleles.add(alleleA);
				myAlleles.add(alleleB);

			} else {
				if (alleleA.isReference()) {
					myAlleles.add(alleleB);
					myAlleles.add(alleleB);
				} else if (alleleB.isReference()) {
					myAlleles.add(alleleA);
					myAlleles.add(alleleA);
				}
				else{
					myAlleles.add(alleleB);
					myAlleles.add(alleleB);
				}

			}
			calls.add(new Genotype(sampleName, myAlleles, logRatio, null, g.getAttributes(), false));
		}

	}

	private int calculateEndPos(List<Allele> alleles, Allele refAllele, GenomeLoc loc) {
		boolean isSNP = true;
		for (Allele a : alleles) {
			if (a.getBaseString().length() != 1) {
				isSNP = false;
				break;
			}
		}

		int endLoc = loc.getStart();
		if (!isSNP)
			endLoc += refAllele.length();

		return endLoc;
	}
	
	private BisulfiteVariantCallContext createBVCfromRef(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext,
			THashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs) {
		
		return null;
	}

	private BisulfiteVariantCallContext calculateGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, THashMap<String, AlignmentContext> stratifiedContexts,
			THashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs, VariantContext vc) {
		// initialize the data for this thread if that hasn't been done yet
		if (bglcms.get() == null) {
			return null;
		}
		// System.err.println(vc.toString());
		THashMap<String, Object> attributes = new THashMap<String, Object>();
		int totalDepth = 0;
		int[] alleleCount = new int[2];
		alleleCount[1] = -1;
		Allele bestAlternateAlelle = Allele.NO_CALL;

		double logRatio = Double.NEGATIVE_INFINITY;
		int bestAF = 0;
		int numOfC = 0;
		int numOfT = 0;
		char Strand = '+';
		double maxSB = Double.NEGATIVE_INFINITY;
		int maxMq0 = 0;
		//double minQD = Double.POSITIVE_INFINITY;
		String cType = null;
		THashSet<String> cytosineConfirmed = new THashSet<String>(10);
		THashSet<String> cytosinePatternInRef = new THashSet<String>(10);
		THashMap<String, Double> logRatiosMap = new THashMap<String, Double>();
		THashMap<String, Integer> bestAFsMap = new THashMap<String, Integer>();
		//for (String sample : BCGLs.keySet()) {
		for ( Entry<String, BisulfiteContextsGenotypeLikelihoods> BCGLEntrySet: BCGLs.entrySet()) {
			String sample = BCGLEntrySet.getKey();
			BisulfiteContextsGenotypeLikelihoods GL = BCGLEntrySet.getValue();

			totalDepth += GL.getDepth();
			
			//for(String cp : BAC.cytosineDefined.getContextDefined().keySet()){
			for( Entry<String, CytosineParameters> cytosineDefinedEntrySet: BAC.cytosineDefined.getContextDefined().entrySet()){
				String cp = cytosineDefinedEntrySet.getKey();
				if(GL.getCytosineParameters().get(cp).isReferenceCytosinePattern){
					
					cytosinePatternInRef.add(cp);
				}
			}


			double[] log10AlleleFrequencyPosteriors = new double[3];
			assignAFPosteriors(GL.getLikelihoods(), log10AlleleFrequencyPosteriors);
			int bestAFguess = MathUtils.maxElementIndex(log10AlleleFrequencyPosteriors);
			int secondAFguess = MathUtils.minElementIndex(log10AlleleFrequencyPosteriors);
			for (int i = 0; i < log10AlleleFrequencyPosteriors.length; i++) {
				if (i != bestAFguess) {
					if (log10AlleleFrequencyPosteriors[i] >= log10AlleleFrequencyPosteriors[secondAFguess]) {
						secondAFguess = i;
					}
				}
			}

			double[] normalizedPosteriors = log10AlleleFrequencyPosteriors;
			if (bestAFguess != 0) {
				bestAF = bestAFguess;
				if (bestAlternateAlelle.isNoCall()) {
					bestAlternateAlelle = GL.getAlleleB();
					alleleCount[0] += bestAFguess == 1 ? 1 : 2;
				} else {
					if (GL.getAlleleB().basesMatch(bestAlternateAlelle)) {
						alleleCount[0] += bestAFguess == 1 ? 1 : 2;
					} else {
						alleleCount[1] += bestAFguess == 1 ? 1 : 2;
					}
				}

			}

			double logRatioTmp = 10 * (normalizedPosteriors[bestAFguess] - normalizedPosteriors[secondAFguess]);
			logRatiosMap.put(sample, logRatioTmp);
			bestAFsMap.put(sample, bestAFguess);
			if (logRatioTmp > logRatio)
				logRatio = logRatioTmp;

			if (passesCallThreshold(logRatioTmp) && GL.getBestMatchedCytosinePattern() != null) {
				numOfC += GL.getNumOfCReadsInBisulfiteCStrand();
				numOfT += GL.getNumOfTReadsInBisulfiteCStrand();
				Strand = GL.getCytosineParameters().get(GL.getBestMatchedCytosinePattern()).cytosineStrand;
				//for (String cytosinePattern : GL.getCytosineParameters().keySet()) {
				for ( Entry<String, CytosineParameters> cytosineParametersEntrySet : GL.getCytosineParameters().entrySet()) {
					String cytosinePattern = cytosineParametersEntrySet.getKey();
					
					if (cytosineParametersEntrySet.getValue().isCytosinePattern) {

						if (GL.getCytosineParameters().get(GL.getBestMatchedCytosinePattern()).isHeterozygousCytosinePattern
								|| GL.getCytosineParameters().get(GL.getBestMatchedCytosinePattern()).isHeterozygousInContextPosition) {
							cytosineConfirmed.add(BaseUtilsMore.makeIupacCodeFrom2String(GL.getCytosineParameters().get(GL.getBestMatchedCytosinePattern()).patternOfAlleleA, GL
									.getCytosineParameters().get(GL.getBestMatchedCytosinePattern()).patternOfAlleleB));
						} else {
							cytosineConfirmed.add(cytosinePattern);
						}

					}
				}
			}
			if (bestAFguess != 0) {
				if (!BAC.NO_SLOD) {

					double overallLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors, 1);

					// the forward lod
					THashMap<String, BisulfiteContextsGenotypeLikelihoods> tmpGLs = new THashMap<String, BisulfiteContextsGenotypeLikelihoods>();
					VariantContext vcForward = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.FORWARD, tmpGLs);
					for (int i = 0; i < log10AlleleFrequencyPosteriors.length; i++)
						log10AlleleFrequencyPosteriors[i] = -1.0 * Double.MAX_VALUE;
					;

					if (tmpGLs.containsKey(sample)) {
						assignAFPosteriors(tmpGLs.get(sample).getLikelihoods(), log10AlleleFrequencyPosteriors);
					}

					double forwardLog10PofNull = log10AlleleFrequencyPosteriors[0];
					double forwardLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors, 1);

					// the reverse lod
					tmpGLs = new THashMap<String, BisulfiteContextsGenotypeLikelihoods>();
					VariantContext vcReverse = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.REVERSE, tmpGLs);

					for (int i = 0; i < log10AlleleFrequencyPosteriors.length; i++)
						log10AlleleFrequencyPosteriors[i] = -1.0 * Double.MAX_VALUE;
					;

					if (tmpGLs.containsKey(sample)) {
						assignAFPosteriors(tmpGLs.get(sample).getLikelihoods(), log10AlleleFrequencyPosteriors);
					}

					double reverseLog10PofNull = log10AlleleFrequencyPosteriors[0];
					double reverseLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors, 1);

					double forwardLod = forwardLog10PofF + reverseLog10PofNull - overallLog10PofF;
					double reverseLod = reverseLog10PofF + forwardLog10PofNull - overallLog10PofF;

					double strandScore = Math.max(forwardLod, reverseLod);

					strandScore *= 10.0;
					if(strandScore > maxSB){
						maxSB = strandScore;
						attributes.put(VCFConstants.STRAND_BIAS_KEY, String.format("%.4f", strandScore));
					}
					

				}

				double QD = logRatio / GL.getDepth();
				//if(QD < minQD){
				//	minQD = QD;
					//attributes.put(VCFConstants, minQD);
				//}

			}
			if(GL.getMQ0() >= maxMq0){
				maxMq0 = GL.getMQ0();
				attributes.put(VCFConstants.MAPPING_QUALITY_ZERO_KEY, maxMq0);
			}
			

		}

		GenotypesContext genotypes = GenotypesContext.create();

		VariantContext dbsnp = null;
		for (RODRecordList rods : tracker.getBoundRodTracks()) {

			for (GATKFeature vc_input : rods) {
				if (vc_input != null ) {
					if(vc_input.getUnderlyingObject() instanceof VariantContext){
						dbsnp = (VariantContext) vc_input.getUnderlyingObject();
						break;
					}
					
				}
			}
		}
		String rsID = null;
		if (dbsnp != null && dbsnp.hasID()) {
			rsID = dbsnp.getID();

			attributes.put(VCFConstants.DBSNP_KEY, true);
		}

		// if the site was downsampled, record that fact
		if (rawContext.hasPileupBeenDownsampled())
			attributes.put(VCFConstants.DOWNSAMPLED_KEY, true);

		List<Allele> myAlleles = vc.getAlleles();
		List<Allele> newAlleles = new ArrayList<Allele>();
		if (bestAF == 0 && myAlleles.size() <= 2) { // to get rid of second best alleles when
													// genotype is 0/0
			for (Allele a : myAlleles) {
				if (a.isReference()) {
					newAlleles.add(a);
				}

			}
		} else {
			newAlleles = myAlleles;
		}

		GenomeLoc loc = refContext.getLocus();

		int endLoc = calculateEndPos(vc.getAlleles(), vc.getReference(), loc);
		THashSet<String> keySets = new THashSet<String>();
		keySets.addAll(BCGLs.keySet());
		assignGenotypes(vc, bestAFsMap, genotypes, logRatiosMap, BCGLs, keySets);

		attributes.put(BisulfiteVCFConstants.NUM_OF_SAMPLES, BCGLs.keySet().size());
		attributes.put(VCFConstants.DEPTH_KEY, totalDepth);
		
		//if context is in ref pattern, just show REF=0, if context not in ref pattern (or no context but it is in reference genome), then show as REF=CG,C;
		if(!cytosinePatternInRef.isEmpty()){
			//System.err.println(cytosinePatternInRef.toString());
			//if((BAC.OutputMode == BisulfiteEnums.OUTPUT_MODE.DEFAULT_FOR_TCGA ||BAC.OutputMode == BisulfiteEnums.OUTPUT_MODE.EMIT_ALL_CPG) 
			//		&& cytosinePatternInRef.contains("CG") ){ //temporarily, just to remove those CCG postion, which are CG in + strand, but CH in H position 

					//cytosinePatternInRef.remove("CH");
					//cytosinePatternInRef.remove("CHH");
					//cytosinePatternInRef.remove("CHG");
			//	}
				
			
			if(cytosineConfirmed.containsAll(cytosinePatternInRef)){
				attributes.put(BisulfiteVCFConstants.REF_C_PATTERN, 0);
			}
			else{
					String tmpC = null;
					for (String c : cytosinePatternInRef) {
						if (tmpC == null) {
							tmpC = c;
						} else {
							tmpC = tmpC + "," + c;
						}
					}
					
					attributes.put(BisulfiteVCFConstants.REF_C_PATTERN, tmpC);
			}
		}
		
		
		
		
		
		// add INFO colum about methylation summary across Read Group
		if (passesCallThreshold(logRatio) & !cytosineConfirmed.isEmpty()) {
			for (String c : cytosineConfirmed) {
				if (cType == null) {
					cType = c;
				} else {
					cType = cType + "," + c;
				}
			}

			attributes.put(BisulfiteVCFConstants.C_STRAND_KEY, Strand);
			attributes.put(BisulfiteVCFConstants.CYTOSINE_TYPE, cType);
		}
		else if(!cytosinePatternInRef.isEmpty()){ //if reference genome is some cytosine pattern, but not confidently called, then it will estimate strand by reference genome
			if(BaseUtils.basesAreEqual(refContext.getBase(), BaseUtils.G)){
				attributes.put(BisulfiteVCFConstants.C_STRAND_KEY, '-');
			}
			else if(BaseUtils.basesAreEqual(refContext.getBase(), BaseUtils.C)){
				attributes.put(BisulfiteVCFConstants.C_STRAND_KEY, '+');
			}
		}

		VariantContextBuilder vcb = new VariantContextBuilder("BG_call", loc.getContig(), loc.getStart(), endLoc, newAlleles);
		vcb.attributes(attributes);
		vcb.filters(passesCallThreshold(logRatio) ? null : filter);
		vcb.genotypes(genotypes);
		if (rsID != null) {
			vcb.id(rsID);
		}

		vcb.log10PError(-logRatio / 10.0);
		vcb.referenceBaseForIndel(refContext.getBase());
		VariantContext vcCall = vcb.make();
		BisulfiteVariantCallContext call = new BisulfiteVariantCallContext(BCGLs, vcCall, rawContext, refContext);
		call.confidentlyCalled = passesCallThreshold(logRatio);
		call.shouldEmit = passesEmitThreshold(logRatio);
		return call;
	}

	private VariantContext createVariantContextFromLikelihoods(ReferenceContext refContext, Allele refAllele, THashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs) {

		List<Allele> noCall = new ArrayList<Allele>();
		noCall.add(Allele.NO_CALL);

		TLinkedHashSet<Allele> allelesSet = new TLinkedHashSet<Allele>(10);
		allelesSet.add(refAllele);

		GenotypesContext genotypes = GenotypesContext.create();

		//for (BisulfiteContextsGenotypeLikelihoods BCGL : BCGLs.values()) {
		for (Entry<String, BisulfiteContextsGenotypeLikelihoods> BCGLEntrySet: BCGLs.entrySet()) {
			BisulfiteContextsGenotypeLikelihoods BCGL = BCGLEntrySet.getValue();
			allelesSet.add(BCGL.getAlleleA());
			allelesSet.add(BCGL.getAlleleB());
			
			/*
			if (!alleles.contains(BCGL.getAlleleA()))
				alleles.add(BCGL.getAlleleA());
			if (!alleles.contains(BCGL.getAlleleB()))
				alleles.add(BCGL.getAlleleB());
			*/
			THashMap<String, Object> attributes = new THashMap<String, Object>();
			GenotypeLikelihoods likelihoods = GenotypeLikelihoods.fromLog10Likelihoods(BCGL.getLikelihoods());
			attributes.put(VCFConstants.DEPTH_KEY, BCGL.getDepth());
			attributes.put(VCFConstants.GENOTYPE_POSTERIORS_KEY, likelihoods.getAsString());

			if (BCGL.getBestMatchedCytosinePattern() != null) {
				//if(BAC.nonDirectional){
				//	attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, BCGL.getNumOfCReadsInBisulfiteCStrand() + BCGL.getNumOfGReadsInGenotypeGStrand());
				//	attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, BCGL.getNumOfTReadsInBisulfiteCStrand() + BCGL.getNumOfAReadsInGenotypeGStrand());
				//	attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, BCGL.getNumOfCReadsInBisulfiteCStrand());
				//	attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, BCGL.getNumOfTReadsInBisulfiteCStrand());
				//}
				//else{
					attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, BCGL.getNumOfCReadsInBisulfiteCStrand());
					attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, BCGL.getNumOfTReadsInBisulfiteCStrand());
				//}
				
				if (BCGL.getCytosineParameters().get(BCGL.getBestMatchedCytosinePattern()).isHeterozygousCytosinePattern
						|| BCGL.getCytosineParameters().get(BCGL.getBestMatchedCytosinePattern()).isHeterozygousInContextPosition) {
					attributes.put(
							BisulfiteVCFConstants.BEST_C_PATTERN,
							BaseUtilsMore.makeIupacCodeFrom2String(BCGL.getCytosineParameters().get(BCGL.getBestMatchedCytosinePattern()).patternOfAlleleA,
									BCGL.getCytosineParameters().get(BCGL.getBestMatchedCytosinePattern()).patternOfAlleleB));
				} else {
					attributes.put(BisulfiteVCFConstants.BEST_C_PATTERN, BCGL.getBestMatchedCytosinePattern());
				}
				attributes.put(BisulfiteVCFConstants.C_STATUS, BCGL.getBaseCountStatusAsString());
			}
			else {
				attributes.put(BisulfiteVCFConstants.BEST_C_PATTERN, VCFConstants.MISSING_VALUE_v4);

				attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, VCFConstants.MISSING_VALUE_v4);
				attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, VCFConstants.MISSING_VALUE_v4);
				attributes.put(BisulfiteVCFConstants.C_STATUS, BCGL.getBaseCountStatusAsString());
			}
			attributes.put(BisulfiteVCFConstants.SOMATIC_STAT_VAR, SOMATIC_STATS);
			attributes.put(BisulfiteVCFConstants.READS_SUPPORT_ALT, BCGL.getDP4AsString());
			attributes.put(BisulfiteVCFConstants.AVE_BASE_QUALITY_KEY, BCGL.getAveBaseQualAsString());
			// Reads supporting ALT. Number of 1) forward ref alleles; 2)
			// reverse ref; 3) forward non-ref; 4) reverse non-ref alleles"
			// Average base quality for reads supporting alleles. For each
			// allele, in the same order as listed
			List<Allele> thisSampleAllele = new ArrayList<Allele>(10);
			thisSampleAllele.addAll(allelesSet);
			genotypes.add(new Genotype(BCGL.getSample(), thisSampleAllele, Genotype.NO_LOG10_PERROR, null, attributes, false));

		}

		List<Allele> alleleList = new ArrayList<Allele>(10);
		alleleList.addAll(allelesSet);
		
		GenomeLoc loc = refContext.getLocus();
		int endLoc = calculateEndPos(alleleList, refAllele, loc);

		VariantContextBuilder vcb = new VariantContextBuilder("BG_call", loc.getContig(), loc.getStart(), endLoc, alleleList);
		vcb.referenceBaseForIndel(refContext.getBase());
		vcb.log10PError(VariantContext.NO_LOG10_PERROR);
		vcb.genotypes(genotypes);

		return vcb.make();
	}
	/*
	public AlignmentContext filterBadReadsInAlignmentContext(AlignmentContext rawContext, ReferenceContext ref, BisulfiteArgumentCollection BAC){
		GenomeLoc loc = rawContext.getLocation();
		ReadBackedPileup oldPileup = rawContext.getBasePileup().getOverlappingFragmentFilteredPileup();
		List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
		List<Integer> elementOffsets = new ArrayList<Integer>();
		for (PileupElement p : oldPileup) {
			
				if(!goodPileupElement(p, BAC, ref))
					continue;


			int elementOffset = p.getOffset();
			if (elementOffset < 0 || elementOffset > p.getRead().getReadLength() - 1)
				continue;
			elementOffsets.add(elementOffset);
			reads.add(p.getRead());
		}
		ReadBackedPileup newPileup = new ReadBackedPileupImpl(loc, reads, elementOffsets);
		return new AlignmentContext( loc, newPileup);
	}
	
	public synchronized boolean goodPileupElement(PileupElement p, BisulfiteArgumentCollection BAC, ReferenceContext refContext) {
		short cPos=1;
		return goodPileupElement(p, BAC, refContext, "C", cPos);
	}
	
	public synchronized boolean goodPileupElement(PileupElement p, BisulfiteArgumentCollection BAC, ReferenceContext refContext, String patConv5, short posCinPatConv5) {
		BitSet bitSetRead;
		int offSet = p.getOffset();
		GATKSAMRecord read = p.getRead();

		if (read.containsTemporaryAttribute("Xn")) {
			//synchronized(read){
				bitSetRead = (BitSet)read.getTemporaryAttribute("Xn");
			//}
			//System.err.println(bitSetRead + " BIT_WISE_TAG: " + p.getRead().getTemporaryAttribute(BIT_WISE_TAG));
		}else{
			
			bitSetRead = new BitSet(read.getReadLength());
			 if(BisSNPUtils.badReads(read,refContext, BAC) || BisSNPUtils.badReadsByMismatchesInWindow(read,refContext, BAC, p) || BisSNPUtils.badReadsByBisulfiteConversion(read,refContext, BAC,p)){
				 //System.err.println("no good reads " + read.getReadName());
			 }else{
				 bitSetRead = BisSNPUtils.notTrimmedBase(offSet, read, BAC, bitSetRead);
				 BitSet bitSetRead2 = BisSNPUtils.goodBasePassBisulfiteConversionFivePrime(offSet, read, BAC,refContext, patConv5, posCinPatConv5);
				 bitSetRead.and(bitSetRead2);
				// System.err.println("BitWiseTag" + bitSetRead.toString() + "\t" + bitSetRead.cardinality());
			 }
			 
			read.setTemporaryAttribute("Xn", bitSetRead);

			 
		}

		if(bitSetRead.get(offSet)){
			//System.err.println(goodBase(p, BAC, refContext, offSet));
			return BisSNPUtils.goodBase(p, BAC, refContext, offSet);

		}else{
			return false;
		}
		
		
	}
*/
}
