/**
 * 
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp;

import java.util.HashMap;
import java.util.HashSet;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.variantcontext.Allele;

import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;




/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 19, 2012 3:38:42 PM
 * 
 */
public class BisulfiteContextsGenotypeLikelihoods {

	private Allele A, B;
	private BisulfiteArgumentCollection BAC;
	private String bestMatchedCytosinePattern;
	// private int cytosinePos;
	private HashSet<String> cytosineContexts;
	private HashMap<String, CytosineParameters> cytosineParameters;
	private int goodReadsCov = 0;
	private HashMap<Integer, double[]> GPsAfterCytosineTenGenotypes;
	private double[] GPsAtCytosineNormalizedByThreeGenotypes;
	private HashMap<String, Double[]> GPsAtCytosineTenGenotypes;
	private HashMap<Integer, double[]> GPsBeforeCytosineTenGenotypes;
	private int mapQual0 = 0;

	private int numOfAReadsInGenotypeGStrand = 0;
	private int numOfCReadsInBisulfiteCStrand = 0;
	private int numOfGReadsInGenotypeGStrand = 0;
	private int numOfOtherReadsInBisulfiteCStrand = 0;
	private int numOfOtherReadsInGenotypeGStrand = 0;
	private int numOfTReadsInBisulfiteCStrand = 0;
	private ReadBackedPileup pileup;
	private int readsAlleleA = 0;
	private int readsAlleleB = 0;
	private int readsFwdAlt = 0;
	private int readsFwdRef = 0;

	private int readsRevAlt = 0;
	private int readsRevRef = 0;
	private ReferenceContext ref = null;
	private double rmsBaseQual = Double.NaN;
	private double rmsMapQual = Double.NaN;
	private String sample;
	private double totalBaseQualAlleleA = Double.NaN;

	private double totalBaseQualAlleleB = Double.NaN;
	private int totalDepth = 0;

	/**
	 * 
	 */
	public BisulfiteContextsGenotypeLikelihoods(String sample, Allele A, Allele B, double log10aaPosteriorLikelihoods, double log10abPosteriorLikelihoods, double log10bbPosteriorLikelihoods,
			HashSet<String> cytosineContexts, int numOfCReadsInBisulfiteCStrand, int numOfTReadsInBisulfiteCStrand, int numOfOtherReadsInBisulfiteCStrand, int numOfGReadsInGenotypeGStrand,
			int numOfAReadsInGenotypeGStrand, int numOfOtherReadsInGenotypeGStrand, int totalDepth, HashMap<String, CytosineParameters> cytosineParameters, String bestMatchedCytosinePattern,
			HashMap<Integer, double[]> GPsBeforeCytosineTenGenotypes, HashMap<Integer, double[]> GPsAfterCytosineTenGenotypes, HashMap<String, Double[]> GPsAtCytosineTenGenotypes,
			ReadBackedPileup pileup, BisulfiteArgumentCollection BAC) {

		this.sample = sample;
		this.numOfCReadsInBisulfiteCStrand = numOfCReadsInBisulfiteCStrand;
		this.numOfTReadsInBisulfiteCStrand = numOfTReadsInBisulfiteCStrand;
		this.numOfOtherReadsInBisulfiteCStrand = numOfOtherReadsInBisulfiteCStrand;
		this.numOfGReadsInGenotypeGStrand = numOfGReadsInGenotypeGStrand;
		this.numOfAReadsInGenotypeGStrand = numOfAReadsInGenotypeGStrand;
		this.numOfOtherReadsInGenotypeGStrand = numOfOtherReadsInGenotypeGStrand;
		this.totalDepth = totalDepth;
		this.bestMatchedCytosinePattern = bestMatchedCytosinePattern;
		this.cytosineParameters = cytosineParameters;

		this.cytosineContexts = cytosineContexts;

		this.GPsBeforeCytosineTenGenotypes = GPsBeforeCytosineTenGenotypes;
		this.GPsAfterCytosineTenGenotypes = GPsAfterCytosineTenGenotypes;
		this.GPsAtCytosineTenGenotypes = GPsAtCytosineTenGenotypes;

		this.GPsAtCytosineNormalizedByThreeGenotypes = new double[] { log10aaPosteriorLikelihoods, log10abPosteriorLikelihoods, log10bbPosteriorLikelihoods };
		this.A = A;
		this.B = B;
		this.pileup = pileup;
		this.BAC = BAC;
		setupSupportingReadsInfo();
	}

	public double getAALikelihoods() {
		return GPsAtCytosineNormalizedByThreeGenotypes[0];
	}

	public double getABLikelihoods() {
		return GPsAtCytosineNormalizedByThreeGenotypes[1];
	}

	public Allele getAlleleA() {
		return A;
	}

	public Allele getAlleleB() {
		return B;
	}

	public String getAveBaseQualAsString() {
		String aveBaseQual = null;
		if (!Double.isNaN(totalBaseQualAlleleA / readsAlleleA) && !Double.isNaN(totalBaseQualAlleleA / readsAlleleA)) {
			aveBaseQual = String.format("%.1f", totalBaseQualAlleleA / readsAlleleA) + "," + String.format("%.1f", totalBaseQualAlleleB / readsAlleleB);
		} else if (!Double.isNaN(totalBaseQualAlleleA / readsAlleleA)) {
			aveBaseQual = String.format("%.1f", totalBaseQualAlleleA / readsAlleleA);
		} else if (!Double.isNaN(totalBaseQualAlleleB / readsAlleleB)) {
			aveBaseQual = String.format("%.1f", totalBaseQualAlleleB / readsAlleleB);
		} else {
			aveBaseQual = VCFConstants.MISSING_VALUE_v4;
		}
		return aveBaseQual;
	}

	public String getBaseCountStatusAsString() {

		String cStatus = numOfCReadsInBisulfiteCStrand + "," + numOfTReadsInBisulfiteCStrand + "," + numOfOtherReadsInBisulfiteCStrand + "," + numOfGReadsInGenotypeGStrand + ","
				+ numOfAReadsInGenotypeGStrand + "," + numOfOtherReadsInGenotypeGStrand;
		return cStatus;
	}

	public double getBBLikelihoods() {
		return GPsAtCytosineNormalizedByThreeGenotypes[2];
	}

	public String getBestMatchedCytosinePattern() {
		return bestMatchedCytosinePattern;
	}

	// Average base quality for reads supporting alleles. For each allele, in
	// the same order as listed
	public double getBQ() {
		return rmsBaseQual;
	}

	public HashMap<String, CytosineParameters> getCytosineParameters() {
		return cytosineParameters;
	}

	public int getDepth() {
		return totalDepth;
	}

	// Reads supporting ALT. Number of 1) forward ref alleles; 2) reverse ref;
	// 3) forward non-ref; 4) reverse non-ref alleles"
	public String getDP4AsString() {
		String DP4String = readsFwdRef + "," + readsRevRef + "," + readsFwdAlt + "," + readsRevAlt;
		return DP4String;
	}

	public int getGoodReadsCoverage() {
		return goodReadsCov;
	}

	public HashMap<Integer, double[]> getGPsAfterCytosineTenGenotypes() {
		return GPsAfterCytosineTenGenotypes;
	}

	public HashMap<String, Double[]> getGPsAtCytosineTenGenotypes() {
		return GPsAtCytosineTenGenotypes;
	}

	public HashMap<Integer, double[]> getGPsBeforeCytosineTenGenotypes() {
		return GPsBeforeCytosineTenGenotypes;
	}

	public double[] getLikelihoods() {
		return GPsAtCytosineNormalizedByThreeGenotypes;
	}

	public double getMethylationLevel() {
		double methy = (double) numOfCReadsInBisulfiteCStrand / (double) (numOfCReadsInBisulfiteCStrand + numOfTReadsInBisulfiteCStrand);
		//if(BAC.nonDirectional){
		//	methy = (double) (numOfCReadsInBisulfiteCStrand + numOfGReadsInGenotypeGStrand)/ (double) (numOfCReadsInBisulfiteCStrand + numOfTReadsInBisulfiteCStrand + numOfGReadsInGenotypeGStrand + numOfAReadsInGenotypeGStrand);
		//}
			
		return methy;
	}

	public double getMQ() {
		return rmsMapQual;
	}

	public int getMQ0() {
		return mapQual0;
	}

	public int getNumOfAReadsInGenotypeGStrand() {
		return numOfAReadsInGenotypeGStrand;
	}

	public int getNumOfCReadsInBisulfiteCStrand() {
		return numOfCReadsInBisulfiteCStrand;
	}

	public int getNumOfGReadsInGenotypeGStrand() {
		return numOfGReadsInGenotypeGStrand;
	}

	public int getNumOfOtherReadsInBisulfiteCStrand() {
		return numOfOtherReadsInBisulfiteCStrand;
	}

	public int getNumOfOtherReadsInGenotypeGStrand() {
		return numOfOtherReadsInGenotypeGStrand;
	}

	public int getNumOfTReadsInBisulfiteCStrand() {
		return numOfTReadsInBisulfiteCStrand;
	}

	public String getSample() {
		return sample;
	}

	public void setGPsAfterCytosineTenGenotypes(HashMap<Integer, double[]> GPsAfterCytosineTenGenotypes) {
		this.GPsAfterCytosineTenGenotypes = GPsAfterCytosineTenGenotypes;
	}

	public void setGPsBeforeCytosineTenGenotypes(HashMap<Integer, double[]> GPsBeforeCytosineTenGenotypes) {
		this.GPsBeforeCytosineTenGenotypes = GPsBeforeCytosineTenGenotypes;
	}

	private void setupSupportingReadsInfo() {
		int[] rmsBaseQualTotal = new int[pileup.depthOfCoverage()];
		int[] rmsMapQualTotal = new int[pileup.depthOfCoverage()];
		int index = 0;

		for (PileupElement p : pileup) {
			SAMRecord samRecord = p.getRead();

			int offset = p.getOffset();
			if (offset < 0)// is deletion
				continue;
			boolean paired = samRecord.getReadPairedFlag();
			

			boolean negStrand = samRecord.getReadNegativeStrandFlag();
			if(paired && samRecord.getSecondOfPairFlag())
				negStrand = !negStrand;
			// Reads supporting ALT. Number of 1) forward ref alleles; 2)
			// reverse ref; 3) forward non-ref; 4) reverse non-ref alleles"
			// Average base quality for reads supporting alleles. For each
			// allele, in the same order as listed
			rmsMapQualTotal[index] = p.getMappingQual();
			rmsBaseQualTotal[index] = p.getQual();
			index++;
			//GATKSAMRecordFilterStorage GATKrecordFilterStor = new GATKSAMRecordFilterStorage(p.getRead(), BAC, ref, p.getOffset());

			if (p.getMappingQual() == 0) {

				mapQual0++;

			}
			if (BisSNPUtils.goodBaseInPileupElement(p, BAC, ref)) {
				// should different by GPs
				goodReadsCov++;
				if (negStrand) {

					if (BaseUtilsMore.iupacCodeEqual(p.getBase(), A.getBases()[0], negStrand)) {
						readsAlleleA++;

						if (Double.isNaN(totalBaseQualAlleleA)) {
							totalBaseQualAlleleA = p.getQual();
						} else {
							totalBaseQualAlleleA += p.getQual();
						}
						if (A.isReference()) {
							readsRevRef++;
						} else {
							readsRevAlt++;
						}
					} else if (BaseUtilsMore.iupacCodeEqual(p.getBase(), B.getBases()[0], negStrand)) {
						readsAlleleB++;

						if (Double.isNaN(totalBaseQualAlleleB)) {
							totalBaseQualAlleleB = p.getQual();
						} else {
							totalBaseQualAlleleB += p.getQual();
						}
						if (B.isReference()) {
							readsRevRef++;
						} else {
							readsRevAlt++;
						}
					}

				} else {
					if (BaseUtilsMore.iupacCodeEqual(p.getBase(), A.getBases()[0], negStrand)) {
						readsAlleleA++;

						if (Double.isNaN(totalBaseQualAlleleA)) {
							totalBaseQualAlleleA = p.getQual();
						} else {
							totalBaseQualAlleleA += p.getQual();
						}
						if (A.isReference()) {
							readsFwdRef++;
						} else {
							readsFwdAlt++;
						}
					} else if (BaseUtilsMore.iupacCodeEqual(p.getBase(), B.getBases()[0], negStrand)) {
						readsAlleleB++;

						if (Double.isNaN(totalBaseQualAlleleB)) {
							totalBaseQualAlleleB = p.getQual();
						} else {
							totalBaseQualAlleleB += p.getQual();
						}
						if (B.isReference()) {
							readsFwdRef++;
						} else {
							readsFwdAlt++;
						}
					}
				}
			}

		}

		rmsBaseQual = MathUtils.rms(rmsBaseQualTotal);
		rmsMapQual = MathUtils.rms(rmsMapQualTotal);
	}

}
