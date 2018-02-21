package main.java.edu.usc.epigenome.uecgatk.bissnp;

import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.refdata.utils.GATKFeature;
import org.broadinstitute.gatk.utils.refdata.utils.RODRecordList;
import org.broadinstitute.gatk.utils.genotyper.DiploidGenotype;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypePriors;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.MathUtils;

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

public class BisulfiteDiploidSNPGenotypePriors implements GenotypePriors {
	public static double DBSNP_NOVAL_HETEROZYGOSITY = 0.02;

	public static double DBSNP_VALIDATE_HETEROZYGOSITY = 0.1;
	
	public static double GLOBAL_MINOR_ALLELE_FREQUECNY = 0.05;

	// --------------------------------------------------------------------------------------------------------------
	//
	// Constants and static information
	//
	// --------------------------------------------------------------------------------------------------------------
	public static double HETEROZYGOSITY = 1e-3;
	public static double TRANSITION_VS_TRANSVERSION = 2.0;

	protected static double BISULFITE_CONVERSION_RATE = 0;
	protected static double OVER_CONVERSION_RATE = 0;

	private final static double[] flatPriors = new double[DiploidGenotype.values().length];

	private double[] priors = null;

	/*
	 * made flat prior for each genotype
	 */
	static {
		for (DiploidGenotype g : DiploidGenotype.values()) {
			flatPriors[g.ordinal()] = Math.log10(1.0 / DiploidGenotype.values().length);
		}
	}

	/**
	 * Create a new DiploidGenotypePriors object with flat priors for each diploid genotype
	 */
	public BisulfiteDiploidSNPGenotypePriors() {
		priors = flatPriors.clone();
	}

	/**
	 * Create a new GenotypeLikelihoods object with priors for a diploid with heterozygosity and
	 * reference base ref
	 * 
	 * @param ref
	 * @param heterozygosity
	 * @param probOfTriStateGenotype
	 *            The prob of seeing a true B/C het when the reference is A
	 */
	public BisulfiteDiploidSNPGenotypePriors(byte ref, double heterozygosity, double probOfTriStateGenotype) {

		priors = flatPriors.clone();
		HETEROZYGOSITY = heterozygosity;
	}

	public BisulfiteDiploidSNPGenotypePriors(byte ref, double heterozygosity, double probOfTriStateGenotype, double bisulfiteConversionRate) {
		priors = flatPriors.clone();
		HETEROZYGOSITY = heterozygosity;

	}
	

	public BisulfiteDiploidSNPGenotypePriors(double[] log10Priors) {
		priors = log10Priors.clone();
	}

	@Override
	public double getHeterozygosity() {
		// TODO Auto-generated method stub
		return HETEROZYGOSITY;
	}

	/**
	 * Returns the prior associated with DiploidGenotype g
	 * 
	 * @param g
	 * @return log10 prior as a double
	 */
	public double getPrior(DiploidGenotype g) {
		return getPriors()[g.ordinal()];
	}

	@Override
	public double[] getPriors() {
		return priors;
	}

	/**
	 * Takes reference base, and three priors for hom-ref, het, hom-var, and fills in the priors
	 * vector appropriately.
	 * 
	 * Suppose A is the reference base, and we are given the probability of being hom-ref, het, and
	 * hom-var, and that pTriSateGenotype (product of reference genome error and heterozygousity) is
	 * the true probability of observing reference A and a true genotype of B/C then this sets the
	 * priors to:
	 * 
	 * AA = pHomRef AG = (pHet - pTriStateGenotype) * transitionRate AC = AT = (pHet -
	 * pTriStateGenotype) * transversionRate GG = pHomVar * transitionRate CC = TT = pHomVar *
	 * transversionRate CT = pTriStateGenotype * transversionRate CG = GT = pTriStateGenotype *
	 * transitionRate When i use switch clause, the program turns to be very slow.,so even this code
	 * is ugly for a lot of if, else, but it is much quicker...
	 * 
	 * @param ref
	 * @param heterozyosity
	 * @param pRefError
	 */

	public double[] getReferencePolarizedPriors(byte ref, double heterozyosity, double pRefError) {
		if (!MathUtils.isBounded(pRefError, 0.0, 0.01)) {
			throw new RuntimeException(String.format("BUG: p Reference error is out of bounds (0.0 - 0.01) is allow range %f", pRefError));
		}

		double pTriStateGenotype = heterozyosity * pRefError;

		double pHomRef = heterozygosity2HomRefProbability(heterozyosity);
		double pHet = heterozygosity2HetProbability(heterozyosity);
		double pHomVar = heterozygosity2HomVarProbability(heterozyosity);

		if (MathUtils.compareDoubles(pHomRef + pHet + pHomVar, 1.0) != 0) {
			throw new RuntimeException(String.format("BUG: Prior probabilities don't sum to one => %f, %f, %f", pHomRef, pHet, pHomVar));
		}

		double[] priors = new double[DiploidGenotype.values().length];

		for (DiploidGenotype g : DiploidGenotype.values()) {
			double POfG = 1;

			// final double transitionRate = 2.0/3.0;
			// final double transversionRate = 1.0/6.0;

			final double transversionRate = 1.0 / (2*(TRANSITION_VS_TRANSVERSION + 1));
			final double transitionRate = transversionRate * 2 * TRANSITION_VS_TRANSVERSION;

			if (BaseUtils.basesAreEqual(ref, BaseUtils.Base.A.base)) { // **ref = A
				if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.A.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.A.base)) {
						POfG = pHomRef; // genotype: AA
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.G.base)) {
						POfG = (pHet - pTriStateGenotype) * transitionRate; // genotype:
																			// AG
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.C.base)) {
						POfG = (pHet - pTriStateGenotype) * transversionRate; // genotype:
																				// AC
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = ((pHet - pTriStateGenotype) * transversionRate); // genotype:
																				// AT
					} else {

					}
				} else if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.G.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.G.base)) {
						POfG = pHomVar * transitionRate; // genotype: GG
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = pTriStateGenotype * transitionRate; // genotype:
																	// GT
					} else {

					}
				} else if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.C.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.G.base)) {
						POfG = pTriStateGenotype * transitionRate; // genotype:
																	// CG
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.C.base)) {
						POfG = pHomVar * transversionRate; // genotype: CC
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = pTriStateGenotype * transversionRate; // genotype:
																		// CT
					} else {

					}
				} else if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.T.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = pHomVar * transversionRate; // genotype: TT
					} else {

					}
				} else {

				}
			} else if (BaseUtils.basesAreEqual(ref, BaseUtils.Base.G.base)) { // **ref = G
				if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.A.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.A.base)) {
						POfG = pHomVar * transitionRate;
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.G.base)) {
						POfG = (pHet - pTriStateGenotype) * transitionRate;
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.C.base)) {
						POfG = (pTriStateGenotype * transitionRate);
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = pTriStateGenotype * transitionRate;
					} else {

					}
				} else if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.G.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.G.base)) {
						POfG = pHomRef;
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = ((pHet - pTriStateGenotype) * transversionRate);
					} else {

					}
				} else if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.C.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.G.base)) {
						POfG = ((pHet - pTriStateGenotype) * transversionRate);
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.C.base)) {
						POfG = (pHomVar * transversionRate);
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = (pTriStateGenotype * transversionRate);
					} else {

					}
				} else if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.T.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = pHomVar * transversionRate;
					} else {

					}
				} else {

				}
			} else if (BaseUtils.basesAreEqual(ref, BaseUtils.Base.C.base)) { // **ref = C
				if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.A.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.A.base)) {
						POfG = pHomVar * transversionRate;
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.G.base)) {
						POfG = pTriStateGenotype * transversionRate;
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.C.base)) {
						POfG = ((pHet - pTriStateGenotype) * transversionRate);
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = pTriStateGenotype * transitionRate;
					} else {

					}
				} else if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.G.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.G.base)) {
						POfG = pHomVar * transversionRate;
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = pTriStateGenotype * transitionRate;
					} else {

					}
				} else if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.C.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.G.base)) {
						POfG = ((pHet - pTriStateGenotype) * transversionRate);
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.C.base)) {
						POfG = pHomRef;
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = ((pHet - pTriStateGenotype) * transitionRate);
					} else {

					}
				} else if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.T.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = pHomVar * transitionRate;
					} else {

					}
				} else {

				}
			} else if (BaseUtils.basesAreEqual(ref, BaseUtils.Base.T.base)) { // **ref = T
				if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.A.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.A.base)) {
						POfG = pHomVar * transversionRate;
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.G.base)) {
						POfG = pTriStateGenotype * transversionRate;
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.C.base)) {
						POfG = (pTriStateGenotype * transitionRate);
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = (pHet - pTriStateGenotype) * transversionRate;
					} else {

					}
				} else if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.G.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.G.base)) {
						POfG = pHomVar * transversionRate;
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = (pHet - pTriStateGenotype) * transversionRate;
					} else {

					}
				} else if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.C.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.G.base)) {
						POfG = (pTriStateGenotype * transitionRate);
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.C.base)) {
						POfG = (pHomVar * transitionRate);
					} else if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = ((pHet - pTriStateGenotype) * transitionRate);
					} else {

					}
				} else if (BaseUtils.basesAreEqual(g.base1, BaseUtils.Base.T.base)) {
					if (BaseUtils.basesAreEqual(g.base2, BaseUtils.Base.T.base)) {
						POfG = pHomRef;
					} else {

					}
				} else {

				}
			} else {

			}

			priors[g.ordinal()] = Math.log10(POfG);
		}

		return priors;
	}
	
	/**
	 * Takes reference base, and three priors for hom-ref, het, hom-var by dbsnp's GMAF, and fills in the priors
	 * vector appropriately.
	 * 
	 * Suppose A is the reference base, and we are given the probability of being hom-ref, het, and
	 * hom-var, and that pTriSateGenotype (product of reference genome error and heterozygousity) is
	 * the true probability of observing reference A and a true genotype of B/C then this sets the
	 * priors to:
	 * 
	 * AA = pHomRef AG = (pHet - pTriStateGenotype) * transitionRate AC = AT = (pHet -
	 * pTriStateGenotype) * transversionRate GG = pHomVar * transitionRate CC = TT = pHomVar *
	 * transversionRate CT = pTriStateGenotype * transversionRate CG = GT = pTriStateGenotype *
	 * transitionRate When i use switch clause, the program turns to be very slow.,so even this code
	 * is ugly for a lot of if, else, but it is much quicker...
	 */

	public double[] getReferencePriorsInvolveDbsnpAf(byte ref, double minorAf, double pRefError, byte alterAllele) {
		if (!MathUtils.isBounded(pRefError, 0.0, 0.01)) {
			throw new RuntimeException(String.format("BUG: p Reference error is out of bounds (0.0 - 0.01) is allow range %f", pRefError));
		}


		minorAf = minorAf-pRefError <= 0 ? minorAf : minorAf-pRefError;
		double majorAf = 1-minorAf;
		

		if (MathUtils.compareDoubles(majorAf + minorAf + pRefError, 1.0) != 0) {
			throw new RuntimeException(String.format("BUG: Prior probabilities don't sum to one => %f, %f, %f", majorAf, minorAf, pRefError));
		}

		double[] priors = new double[DiploidGenotype.values().length];

		for (DiploidGenotype g : DiploidGenotype.values()) {
			double POfG = 1;

			// final double transitionRate = 2.0/3.0;
			// final double transversionRate = 1.0/6.0;

			final double transversionRate = 1.0 / (TRANSITION_VS_TRANSVERSION + 1);
			final double transitionRate = transversionRate * TRANSITION_VS_TRANSVERSION;
			if(BaseUtils.basesAreEqual(g.base1, ref)){ 
				if(BaseUtils.basesAreEqual(g.base2, ref)){ // genotype: ref-ref
					POfG = majorAf * majorAf;
				}
				else if(BaseUtils.basesAreEqual(g.base2, alterAllele)){ // genotype: ref-altMinorAllele
					POfG = majorAf * minorAf;
				}
				else{ // genotype: ref-otherAllele
					if(isTransition(g.base2, ref)){
						POfG = majorAf * pRefError * transitionRate;
					}
					else{
						POfG = majorAf * pRefError * transversionRate;
					}
				}
			}
			else if(BaseUtils.basesAreEqual(g.base1, alterAllele)){
				if(BaseUtils.basesAreEqual(g.base2, ref)){ // genotype: altMinorAllele-ref
					POfG = majorAf * minorAf;
				}
				else if(BaseUtils.basesAreEqual(g.base2, alterAllele)){ // genotype: altMinorAllele-altMinorAllele
					POfG = minorAf * minorAf;;
				}
				else{// genotype: altMinorAllele-otherAllele
					if(isTransition(g.base2, alterAllele)){
						POfG = minorAf * pRefError * transitionRate;
					}
					else{
						POfG = minorAf * pRefError * transversionRate;
					}
				}
			}
			else{
				if(BaseUtils.basesAreEqual(g.base2, ref)){ // genotype: otherAllele-ref
					POfG = majorAf * pRefError;
				}
				else if(BaseUtils.basesAreEqual(g.base2, alterAllele)){  // genotype: otherAllele-altMinorAllele
					POfG = minorAf * pRefError;
				}
				else{ // genotype: otherAllele-otherAllele
					POfG = pRefError*pRefError;
				}
			}

			priors[g.ordinal()] = Math.log10(POfG);
		}

		return priors;
	}

	public double heterozygosity2HetProbability(double h) {
		if (MathUtils.isNegative(h)) {
			throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
		}

		return h;
	}

	/**
	 * 
	 * @param h
	 * @return
	 */
	public double heterozygosity2HomRefProbability(double h) {
		if (MathUtils.isNegative(h)) {
			throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
		}

		double v = 1.0 - (3 * h / 2.0);
		if (MathUtils.isNegative(v)) {
			System.err.println(h);
			throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
		}

		return v;
	}

	public double heterozygosity2HomVarProbability(double h) {
		if (MathUtils.isNegative(h)) {
			throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
		}

		return h / 2.0;
	}

	/**
	 * check if the base is transition or transversion relative to reference base
	 */
	public boolean isTransition(byte base, byte ref) {
		boolean transition;
		switch (ref) {
		case 'A':
			transition = base == BaseUtils.Base.G.base ? true : false;
			break;
		case 'G':
			transition = base == BaseUtils.Base.A.base ? true : false;
			break;
		case 'C':
			transition = base == BaseUtils.Base.T.base ? true : false;
			break;
		case 'T':
			transition = base == BaseUtils.Base.C.base ? true : false;
			break;
		default:
			throw new RuntimeException(String.format("ref base is bad "));
		}
		return transition;
	}

	// set prior for each locus
	public void setPriors(RefMetaDataTracker tracker, ReferenceContext ref, BisulfiteArgumentCollection bac, GenomeLoc loc) {

		DBSNP_NOVAL_HETEROZYGOSITY = bac.novelDbsnpHet;
		DBSNP_VALIDATE_HETEROZYGOSITY = bac.validateDbsnpHet;
		TRANSITION_VS_TRANSVERSION = bac.tiVsTv;
		double probOfTriStateGenotype = bac.referenceGenomeErr;
		double heterozygosity = bac.heterozygosity;
		byte refBase = ref.getBase();

		int count = 0;
		for (RODRecordList rods : tracker.getBoundRodTracks()) {

			for (GATKFeature vc_input : rods) {
				if (vc_input != null && vc_input.getUnderlyingObject() instanceof VariantContext) {
					VariantContext dbsnp = (VariantContext) vc_input.getUnderlyingObject();
					if (dbsnp.isSNP()) {
						count++;
						// VLD is only for vcf4 format, in 3.3 like mouse dbSNP file,
						// it will not works..
						dbsnp.getAltAlleleWithHighestAlleleCount();
						if(dbsnp.getAttribute(BisulfiteVCFConstants.GLOBAL_MINOR_ALLELE_FRE) != null){
							//System.err.println(dbsnp.getID() + "\t" + dbsnp.getAttributeAsDouble(BisulfiteVCFConstants.GLOBAL_MINOR_ALLELE_FRE, GLOBAL_MINOR_ALLELE_FREQUECNY));
							priors = getReferencePriorsInvolveDbsnpAf(refBase, dbsnp.getAttributeAsDouble(BisulfiteVCFConstants.GLOBAL_MINOR_ALLELE_FRE, GLOBAL_MINOR_ALLELE_FREQUECNY), probOfTriStateGenotype, dbsnp.getAltAlleleWithHighestAlleleCount().getBases()[0]);
						}
						else if (dbsnp.getAttribute(BisulfiteVCFConstants.VALIDATED_SNP) != null) {
							
							priors = getReferencePolarizedPriors(refBase, DBSNP_VALIDATE_HETEROZYGOSITY, probOfTriStateGenotype);
						} else {
							priors = getReferencePolarizedPriors(refBase, DBSNP_NOVAL_HETEROZYGOSITY, probOfTriStateGenotype);
						}

						break;
					}

				}
			}
		}

		if (count == 0)
			priors = getReferencePolarizedPriors(refBase, heterozygosity, probOfTriStateGenotype);

	}

	@Override
	public boolean validate(boolean throwException) {
		// TODO Auto-generated method stub
		return false;
	}

}
