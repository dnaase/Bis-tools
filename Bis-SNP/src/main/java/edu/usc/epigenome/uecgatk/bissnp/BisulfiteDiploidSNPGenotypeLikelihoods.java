package main.java.edu.usc.epigenome.uecgatk.bissnp;

import static java.lang.Math.log10;
import static java.lang.Math.pow;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;

import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.genotyper.DiploidGenotype;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.exceptions.UserException.MalformedBAM;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;



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

public class BisulfiteDiploidSNPGenotypeLikelihoods implements Cloneable {
	// record the reference genome error
	private double PROB_OF_REFERENCE_ERROR = 1e-6;

	public boolean VERBOSE = false;
	private BisulfiteArgumentCollection BAC;
	private double BISULFITE_CONVERSION_RATE;

	private double DETERMINED_CYTOSINE_TYPE_C_METHY_NEG = 0;
	// record cytosine pattern methylation level in positive strand which owns the
	// maximum posterior probability
	private double DETERMINED_CYTOSINE_TYPE_C_METHY_POS = 0;

	private double[] log10Likelihoods = null;

	private double[] log10Posteriors = null;

	private double OVER_CONVERSION_RATE;
	private BisulfiteDiploidSNPGenotypePriors priors = null;

	private ReferenceContext ref;

	public BisulfiteDiploidSNPGenotypeLikelihoods(RefMetaDataTracker tracker, ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors, BisulfiteArgumentCollection BAC,
			double cytosineMethyStatus) {
		this.priors = priors;
		this.BAC = BAC;
		this.ref = ref;
		this.BISULFITE_CONVERSION_RATE = BAC.bsRate;
		this.OVER_CONVERSION_RATE = BAC.overRate;
		this.PROB_OF_REFERENCE_ERROR = BAC.referenceGenomeErr;
		this.DETERMINED_CYTOSINE_TYPE_C_METHY_POS = cytosineMethyStatus;
		this.DETERMINED_CYTOSINE_TYPE_C_METHY_NEG = cytosineMethyStatus;

	}


	/**
	 * add pileupelement into likelihood calculation. and call likelihood calculation function
	 */

	public int add(PileupElement p, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, boolean negStrand) {

		byte observedBase = 0, qualityScore = 0;
		
		observedBase = p.getBase();
		SAMRecord samRecord = p.getRead();
		boolean paired = samRecord.getReadPairedFlag();
		if (paired && samRecord.getSecondOfPairFlag() && !BAC.nonDirectional) {
			negStrand = !negStrand;

		}
		//synchronized(p.getRead()) {
		//	if(!BisSNPUtils.goodPileupElement(p, BAC, ref))
		//		return 0;
		//}
		if(!BisSNPUtils.goodBaseInPileupElement(p, BAC, ref))
			return 0;

		byte qual = p.getQual();
		if (qual > SAMUtils.MAX_PHRED_SCORE)
			throw new MalformedBAM(
					p.getRead(),
					String.format(
							"the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details",
							SAMUtils.MAX_PHRED_SCORE, qual, p.getRead().getReadName()));
		if (capBaseQualsAtMappingQual)
			qual = (byte) Math.min(p.getQual(), p.getMappingQual());

		qualityScore = qual;



		BisulfiteDiploidSNPGenotypeLikelihoods gl;
		if(BAC.nonDirectional){
			gl = getCalculateGenotypeLikelihoodsNonDirectional(observedBase, qualityScore, negStrand);
		}
		else{
			gl = getCalculateGenotypeLikelihoods(observedBase, qualityScore, negStrand);
		}
		

		// for bad bases, there are no likelihoods
		if (gl == null)
			return 0;

		double[] likelihoods = gl.getLikelihoods();

		for (DiploidGenotype g : DiploidGenotype.values()) {
			double likelihood = likelihoods[g.ordinal()];

			log10Likelihoods[g.ordinal()] += likelihood;
			log10Posteriors[g.ordinal()] += likelihood;
		}

		return 1;
	}

	/**
	 * add pileup into likelihood calculation.
	 */
	public int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual) {
		int n = 0;
		

		ReadBackedPileup pileupPositiveStrand = pileup.getPositiveStrandPileup();
		for (PileupElement p : pileupPositiveStrand)
			n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, false);

		ReadBackedPileup pileupNegativeStrand = pileup.getNegativeStrandPileup();
		for (PileupElement p : pileupNegativeStrand)
			n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, true);

		if (VERBOSE) {
			System.out.println("summary:");
			for (DiploidGenotype g : DiploidGenotype.values()) {
				System.out.printf("%s\t", g);
			}
			System.out.println();
			for (DiploidGenotype g : DiploidGenotype.values()) {
				System.out.printf("%.2f\t", log10Likelihoods[g.ordinal()]);
			}
			System.out.println();
			for (DiploidGenotype g : DiploidGenotype.values()) {
				System.out.printf("%.2f\t", getPriors()[g.ordinal()]);
			}
			System.out.println();
		}

		return n;
	}

	/**
	 * clear the likelihood to primary status. this could save time on prior calculation.
	 */
	public void clearLikelihoods(double cytosineMethyStatus) {

		this.DETERMINED_CYTOSINE_TYPE_C_METHY_POS = cytosineMethyStatus;
		this.DETERMINED_CYTOSINE_TYPE_C_METHY_NEG = cytosineMethyStatus;
		setToZeroBs();
	}

	public double[] getLikelihoods() {
		return log10Likelihoods;
	}

	public double[] getPosteriors() {
		return log10Posteriors;
	}

	public double[] getPriors() {
		return this.priors.getPriors();
	}

	/**
	 * set prior calculated
	 */
	public void setPriors(RefMetaDataTracker tracker, ReferenceContext ref, double heterozygousity, double novelDbsnpHet, double validateDbsnpHet, GenomeLoc loc) {
		this.priors.setPriors(tracker, ref, BAC, loc);
		setToZeroBs();
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		BisulfiteDiploidSNPGenotypeLikelihoods c = (BisulfiteDiploidSNPGenotypeLikelihoods) super.clone();
		c.priors = priors;
		c.log10Likelihoods = log10Likelihoods.clone();
		c.log10Posteriors = log10Posteriors.clone();
		return c;
	}

	/**
	 * the main likelihood calculation function for any input base with quality information and
	 * strand information
	 */
	protected BisulfiteDiploidSNPGenotypeLikelihoods getCalculateGenotypeLikelihoods(byte observedBase, byte qualityScore, boolean negStrand) {

		try {
			BisulfiteDiploidSNPGenotypeLikelihoods gl = (BisulfiteDiploidSNPGenotypeLikelihoods) this.clone();
			gl.setToZeroBs();
			double error = pow(10, -qualityScore / 10.0);
			double pOfBase1 = 0.5, pOfBase2 = 1.0 - pOfBase1;
			for (DiploidGenotype g : DiploidGenotype.values()) {
				double likelihood = 0.0;
				if (!negStrand) {
					if (BaseUtils.basesAreEqual(observedBase, BaseUtils.Base.C.base)) {
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0 - error)
								* ((1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * (1.0 - BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * (1.0 - OVER_CONVERSION_RATE)) : pOfBase1
								* (error / 3.0);
						likelihood += observedBase == g.base2 ? pOfBase2 * (1.0 - error)
								* ((1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * (1.0 - BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * (1.0 - OVER_CONVERSION_RATE)) : pOfBase2
								* (error / 3.0);
						if (VERBOSE) {
							// System.out.println("flag1: observedBase-" +
							// observedBase + "\t" + "g.base1-" + g.base1 + "\t"
							// + "g.base2-" + g.base2 + "\t" + "likelihood-" +
							// log10(likelihood));
						}
					} else if (BaseUtils.basesAreEqual(observedBase, BaseUtils.Base.T.base)) {
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0 - error) : (g.base1 == BaseUtils.Base.C.base ? pOfBase1
								* (error / 3.0 + (1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * OVER_CONVERSION_RATE) : pOfBase1
								* (error / 3.0));
						likelihood += observedBase == g.base2 ? pOfBase1 * (1.0 - error) : (g.base2 == BaseUtils.Base.C.base ? pOfBase2
								* (error / 3.0 + (1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * OVER_CONVERSION_RATE) : pOfBase2
								* (error / 3.0));

						if (VERBOSE) {
							// System.out.println("flag3: observedBase-" +
							// observedBase + "\t" + "g.base1-" + g.base1 + "\t"
							// + "g.base2-" + g.base2 + "\t" + "likelihood-" +
							// log10(likelihood));
						}
					} else {
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0 - error) : pOfBase1 * (error / 3.0);
						likelihood += observedBase == g.base2 ? pOfBase2 * (1.0 - error) : pOfBase2 * (error / 3.0);
						if (VERBOSE) {
							// System.out.println("flag5: observedBase-" +
							// observedBase + "\t" + "g.base1-" + g.base1 + "\t"
							// + "g.base2-" + g.base2 + "\t" + "likelihood-" +
							// log10(likelihood));
						}
					}
				} else {
					if (BaseUtils.basesAreEqual(observedBase, BaseUtils.Base.G.base)) {
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0 - error)
								* ((1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * (1.0 - BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * (1.0 - OVER_CONVERSION_RATE)) : pOfBase1
								* (error / 3.0);
						likelihood += observedBase == g.base2 ? pOfBase2 * (1.0 - error)
								* ((1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * (1.0 - BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * (1.0 - OVER_CONVERSION_RATE)) : pOfBase2
								* (error / 3.0);
						if (VERBOSE) {
							// System.out.println("flag6: observedBase-" +
							// observedBase + "\t" + "g.base1-" + g.base1 + "\t"
							// + "g.base2-" + g.base2 + "\t" + "likelihood-" +
							// log10(likelihood));
						}

					} else if (BaseUtils.basesAreEqual(observedBase, BaseUtils.Base.A.base)) {
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0 - error) : (g.base1 == BaseUtils.Base.G.base ? pOfBase1
								* (error / 3.0 + (1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * OVER_CONVERSION_RATE) : pOfBase1
								* (error / 3.0));
						likelihood += observedBase == g.base2 ? pOfBase1 * (1.0 - error) : (g.base2 == BaseUtils.Base.G.base ? pOfBase2
								* (error / 3.0 + (1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * OVER_CONVERSION_RATE) : pOfBase2
								* (error / 3.0));
						if (VERBOSE) {
							// System.out.println("flag9: observedBase-" +
							// observedBase + "\t" + "g.base1-" + g.base1 + "\t"
							// + "g.base2-" + g.base2 + "\t" + "likelihood-" +
							// log10(likelihood));
						}
					} else {
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0 - error) : pOfBase1 * (error / 3.0);
						likelihood += observedBase == g.base2 ? pOfBase2 * (1.0 - error) : pOfBase2 * (error / 3.0);
						if (VERBOSE) {
							// System.out.println("flag9: observedBase-" +
							// observedBase + "\t" + "g.base1-" + g.base1 + "\t"
							// + "g.base2-" + g.base2 + "\t" + "likelihood-" +
							// log10(likelihood) + "\t" + qualityScore + "\t" +
							// error);
						}
					}
				}
				double logLikelihood = log10(likelihood);
				gl.log10Likelihoods[g.ordinal()] += logLikelihood;
				gl.log10Posteriors[g.ordinal()] += logLikelihood;
			}

			if (VERBOSE) {
				for (DiploidGenotype g : DiploidGenotype.values()) {
					System.out.printf("%s\t", g);
				}
				System.out.println();
				for (DiploidGenotype g : DiploidGenotype.values()) {
					System.out.printf("%.2f\t", gl.log10Likelihoods[g.ordinal()]);
				}
				System.out.println();
				for (DiploidGenotype g : DiploidGenotype.values()) {
					System.out.printf("%.2f\t", gl.getPriors()[g.ordinal()]);
				}
				System.out.println();
			}

			return gl;

		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException(e);
		}

	}

	/**
	 * the likelihood calculation function for NON-direcitonal protocol only !!! 
	 * any input base with quality information and strand information
	 */
	protected BisulfiteDiploidSNPGenotypeLikelihoods getCalculateGenotypeLikelihoodsNonDirectional(byte observedBase, byte qualityScore, boolean negStrand) {

		try {
			BisulfiteDiploidSNPGenotypeLikelihoods gl = (BisulfiteDiploidSNPGenotypeLikelihoods) this.clone();
			gl.setToZeroBs();
			double error = pow(10, -qualityScore / 10.0);
			double pOfBase1 = 0.5, pOfBase2 = 1.0 - pOfBase1;
			for (DiploidGenotype g : DiploidGenotype.values()) {
				double likelihood = 0.0;
				if (!negStrand) {
					if (BaseUtils.basesAreEqual(observedBase, BaseUtils.Base.C.base) || BaseUtils.basesAreEqual(observedBase, BaseUtils.Base.G.base)) {
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0 - error)
								* ((1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * (1.0 - BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * (1.0 - OVER_CONVERSION_RATE)) : pOfBase1
								* (error / 3.0);
						likelihood += observedBase == g.base2 ? pOfBase2 * (1.0 - error)
								* ((1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * (1.0 - BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * (1.0 - OVER_CONVERSION_RATE)) : pOfBase2
								* (error / 3.0);
						if (VERBOSE) {
							// System.out.println("flag1: observedBase-" +
							// observedBase + "\t" + "g.base1-" + g.base1 + "\t"
							// + "g.base2-" + g.base2 + "\t" + "likelihood-" +
							// log10(likelihood));
						}
					} else if (BaseUtils.basesAreEqual(observedBase, BaseUtils.Base.T.base) || BaseUtils.basesAreEqual(observedBase, BaseUtils.Base.A.base)) {
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0 - error) : (g.base1 == BaseUtils.Base.C.base || g.base1 == BaseUtils.Base.G.base ? pOfBase1
								* (error / 3.0 + (1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * OVER_CONVERSION_RATE) : pOfBase1
								* (error / 3.0));
						likelihood += observedBase == g.base2 ? pOfBase1 * (1.0 - error) : (g.base2 == BaseUtils.Base.C.base || g.base2 == BaseUtils.Base.G.base? pOfBase2
								* (error / 3.0 + (1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_POS) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_POS * OVER_CONVERSION_RATE) : pOfBase2
								* (error / 3.0));

						if (VERBOSE) {
							// System.out.println("flag3: observedBase-" +
							// observedBase + "\t" + "g.base1-" + g.base1 + "\t"
							// + "g.base2-" + g.base2 + "\t" + "likelihood-" +
							// log10(likelihood));
						}
					} else {
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0 - error) : pOfBase1 * (error / 3.0);
						likelihood += observedBase == g.base2 ? pOfBase2 * (1.0 - error) : pOfBase2 * (error / 3.0);
						if (VERBOSE) {
							// System.out.println("flag5: observedBase-" +
							// observedBase + "\t" + "g.base1-" + g.base1 + "\t"
							// + "g.base2-" + g.base2 + "\t" + "likelihood-" +
							// log10(likelihood));
						}
					}
				} else {
					if (BaseUtils.basesAreEqual(observedBase, BaseUtils.Base.G.base) || BaseUtils.basesAreEqual(observedBase, BaseUtils.Base.C.base)) {
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0 - error)
								* ((1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * (1.0 - BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * (1.0 - OVER_CONVERSION_RATE)) : pOfBase1
								* (error / 3.0);
						likelihood += observedBase == g.base2 ? pOfBase2 * (1.0 - error)
								* ((1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * (1.0 - BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * (1.0 - OVER_CONVERSION_RATE)) : pOfBase2
								* (error / 3.0);
						if (VERBOSE) {
							// System.out.println("flag6: observedBase-" +
							// observedBase + "\t" + "g.base1-" + g.base1 + "\t"
							// + "g.base2-" + g.base2 + "\t" + "likelihood-" +
							// log10(likelihood));
						}

					} else if (BaseUtils.basesAreEqual(observedBase, BaseUtils.Base.A.base) || BaseUtils.basesAreEqual(observedBase, BaseUtils.Base.T.base)) {
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0 - error) : (g.base1 == BaseUtils.Base.G.base || g.base1 == BaseUtils.Base.C.base ? pOfBase1
								* (error / 3.0 + (1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * OVER_CONVERSION_RATE) : pOfBase1
								* (error / 3.0));
						likelihood += observedBase == g.base2 ? pOfBase1 * (1.0 - error) : (g.base2 == BaseUtils.Base.G.base || g.base2 == BaseUtils.Base.C.base ? pOfBase2
								* (error / 3.0 + (1.0 - DETERMINED_CYTOSINE_TYPE_C_METHY_NEG) * BISULFITE_CONVERSION_RATE + DETERMINED_CYTOSINE_TYPE_C_METHY_NEG * OVER_CONVERSION_RATE) : pOfBase2
								* (error / 3.0));
						if (VERBOSE) {
							// System.out.println("flag9: observedBase-" +
							// observedBase + "\t" + "g.base1-" + g.base1 + "\t"
							// + "g.base2-" + g.base2 + "\t" + "likelihood-" +
							// log10(likelihood));
						}
					} else {
						likelihood += observedBase == g.base1 ? pOfBase1 * (1.0 - error) : pOfBase1 * (error / 3.0);
						likelihood += observedBase == g.base2 ? pOfBase2 * (1.0 - error) : pOfBase2 * (error / 3.0);
						if (VERBOSE) {
							// System.out.println("flag9: observedBase-" +
							// observedBase + "\t" + "g.base1-" + g.base1 + "\t"
							// + "g.base2-" + g.base2 + "\t" + "likelihood-" +
							// log10(likelihood) + "\t" + qualityScore + "\t" +
							// error);
						}
					}
				}
				double logLikelihood = log10(likelihood);
				gl.log10Likelihoods[g.ordinal()] += logLikelihood;
				gl.log10Posteriors[g.ordinal()] += logLikelihood;
			}

			if (VERBOSE) {
				for (DiploidGenotype g : DiploidGenotype.values()) {
					System.out.printf("%s\t", g);
				}
				System.out.println();
				for (DiploidGenotype g : DiploidGenotype.values()) {
					System.out.printf("%.2f\t", gl.log10Likelihoods[g.ordinal()]);
				}
				System.out.println();
				for (DiploidGenotype g : DiploidGenotype.values()) {
					System.out.printf("%.2f\t", gl.getPriors()[g.ordinal()]);
				}
				System.out.println();
			}

			return gl;

		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException(e);
		}

	}
	
	protected void setToZeroBs() {
		log10Likelihoods = new double[DiploidGenotype.values().length];
		for (DiploidGenotype g : DiploidGenotype.values()) {
			int i = g.ordinal();
			log10Likelihoods[i] = log10(1 - PROB_OF_REFERENCE_ERROR);
		}
		// posteriors are all set to the priors
		log10Posteriors = this.priors.getPriors().clone();
	}

	

}
