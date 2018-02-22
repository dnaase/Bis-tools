package main.java.edu.usc.epigenome.uecgatk.bissnp;

import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.List;
import java.util.HashMap;
import java.util.HashSet;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.Allele;

import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.AlignmentContextUtils;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.genotyper.DiploidGenotype;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;


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

public class BisulfiteSNPGenotypeLikelihoodsCalculationModel {

    protected Byte alternateAllele = null;

    protected Byte bestAllele = null;

    protected int numCNegStrand = 0;
    protected int numCPosStrand = 0;
    protected int numOtherNegStrand = 0;
    protected int numOtherPosStrand = 0;
    protected int numTNegStrand = 0;
    protected int numTPosStrand = 0;
    protected long testLoc;
    private BisulfiteArgumentCollection BAC;
    private HashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs = null;

    private double FLAT_METHY_STATUS = 0.5;
    private Allele refAllele = null;
    private boolean useBAQ = false;

    public BisulfiteSNPGenotypeLikelihoodsCalculationModel(BisulfiteArgumentCollection BAC, boolean useBAQ) {

        this.BAC = BAC;

        this.testLoc = BAC.testLocus;

        this.useBAQ = useBAQ;
        if (BAC.sequencingMode == BisulfiteEnums.MethylSNPModel.NM) {
            FLAT_METHY_STATUS = 0.0;
        }
    }

    public String checkCytosineStatus(ReadBackedPileup pileup, BisulfiteArgumentCollection BAC, RefMetaDataTracker tracker, ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors,
                                      HashMap<Integer, double[]> GPsBeforeCytosineTenGenotypes, HashMap<Integer, double[]> GPsAfterCytosineTenGenotypes, HashMap<String, Double[]> GPsAtCytosineTenGenotypes,
                                      HashMap<String, CytosineParameters> cytosineParametersStatus) {
        String bestCytosinePattern = null;

        GenomeLoc location = pileup.getLocation();
        String contig = location.getContig();
        int position = location.getStart();
        double tmpMethy = FLAT_METHY_STATUS;

        // check adjacent position likelihood
        int maxCytosineLength = BAC.cytosineDefined.getMaxCytosinePatternLen();
        HashMap<Integer, methyStatus> cytosineAndAdjacent = new HashMap<Integer, methyStatus>(maxCytosineLength * 2 + 1);
        byte[] refBytes = new byte[maxCytosineLength * 2 + 1];

        for (int i = 0 - maxCytosineLength, index = 0; i <= maxCytosineLength; i++, index++) {
            GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(contig, position + i);
            if (!ref.getWindow().containsP(loc))
                continue;

            ReferenceContext tmpRef = new ReferenceContext(ref.getGenomeLocParser(), loc, ref.getWindow(), ref.getBases());
            refBytes[index] = tmpRef.getBase();

            if (i == 0)
                continue;
            List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
            ;
            List<Integer> elementOffsets = new ArrayList<Integer>();

            for (PileupElement p : pileup) {
                int elementOffset = i + p.getOffset();
              //  int pos = p.getRead().getReferencePositionAtReadPosition(elementOffset);
                if (elementOffset < 0 || elementOffset > p.getRead().getReadLength() - 1
                    || BisSNPUtils.isInsertionDeletionSoftHardClipByPos(p.getRead(),elementOffset))
                        //|| pos < p.getRead().getAlignmentStart() || pos >= p.getRead().getAlignmentEnd())
                    continue;
                elementOffsets.add(elementOffset);
                reads.add(p.getRead());

               // if(p.getRead().getReadName().equalsIgnoreCase("H32CJADXY160408:2:1216:18985:75061")){
                   // System.err.println(p.getRead().getReadName() + "\t" + p.getRead().getReadNegativeStrandFlag() + "\t" +
                    //        p.getRead().getAlignmentStart()  + "\t" +  p.getRead().getAlignmentEnd()+ "\t" +
                   //        p.getRead().getUnclippedStart()  + "\t" +  p.getRead().getUnclippedEnd() + "\t"
                    //        + elementOffset + "\t" + i + "\t" + p.getOffset() + "\t" + pos + "\t" + p.getRead().getCigarString());
               // }

            }
            ReadBackedPileup tmpPileup = new ReadBackedPileupImpl(loc, reads, elementOffsets);

            tmpMethy = getMethyLevelFromPileup(tmpPileup, ref);
            BisulfiteDiploidSNPGenotypeLikelihoods tmpGL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, tmpRef, priors, BAC, tmpMethy);

            tmpGL.setPriors(tracker, tmpRef, BAC.heterozygosity, BAC.novelDbsnpHet, BAC.validateDbsnpHet, loc);
            if (position == BAC.testLocus) {
                System.err.println("i: " + i + "\ttmpRef: " + tmpRef.getBase());
                tmpGL.VERBOSE = true;
            }
            int nGoodBases = tmpGL.add(tmpPileup, true, true);
            if (nGoodBases == 0)
                continue;
            double[] posteriorNormalized = MathUtils.normalizeFromLog10(tmpGL.getPosteriors(), true, false);
            Integer distanceToCytosine = Math.abs(i);
            if (i < 0) {
                GPsBeforeCytosineTenGenotypes.put(distanceToCytosine, posteriorNormalized.clone());
            } else {
                GPsAfterCytosineTenGenotypes.put(distanceToCytosine, posteriorNormalized.clone());
            }
            getBestGenotypeFromPosterior(posteriorNormalized, cytosineAndAdjacent, i, position);

        }

        tmpMethy = getMethyLevelFromPileup(pileup, ref);
        BisulfiteDiploidSNPGenotypeLikelihoods tmpGL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, ref, priors, BAC, tmpMethy);
        tmpGL.setPriors(tracker, ref, BAC.heterozygosity, BAC.novelDbsnpHet, BAC.validateDbsnpHet, location);
        //boolean firstSeen = true;
        boolean cytosinePatternNegStrand = false;
        boolean heterozygousPattern = false; // heterozygous at cytosine psoition
        byte basesAlelleA = 'C';
        byte basesAlelleB = 'C';
        //check cytosine position first, determine the cytosine pattern strand
        tmpGL.clearLikelihoods(tmpMethy);

        boolean notCytosinePos = false;
        int nGoodBases = tmpGL.add(pileup, true, true);
        if (nGoodBases == 0)
            return bestCytosinePattern;
        double[] posteriorNormalized = MathUtils.normalizeFromLog10(tmpGL.getPosteriors(), true, false);

        getBestGenotypeFromPosterior(posteriorNormalized, cytosineAndAdjacent, 0, position);


        methyStatus tmpMethyStatusCytosinePos = cytosineAndAdjacent.get(0);
        if (tmpMethyStatusCytosinePos == null) {
            notCytosinePos = true;
        } else if (tmpMethyStatusCytosinePos.genotype == null) {
            notCytosinePos = true;
        } else if (tmpMethyStatusCytosinePos.genotype.isHet()) {
            // for CpG, if C is heterozygous,then marked it here.
            if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(BaseUtils.Base.C.base, tmpMethyStatusCytosinePos.genotype.base1)
                    || BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(BaseUtils.Base.C.base, tmpMethyStatusCytosinePos.genotype.base2)) {

                if (tmpMethyStatusCytosinePos.ratio < this.BAC.STANDARD_CONFIDENCE_FOR_CALLING) {
                    notCytosinePos = true;
                } else {
                    heterozygousPattern = true;
                    cytosinePatternNegStrand = false;
                    basesAlelleA = tmpMethyStatusCytosinePos.genotype.base1;
                    basesAlelleB = tmpMethyStatusCytosinePos.genotype.base2;
                }

            } else if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(BaseUtils.Base.C.base, BaseUtilsMore.iupacCodeComplement(tmpMethyStatusCytosinePos.genotype.base1))
                    || BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(BaseUtils.Base.C.base, BaseUtilsMore.iupacCodeComplement(tmpMethyStatusCytosinePos.genotype.base2))) {

                if (tmpMethyStatusCytosinePos.ratio < this.BAC.STANDARD_CONFIDENCE_FOR_CALLING) {
                    notCytosinePos = true;
                } else {
                    heterozygousPattern = true;
                    cytosinePatternNegStrand = true;
                    basesAlelleA = tmpMethyStatusCytosinePos.genotype.base1;
                    basesAlelleB = tmpMethyStatusCytosinePos.genotype.base2;
                }


            } else {
                notCytosinePos = true;
            }

        } else {

            if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(BaseUtils.Base.C.base, tmpMethyStatusCytosinePos.genotype.base1)) {

                if (tmpMethyStatusCytosinePos.ratio < this.BAC.STANDARD_CONFIDENCE_FOR_CALLING) {
                    notCytosinePos = true;
                } else {
                    heterozygousPattern = false;
                    cytosinePatternNegStrand = false;
                }


            } else if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(BaseUtils.Base.C.base, BaseUtilsMore.iupacCodeComplement(tmpMethyStatusCytosinePos.genotype.base1))) {

                if (tmpMethyStatusCytosinePos.ratio < this.BAC.STANDARD_CONFIDENCE_FOR_CALLING) {
                    notCytosinePos = true;
                } else {
                    heterozygousPattern = false;
                    cytosinePatternNegStrand = true;
                }


            } else {
                notCytosinePos = true;
            }

        }
        //look at fwd or rvd for each cytosine type

        //pick up the longest pattern, rather than the max likelihood pattern, since some of the false positive may be due to uneven reads distribution...


        //for (String cytosineType : BAC.cytosineDefined.getContextDefined().keySet()) {
        for (Entry<String, CytosineParameters> cytosineTypeEntrySet : BAC.cytosineDefined.getContextDefined().entrySet()) {
            String cytosineType = cytosineTypeEntrySet.getKey();


            boolean heterozygousAtContext = false; // heterozygous at context

            int cytosinePos = cytosineTypeEntrySet.getValue().cytosinePosition;

            //check the reference cytosine pattern list
            CytosineParameters cps = new CytosineParameters();
            cps.isReferenceCytosinePattern = BisSNPUtils.isRefCytosinePattern(ref, cytosineType, cytosinePos);
            cytosineParametersStatus.put(cytosineType, cps);
            if (notCytosinePos)
                continue;

            int i = 1;
            int countMatched = 1;
            byte[] basesAlelleARev = cytosineType.getBytes();
            byte[] basesAlelleBRev = cytosineType.getBytes();
            byte[] basesAlelleAFwd = cytosineType.getBytes();
            byte[] basesAlelleBFwd = cytosineType.getBytes();

            if (cytosinePatternNegStrand) {
                // reverse strand
                if (heterozygousPattern) {
                    basesAlelleARev[cytosinePos - 1] = basesAlelleA;
                    basesAlelleBRev[cytosinePos - 1] = basesAlelleB;
                }
                for (byte base : cytosineType.getBytes()) {
                    int pos = cytosinePos - i;
                    int index = i - 1;
                    i++;

                    if (pos == 0) //cytosine position
                        continue;
                    methyStatus tmpMethyStatus = cytosineAndAdjacent.get(pos);
                    if (tmpMethyStatus == null) {
                        break;
                    } else if (tmpMethyStatus.genotype == null) {
                        break;
                    } else {
                        if (tmpMethyStatus.ratio < this.BAC.STANDARD_CONFIDENCE_FOR_CALLING) {
                            break;
                        }

                        if (tmpMethyStatus.genotype.isHet()) {
                            if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base1))
                                    || BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base2))) {
                                heterozygousAtContext = true;
                                basesAlelleARev[index] = tmpMethyStatus.genotype.base1;
                                basesAlelleBRev[index] = tmpMethyStatus.genotype.base2;
                                countMatched++;

                            } else {
                                break;
                            }
                        } else {

                            if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base1))) {
                                countMatched++;
                            } else {
                                break;
                            }

                        }
                    }

                }
            } else {
                // forward strand
                if (heterozygousPattern) {
                    basesAlelleAFwd[cytosinePos - 1] = basesAlelleA;
                    basesAlelleBFwd[cytosinePos - 1] = basesAlelleB;
                }
                //	if(cps.isReferenceCytosinePattern)
                //	if(ref.getLocus().getStart() == 7021736)
                //	System.err.println(ref.getLocus()+ "\t" + cytosineType + "\t" + BisSNPUtils.isRefCytosinePattern(ref, cytosineType, false) + "\t" + BisSNPUtils.isRefCytosinePattern(ref, cytosineType, true));

                for (byte base : cytosineType.getBytes()) {
                    int pos = i - cytosinePos;
                    int index = i - 1;

                    i++;
                    if (pos == 0)
                        continue;
                    methyStatus tmpMethyStatus = cytosineAndAdjacent.get(pos);
                    if (tmpMethyStatus == null) {
                        break;
                    } else if (tmpMethyStatus.genotype == null) {
                        break;
                    } else {
                        if (tmpMethyStatus.ratio < this.BAC.STANDARD_CONFIDENCE_FOR_CALLING) {
                            break;
                        }

                        if (tmpMethyStatus.genotype.isHet()) {
                            if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, tmpMethyStatus.genotype.base1)
                                    || BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, tmpMethyStatus.genotype.base2)) {
                                // isHeterozygousInContextPosition
                                heterozygousAtContext = true;
                                basesAlelleAFwd[index] = tmpMethyStatus.genotype.base1;
                                basesAlelleBFwd[index] = tmpMethyStatus.genotype.base2;
                                countMatched++;

                            } else {
                                break;
                            }
                        } else {

                            if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, tmpMethyStatus.genotype.base1)) {
                                countMatched++;
                            } else {
                                break;
                            }

                        }
                    }

                    //if(position == 10784184 && cytosineType.equalsIgnoreCase("HCH")){
                    //	System.err.println(i + "\t" + pos + "\t" + tmpMethyStatus.ratio + "\t" + tmpMethyStatus.genotype.isHet() + "\t" + base + "\t" + tmpMethyStatus.genotype.base1 + "\t" + tmpMethyStatus.genotype.base2 + "\t" + countMatched);
                    //}

                }
            }

            //if(position == 10784184 && cytosineType.equalsIgnoreCase("HCH")){
            //	System.err.println(countMatched + "\t" + bestCytosinePattern + "\t" + cytosinePatternNegStrand + "\t" + cytosinePos);
            //}

            if (countMatched < cytosineType.length() || heterozygousPattern) {
                continue;
            } else {
                if (bestCytosinePattern != null) {
                    if (bestCytosinePattern.length() >= countMatched) {
                        continue;
                    }
                }
                bestCytosinePattern = cytosineType;
                cps.isCytosinePattern = true;
                cps.cytosineMethylation = tmpMethy;
                cps.isHeterozygousCytosinePattern = heterozygousPattern;
                cps.isHeterozygousInContextPosition = heterozygousAtContext;

                if (cytosinePatternNegStrand) {
                    cps.cytosineStrand = '-';
                    if (heterozygousPattern || heterozygousAtContext) {
                        cps.patternOfAlleleA = new String(basesAlelleARev);
                        cps.patternOfAlleleB = new String(basesAlelleBRev);

                        if ((basesAlelleARev[cytosinePos - 1] == BaseUtils.Base.C.base && basesAlelleBRev[cytosinePos - 1] == BaseUtils.Base.T.base)
                                || (basesAlelleARev[cytosinePos - 1] == BaseUtils.Base.T.base && basesAlelleBRev[cytosinePos - 1] == BaseUtils.Base.C.base)
                                || (basesAlelleARev[cytosinePos - 1] == BaseUtils.Base.G.base && basesAlelleBRev[cytosinePos - 1] == BaseUtils.Base.A.base)
                                || (basesAlelleARev[cytosinePos - 1] == BaseUtils.Base.A.base && basesAlelleBRev[cytosinePos - 1] == BaseUtils.Base.G.base)) {
                            cps.isCTHeterozygousLoci = true;
                        }

                    }
                } else {
                    cps.cytosineStrand = '+';
                    if (heterozygousPattern || heterozygousAtContext) {
                        cps.patternOfAlleleA = new String(basesAlelleAFwd);
                        cps.patternOfAlleleB = new String(basesAlelleBFwd);

                        if ((basesAlelleAFwd[cytosinePos - 1] == BaseUtils.Base.C.base && basesAlelleBFwd[cytosinePos - 1] == BaseUtils.Base.T.base)
                                || (basesAlelleAFwd[cytosinePos - 1] == BaseUtils.Base.T.base && basesAlelleBFwd[cytosinePos - 1] == BaseUtils.Base.C.base)
                                || (basesAlelleAFwd[cytosinePos - 1] == BaseUtils.Base.G.base && basesAlelleBFwd[cytosinePos - 1] == BaseUtils.Base.A.base)
                                || (basesAlelleAFwd[cytosinePos - 1] == BaseUtils.Base.A.base && basesAlelleBFwd[cytosinePos - 1] == BaseUtils.Base.G.base)) {
                            cps.isCTHeterozygousLoci = true;
                        }

                    }
                }
                cytosineParametersStatus.put(cytosineType, cps);
            }


            if (position == BAC.testLocus) {
                System.err.println("countMatchedOnFwd: " + countMatched);
            }
        }

        return bestCytosinePattern;

    }


    public HashMap<String, BisulfiteContextsGenotypeLikelihoods> getBsContextGenotypeLikelihoods() {

        return BCGLs;
    }

    public Allele getRefAllele() {
        return refAllele;
    }

    public void setBsLikelihoods(RefMetaDataTracker tracker, ReferenceContext ref, HashMap<String, AlignmentContext> contexts, AlignmentContextUtils.ReadOrientation contextType,
                                 HashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs) {

        byte refBase = ref.getBase();

        if (!Allele.acceptableAlleleBases(new byte[]{refBase}))
            refBase = 'N';

        this.refAllele = Allele.create(refBase, true);

        for (Entry<String, AlignmentContext> sample : contexts.entrySet()) {
            // when pairs of reads overlapped, get rid of both of them when disagree, get rid of bad base qual reads when they agree
            ReadBackedPileup pileup = AlignmentContextUtils.stratify(sample.getValue(), contextType).getBasePileup().getOverlappingFragmentFilteredPileup();
            if (useBAQ) {
                //pileup = createBAQedPileup(pileup);
            }

            int numCNegStrand = 0;
            int numTNegStrand = 0;
            int numANegStrand = 0;
            int numGNegStrand = 0;
            int numCPosStrand = 0;
            int numTPosStrand = 0;
            int numAPosStrand = 0;
            int numGPosStrand = 0;

            int numGGalleleStrand = 0;
            int numAGalleleStrand = 0;
            int numOtherGalleleStrand = 0;
            int numCCalleleStrand = 0;
            int numTCalleleStrand = 0;
            int numOtherCalleleStrand = 0;

            for (PileupElement p : pileup) {
                SAMRecord samRecord = p.getRead();
                boolean negStrand = samRecord.getReadNegativeStrandFlag();

                int offset = p.getOffset();
                if (offset < 0)// is deletion
                    continue;

                boolean paired = samRecord.getReadPairedFlag();
                //synchronized(p.getRead()) {
                //	boolean goodBase = BisSNPUtils.goodPileupElement(p, BAC, ref);
                //	if(!goodBase)
                //		continue;
                //}


                if (paired && samRecord.getSecondOfPairFlag() && !BAC.nonDirectional) {

                    negStrand = !negStrand;
                }


                // summary number of C,T in the positive and negative strand

                if (BisSNPUtils.goodBaseInPileupElement(p, BAC, ref)) {
                    if (negStrand) {
                        if (p.getBase() == BaseUtils.Base.G.base) {
                            numCNegStrand++;
                        } else if (p.getBase() == BaseUtils.Base.A.base) {
                            numTNegStrand++;
                        } else if (p.getBase() == BaseUtils.Base.C.base) {
                            numGNegStrand++;
                        } else if (p.getBase() == BaseUtils.Base.T.base) {
                            numANegStrand++;
                        }

                    } else {
                        if (p.getBase() == BaseUtils.Base.C.base) {
                            numCPosStrand++;
                        } else if (p.getBase() == BaseUtils.Base.T.base) {
                            numTPosStrand++;
                        } else if (p.getBase() == BaseUtils.Base.G.base) {
                            numGPosStrand++;
                        } else if (p.getBase() == BaseUtils.Base.A.base) {
                            numAPosStrand++;
                        }
                    }

                    //System.err.println(CigarUtil.cigarArrayFromString(p.getRead().getCigarString())[offset]);
                    //if(CigarUtil.cigarArrayFromString(p.getRead().getCigarString())[offset]=='S'){
                    //	System.err.println(p.getRead().getCigarString() + "\t" + p.getRead().getReadName());
                    //}

                }


                if ((pileup.getLocation().getStart()) == testLoc) {
                    //System.err.println("before filter:\t" + onRefCoord + "\t" + offset + "\t" + negStrand + "\t" + pileup.getLocation().getStart() + "\t" + (char) p.getBase() + "\t" + p.getQual());
                    //System.err.println("refBase: " + refBase + "\tGoodBase: " + goodBase);

                    //if (paired)
                    //	System.err.println("isGoodBase: " + goodBase + "\tsecondOfPair: " + "\t" + samRecord.getSecondOfPairFlag());

                }
            }

            BisulfiteDiploidSNPGenotypePriors priors = new BisulfiteDiploidSNPGenotypePriors();

            HashMap<Integer, double[]> GPsBeforeCytosineTenGenotypes = new HashMap<Integer, double[]>();
            HashMap<Integer, double[]> GPsAfterCytosineTenGenotypes = new HashMap<Integer, double[]>();
            HashMap<String, Double[]> GPsAtCytosineTenGenotypes = new HashMap<String, Double[]>();
            HashMap<String, CytosineParameters> cytosineParametersStatus = new HashMap<String, CytosineParameters>();

            String bestMatchedCytosinePattern = checkCytosineStatus(pileup, BAC, tracker, ref, priors, GPsBeforeCytosineTenGenotypes, GPsAfterCytosineTenGenotypes, GPsAtCytosineTenGenotypes,
                    cytosineParametersStatus);


            char cytosineStrand;
            BisulfiteDiploidSNPGenotypeLikelihoods GL;

            GL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, ref, priors, BAC, getMethyLevelFromPileup(pileup, ref));
            if (cytosineParametersStatus.containsKey(bestMatchedCytosinePattern)) {
                cytosineStrand = cytosineParametersStatus.get(bestMatchedCytosinePattern).cytosineStrand;
            } else {
                cytosineStrand = '+';
                if ((pileup.getLocation().getStart()) == testLoc && bestMatchedCytosinePattern != null)
                    System.err.println(BAC.cytosineDefined.getContextDefined().get(bestMatchedCytosinePattern).cytosineMethylation + "\t" + bestMatchedCytosinePattern);
                bestMatchedCytosinePattern = null;
            }


            if ((pileup.getLocation().getStart()) == testLoc)
                GL.VERBOSE = true;

            GL.setPriors(tracker, ref, BAC.heterozygosity, BAC.novelDbsnpHet, BAC.validateDbsnpHet, ref.getLocus());

            int nGoodBases = GL.add(pileup, true, true);

            if (nGoodBases == 0)
                continue;
            if (cytosineStrand == '+') {
                numGGalleleStrand = numGNegStrand;
                numAGalleleStrand = numANegStrand;
                //if(BAC.nonDirectional){
                //	numOtherGalleleStrand = 0;
                //	numCCalleleStrand = numCPosStrand;
                //	numTCalleleStrand = numTPosStrand;
                //}
                //	else{
                numOtherGalleleStrand = numCNegStrand + numTNegStrand;
                numCCalleleStrand = numCPosStrand;
                numTCalleleStrand = numTPosStrand;
                //}

                numOtherCalleleStrand = numGPosStrand + numAPosStrand;
            } else {
                numGGalleleStrand = numGPosStrand;
                numAGalleleStrand = numAPosStrand;
                //	if(BAC.nonDirectional){
                //		numOtherGalleleStrand = 0;
                //		numCCalleleStrand = numCNegStrand;
                //		numTCalleleStrand = numTNegStrand;
                //	}
                //	else{
                numOtherGalleleStrand = numCPosStrand + numTPosStrand;
                numCCalleleStrand = numCNegStrand;
                numTCalleleStrand = numTNegStrand;
                //	}

                numOtherCalleleStrand = numGNegStrand + numANegStrand;
            }

            double[] prio = GL.getPriors();
            double[] likelihoods = GL.getLikelihoods();
            double[] posterior = MathUtils.normalizeFromLog10(GL.getPosteriors(), true, false);

            initializeBestAndAlternateAlleleFromPosterior(posterior, pileup.getLocation().getStart());

            if ((alternateAllele == null && bestAllele == refBase) || (bestAllele == null)) {

                if (BAC.OutputMode == BisulfiteEnums.OUTPUT_MODE.EMIT_VARIANTS_ONLY)

                    return;
            }

            Allele AlleleA, AlleleB;

            if (alternateAllele == null || BaseUtils.basesAreEqual(alternateAllele, refBase) || alternateAllele == bestAllele) {
                AlleleA = Allele.create(refBase, true);

                if (BaseUtils.basesAreEqual(bestAllele, refBase)) {
                    for (byte base : BaseUtils.BASES) {
                        if (base != refBase) {
                            bestAllele = base;
                            break;
                        }
                    }
                }

                AlleleB = Allele.create(bestAllele, false);

                alternateAllele = bestAllele;
                bestAllele = refBase;

            } else if (BaseUtils.basesAreEqual(bestAllele, refBase)) {
                AlleleA = Allele.create(bestAllele, true);
                AlleleB = Allele.create(alternateAllele, false);

            } else {
                AlleleA = Allele.create(bestAllele, false);
                AlleleB = Allele.create(alternateAllele, false);
                if (AlleleA.equals(refAllele, true)) {
                    AlleleA = Allele.create(bestAllele, true);
                }

                if (AlleleB.equals(refAllele, true)) {
                    AlleleB = Allele.create(alternateAllele, true);

                }
            }

            DiploidGenotype AAGenotype = DiploidGenotype.createHomGenotype(bestAllele);
            DiploidGenotype ABGenotype = DiploidGenotype.createDiploidGenotype(bestAllele, alternateAllele);
            DiploidGenotype BBGenotype = DiploidGenotype.createHomGenotype(alternateAllele);

            if ((pileup.getLocation().getStart()) == testLoc) {
                System.err.println("sample: " + sample.getKey());
                System.err.println("sample location: " + pileup.getPileupString((char) refBase));
                System.err.println("sample: " + sample.getValue().getLocation().getStart());
                System.err.println("refBase: " + refBase + " bestAllele: " + bestAllele + " alternateAllele: " + alternateAllele);
                System.err.println("AAGenotype: " + AAGenotype.toString() + " ABGenotype: " + ABGenotype.toString() + " BBGenotype: " + BBGenotype.toString());
                System.err.println("AlleleA: " + AlleleA.toString() + " AlleleB: " + AlleleB.toString());
                System.err.println("AAGenotype " + likelihoods[AAGenotype.ordinal()] + "\t" + prio[AAGenotype.ordinal()] + "\t" + posterior[AAGenotype.ordinal()]);
                System.err.println("ABGenotype " + likelihoods[ABGenotype.ordinal()] + "\t" + prio[ABGenotype.ordinal()] + "\t" + posterior[ABGenotype.ordinal()]);
                System.err.println("BBGenotype " + likelihoods[BBGenotype.ordinal()] + "\t" + prio[BBGenotype.ordinal()] + "\t" + posterior[BBGenotype.ordinal()]);
                System.err.println("Cytosine status: C-neg: " + numCNegStrand + "\tC-pos: " + numCPosStrand + "\tT-neg: " + numTNegStrand + "\tT-pos: " + numTPosStrand);

            }
            HashSet<String> cytosineContexts = new HashSet<String>(10);
            cytosineContexts.addAll(BAC.cytosineDefined.getContextDefined().keySet());

            BCGLs.put(sample.getKey(),
                    new BisulfiteContextsGenotypeLikelihoods(sample.getKey(), AlleleA, AlleleB, posterior[AAGenotype.ordinal()], posterior[ABGenotype.ordinal()], posterior[BBGenotype.ordinal()],
                            cytosineContexts, numCCalleleStrand, numTCalleleStrand, numOtherCalleleStrand, numGGalleleStrand, numAGalleleStrand, numOtherGalleleStrand, getFilteredDepth(pileup),
                            cytosineParametersStatus, bestMatchedCytosinePattern, GPsBeforeCytosineTenGenotypes, GPsAfterCytosineTenGenotypes, GPsAtCytosineTenGenotypes, pileup, BAC));
        }
        this.BCGLs = BCGLs;

    }

    private void getBestGenotypeFromPosterior(double[] posterior, HashMap<Integer, methyStatus> cytosineAdjacent, int key, int location) {
        double maxCount = Double.NEGATIVE_INFINITY;
        double secondMaxCount = Double.NEGATIVE_INFINITY;
        methyStatus tmpMethyStatus = new methyStatus();
        tmpMethyStatus.genotype = null;
        tmpMethyStatus.ratio = 0.0;
        DiploidGenotype bestGenotype = DiploidGenotype.createHomGenotype(BaseUtils.Base.A.base);

        for (DiploidGenotype g : DiploidGenotype.values()) {
            if (posterior[g.ordinal()] > maxCount) {
                secondMaxCount = maxCount;
                maxCount = posterior[g.ordinal()];

                bestGenotype = g;
            } else if (posterior[g.ordinal()] > secondMaxCount && posterior[g.ordinal()] <= maxCount) {
                secondMaxCount = posterior[g.ordinal()];
            }
        }
        tmpMethyStatus.ratio = 10 * (maxCount - secondMaxCount);
        if (location == BAC.testLocus) {
            System.err.println("maxCount: " + maxCount + "\tsecondMaxCount: " + secondMaxCount + "\tratio: " + tmpMethyStatus.ratio + "\tgenotype: " + bestGenotype);
            for (double poster : posterior) {
                System.err.println(poster);
            }
        }
        tmpMethyStatus.genotype = bestGenotype;
        cytosineAdjacent.put(key, tmpMethyStatus);

    }

    private int getFilteredDepth(ReadBackedPileup pileup) {
        int count = 0;
        for (PileupElement p : pileup) {
            if (BaseUtils.isRegularBase(p.getBase()))
                //count += p.getRepresentativeCount();
                count++;
        }

        return count;
    }

    private double getMethyLevelFromPileup(ReadBackedPileup pileup, ReferenceContext ref) {
        int numCNegStrand = 0;
        int numTNegStrand = 0;
        int numCPosStrand = 0;
        int numTPosStrand = 0;

        // for non directional protocol..
        int numGNegStrand = 0;
        int numANegStrand = 0;
        int numGPosStrand = 0;
        int numAPosStrand = 0;
        for (PileupElement p : pileup) {
            SAMRecord samRecord = p.getRead();
            boolean negStrand = samRecord.getReadNegativeStrandFlag();


            int offset = p.getOffset();
            if (offset < 0)// is deletion
                continue;
            boolean paired = samRecord.getReadPairedFlag();
            //	synchronized(p.getRead()) {
            //		boolean goodBase = BisSNPUtils.goodPileupElement(p, BAC, ref);
            //		if(!goodBase)
            //			continue;
            //	}


            if (paired && samRecord.getSecondOfPairFlag() && !BAC.nonDirectional) {

                negStrand = !negStrand;
            }
            if (BisSNPUtils.goodBaseInPileupElement(p, BAC, ref)) {
                if (negStrand) {
                    if (p.getBase() == BaseUtils.Base.G.base) {
                        numCNegStrand++;
                    } else if (p.getBase() == BaseUtils.Base.A.base) {
                        numTNegStrand++;
                    } else if (p.getBase() == BaseUtils.Base.C.base) {
                        numGNegStrand++;
                    } else if (p.getBase() == BaseUtils.Base.T.base) {
                        numANegStrand++;
                    }

                } else {
                    if (p.getBase() == BaseUtils.Base.C.base) {
                        numCPosStrand++;
                    } else if (p.getBase() == BaseUtils.Base.T.base) {
                        numTPosStrand++;
                    } else if (p.getBase() == BaseUtils.Base.G.base) {
                        numGPosStrand++;
                    } else if (p.getBase() == BaseUtils.Base.A.base) {
                        numAPosStrand++;
                    }

                }
            }


        }
        int methy = numCNegStrand + numCPosStrand;
        int unmethyPlusMethy = numCNegStrand + numTNegStrand + numCPosStrand + numTPosStrand;
        //if(BAC.nonDirectional){
        //	methy = numCNegStrand + numCPosStrand + numGNegStrand + numGPosStrand;
        //	unmethyPlusMethy = numCNegStrand + numTNegStrand + numCPosStrand + numTPosStrand + numGNegStrand + numANegStrand + numGPosStrand + numAPosStrand;
        //}
        double methyLevel = (double) methy / (double) (unmethyPlusMethy);

        if (unmethyPlusMethy == 0)
            return 0;
        else
            return methyLevel;

    }

    private void initializeBestAndAlternateAlleleFromPosterior(double[] posterior, int location) {
        double maxCount = Double.NEGATIVE_INFINITY;
        double secondMaxCount = Double.NEGATIVE_INFINITY;
        DiploidGenotype bestGenotype = DiploidGenotype.createHomGenotype(BaseUtils.Base.A.base);
        DiploidGenotype secondGenotype = DiploidGenotype.createHomGenotype(BaseUtils.Base.A.base);
        bestAllele = null;
        alternateAllele = null;

        for (DiploidGenotype g : DiploidGenotype.values()) {
            if (posterior[g.ordinal()] > maxCount) {
                secondMaxCount = maxCount;
                maxCount = posterior[g.ordinal()];
                if (bestGenotype.base1 != secondGenotype.base1) {
                    secondGenotype = bestGenotype;
                }
                bestGenotype = g;
            } else if (posterior[g.ordinal()] > secondMaxCount && posterior[g.ordinal()] <= maxCount) {
                secondMaxCount = posterior[g.ordinal()];
                secondGenotype = g;
            }
        }
        if (bestGenotype.isHom()) {
            bestAllele = bestGenotype.base1;
            if (secondGenotype.isHom()) {
                alternateAllele = secondGenotype.base1;
            } else {
                if (secondGenotype.base1 == bestAllele) {
                    alternateAllele = secondGenotype.base2;
                } else {
                    alternateAllele = secondGenotype.base1;
                }
            }

        } else {
            DiploidGenotype temp1 = DiploidGenotype.createHomGenotype(bestGenotype.base1);
            DiploidGenotype temp2 = DiploidGenotype.createHomGenotype(bestGenotype.base2);
            if (posterior[temp1.ordinal()] > posterior[temp2.ordinal()]) {
                bestAllele = bestGenotype.base1;
                alternateAllele = bestGenotype.base2;
            } else {
                bestAllele = bestGenotype.base2;
                alternateAllele = bestGenotype.base1;
            }
        }

        if (location == testLoc) {
            for (DiploidGenotype g : DiploidGenotype.values()) {
                System.err.println(g.base1 + "-" + g.base2 + ": " + posterior[g.ordinal()]);
            }
            System.err.println("bestAllele: " + bestAllele + "\t" + maxCount);
            if (alternateAllele != null) {
                System.err.println("AlternateAllele: " + "\t" + alternateAllele + "\t" + secondMaxCount);
            }
        }
    }

    // inner class to record genotype and posterior ratio of best and second
    // best genotype
    private class methyStatus {
        DiploidGenotype genotype;
        double ratio;

        methyStatus() {

        }
    }

}
