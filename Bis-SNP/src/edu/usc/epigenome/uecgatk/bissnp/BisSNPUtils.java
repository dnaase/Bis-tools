package edu.usc.epigenome.uecgatk.bissnp;

import gnu.trove.map.hash.THashMap;

import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;


import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

import org.broadinstitute.sting.utils.genotyper.DiploidGenotype;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;




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

/*
 * provide easy access for some function inside BisSNP. e.g. judge if it is a type of cytosine
 * pattern(HCG-2, or GCH-1) from the given pileup and corresponding reference seq, dbSNP. should
 * provide methylation value, otherwise will use flat methylation value; and also likelihood ratio
 * criteria, bisulfite conversion rate
 */

public class BisSNPUtils {
	
	private final static Pattern mdPat = Pattern.compile("\\G(?:([0-9]+)|([ACTGNactgn])|(\\^[ACTGNactgn]+))");
	
	
	
	public static boolean isCytosine(byte base, boolean bisulfiteConversionSpace) {
		char refC = (char) base;
		boolean out;

		if (bisulfiteConversionSpace) {
			out = ((refC == 'C') || (refC == 'T'));
		} else {
			out = (refC == 'C');
		}

		return out;
	}

	public static boolean isCytosine(byte base, boolean bisulfiteConversionSpace, boolean secondPair) {
		char refC = (char) base;
		boolean out;

		if (bisulfiteConversionSpace) {
			if (secondPair) {
				out = ((refC == 'G') || (refC == 'A'));
			} else {
				out = ((refC == 'C') || (refC == 'T'));
			}

		} else {
			if (secondPair) {
				out = (refC == 'G');
			} else {
				out = (refC == 'C');
			}

		}

		return out;
	}

	public static boolean isCytosine(int pos, String seqStr, boolean bisulfiteConversionSpace) {
		char refC = seqStr.charAt(pos);

		boolean out;

		if (bisulfiteConversionSpace) {
			out = ((refC == 'C') || (refC == 'T'));
		} else {
			out = (refC == 'C');
		}

		return out;
	}

	public static boolean isCytosine(int pos, String seqStr, boolean bisulfiteConversionSpace, boolean secondPair) {
		char refC = seqStr.charAt(pos);

		boolean out;

		if (bisulfiteConversionSpace) {
			if (secondPair) {
				out = ((refC == 'G') || (refC == 'A'));
			} else {
				out = ((refC == 'C') || (refC == 'T'));
			}

		} else {
			if (secondPair) {
				out = (refC == 'G');
			} else {
				out = (refC == 'C');
			}

		}

		return out;
	}

	// het at C sites or G sites.
	public static boolean isHetCpattern(VariantContext vc) {
		if (vc.hasGenotypes()) {
			Iterator<Genotype> it = vc.getGenotypes().iterator();
			while (it.hasNext()) {
				Genotype tmp = it.next();
				if (tmp.getGenotypeString().equalsIgnoreCase("C/C") || tmp.getGenotypeString().equalsIgnoreCase("G/G")) {
					return true;
				}
			}
		}
		return false;
	}

	public static boolean isHetCpg_at_C(VariantContext vc) {
		if (vc.hasGenotypes()) {
			Iterator<Genotype> it = vc.getGenotypes().iterator();
			while (it.hasNext()) {
				Genotype tmp = it.next();
				if (tmp.isHom())
					return false;
				if (tmp.hasExtendedAttribute(BisulfiteVCFConstants.BEST_C_PATTERN)) {
					byte[] bases = tmp.getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".").getBytes();
					if (bases.length < 2)
						return false;
					if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[0], BaseUtilsMore.C) && BaseUtils.basesAreEqual(bases[1], BaseUtilsMore.G)) {
						return true;
					}

				}
			}
		}
		return false;
	}

	public static boolean isHetCpg_at_G(VariantContext vc) {
		if (vc.hasGenotypes()) {
			Iterator<Genotype> it = vc.getGenotypes().iterator();
			while (it.hasNext()) {
				Genotype tmp = it.next();
				if (tmp.isHom()) {
					if (tmp.hasExtendedAttribute(BisulfiteVCFConstants.BEST_C_PATTERN)) {
						byte[] bases = tmp.getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".").getBytes();
						if (bases.length < 2)
							return false;
						if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[1], BaseUtilsMore.G) && BaseUtils.basesAreEqual(bases[0], BaseUtilsMore.C)) {
							return true;
						}

					}
				}
			}
		}
		return false;
	}

	public static boolean isHetCph_at_C(VariantContext vc) {
		if (vc.hasGenotypes()) {
			Iterator<Genotype> it = vc.getGenotypes().iterator();
			while (it.hasNext()) {
				Genotype tmp = it.next();
				if (tmp.isHom())
					return false;
				if (tmp.hasExtendedAttribute(BisulfiteVCFConstants.BEST_C_PATTERN)) {
					byte[] bases = tmp.getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".").getBytes();
					if (bases.length < 2)
						return false;
					if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[0], BaseUtilsMore.C) && bases[1] == BaseUtilsMore.H) {
						return true;
					}

				}
			}
		}
		return false;
	}

	public static boolean isHomoC(VariantContext vc) {
		if (vc.hasGenotypes()) {
			Iterator<Genotype> it = vc.getGenotypes().iterator();
			while (it.hasNext()) {
				Genotype tmp = it.next();
				// System.err.println(tmp.getGenotypeString());
				if (tmp.getGenotypeString().equalsIgnoreCase("C/C") || tmp.getGenotypeString().equalsIgnoreCase("G/G")) {
					return true;
				}
			}
		}
		return false;
	}

	public static boolean isRefCpg(ReferenceContext refContext) {
		byte[] refBytes = new byte[3];
		for (int i = -1, index = 0; i <= 1; i++, index++) {
			GenomeLoc loc = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), refContext.getLocus().getStart() + i);
			if (!refContext.getWindow().containsP(loc)) {
				refBytes[index] = BaseUtilsMore.N;
				continue;
			}

			ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(), loc, refContext.getWindow(), refContext.getBases());
			refBytes[index] = tmpRef.getBase();
		}
		if (BaseUtils.basesAreEqual(refBytes[0], BaseUtilsMore.C) && BaseUtils.basesAreEqual(refBytes[1], BaseUtilsMore.G)) {
			return true;
		} else if (BaseUtils.basesAreEqual(refBytes[1], BaseUtilsMore.C) && BaseUtils.basesAreEqual(refBytes[2], BaseUtilsMore.G)) {
			return true;
		} else {
			return false;
		}
	}

	public static boolean isRefCph(ReferenceContext refContext) {
		byte[] refBytes = new byte[3];
		for (int i = -1, index = 0; i <= 1; i++, index++) {
			GenomeLoc loc = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), refContext.getLocus().getStart() + i);
			if (!refContext.getWindow().containsP(loc)) {
				refBytes[index] = BaseUtilsMore.N;
				continue;
			}

			ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(), loc, refContext.getWindow(), refContext.getBases());
			refBytes[index] = tmpRef.getBase();
		}
		if (BaseUtils.basesAreEqual(refBytes[1], BaseUtilsMore.G)
				&& (!BaseUtils.basesAreEqual(refBytes[0], BaseUtilsMore.C) && !BaseUtils.basesAreEqual(refBytes[0], BaseUtilsMore.N) && !BaseUtils.basesAreEqual(refBytes[0],  BaseUtilsMore.N))) {
			return true;
		} else if (BaseUtils.basesAreEqual(refBytes[1], BaseUtilsMore.C)
				&& (!BaseUtils.basesAreEqual(refBytes[2], BaseUtilsMore.G) && !BaseUtils.basesAreEqual(refBytes[2], BaseUtilsMore.N) && !BaseUtils.basesAreEqual(refBytes[2],  BaseUtilsMore.N))) {
			;
			return true;
		} else {

			return false;
		}
	}

	public static boolean isRefCytosinePattern(ReferenceContext refContext, String cytosinePattern) {
		if (isRefCytosinePattern(refContext, cytosinePattern, false) || isRefCytosinePattern(refContext, cytosinePattern, true))
			return true;
		return false;
	}

	public static boolean isRefCytosinePattern(ReferenceContext refContext, String cytosinePattern, boolean negStrand) {

		byte[] bases = cytosinePattern.getBytes();
		int matches = 0;
		if (negStrand) {
			byte[] refBytesRev = new byte[cytosinePattern.length()];

			for (int i = 0; i < cytosinePattern.length(); i++) {
				GenomeLoc loc = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), refContext.getLocus().getStart() - i);
				if (!refContext.getWindow().containsP(loc)) {
					refBytesRev[i] = BaseUtilsMore.N;
					continue;
				}

				ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(), loc, refContext.getWindow(), refContext.getBases());
				refBytesRev[i] = tmpRef.getBase();
			}
			//if(refContext.getLocus().getStart() == 7021736)
			//	System.err.println(new String(refBytesRev) + "\t" + cytosinePattern);
			
			refBytesRev = BaseUtils.simpleComplement(refBytesRev);
			
			//if(refContext.getLocus().getStart() == 7021736)
			//	System.err.println(new String(refBytesRev) + "\t" + cytosinePattern);
			
			for (int i = 0; i < bases.length; i++) {
				if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[i], refBytesRev[i])) {
					matches++;
				}
			}
			if (matches == bases.length)
				return true;
		} else {
			byte[] refBytesFwd = new byte[cytosinePattern.length()];
			for (int i = 0; i < cytosinePattern.length(); i++) {
				GenomeLoc loc = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), refContext.getLocus().getStart() + i);
				if (!refContext.getWindow().containsP(loc)) {
					refBytesFwd[i] = BaseUtils.NO_CALL_INDEX;
					continue;
				}

				ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(), loc, refContext.getWindow(), refContext.getBases());
				refBytesFwd[i] = tmpRef.getBase();
			}

			for (int i = 0; i < bases.length; i++) {
				if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[i], refBytesFwd[i])) {
					matches++;
				}
			}
			if (matches == bases.length)
				return true;
		}

		return false;
	}
	
	public static boolean isRefCytosinePattern(ReferenceContext refContext, String cytosinePattern, int cytosinePos) {
		if (isRefCytosinePattern(refContext, cytosinePattern, cytosinePos, false) || isRefCytosinePattern(refContext, cytosinePattern, cytosinePos, true))
			return true;
		return false;
	}

	public static boolean isRefCytosinePattern(ReferenceContext refContext, String cytosinePattern, int cytosinePos, boolean negStrand) {

		byte[] bases = cytosinePattern.getBytes();
		int matches = 0;
		if (negStrand) {
			byte[] refBytesRev = new byte[cytosinePattern.length()];

			for (int i = 0; i < cytosinePattern.length(); i++) {
				GenomeLoc loc = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), refContext.getLocus().getStart() - (i-cytosinePos+1));
				if (!refContext.getWindow().containsP(loc)) {
					refBytesRev[i] = BaseUtils.NO_CALL_INDEX;
					continue;
				}

				ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(), loc, refContext.getWindow(), refContext.getBases());
				refBytesRev[i] = tmpRef.getBase();
			}
			
			refBytesRev = BaseUtilsMore.simpleIupacCodeComplement(refBytesRev);
			
			for (int i = 0; i < bases.length; i++) {
				if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[i], refBytesRev[i])) {
					matches++;
				}
			}
			if (matches == bases.length)
				return true;
		} else {
			byte[] refBytesFwd = new byte[cytosinePattern.length()];
			for (int i = 0; i < cytosinePattern.length(); i++) {
				GenomeLoc loc = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), refContext.getLocus().getStart() + (i-cytosinePos+1));
				if (!refContext.getWindow().containsP(loc)) {
					refBytesFwd[i] = BaseUtils.NO_CALL_INDEX;
					continue;
				}

				ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(), loc, refContext.getWindow(), refContext.getBases());
				refBytesFwd[i] = tmpRef.getBase();
			}
			
			for (int i = 0; i < bases.length; i++) {
				if (BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[i], refBytesFwd[i])) {
					matches++;
				}
			}
			if (matches == bases.length)
				return true;
		}

		return false;
	}

	public static boolean isThymine(int pos, String seqStr) {
		char seqC = seqStr.charAt(pos);

		return (seqC == 'T');
	}

	public static boolean isThymineInCytosinePos(int pos, String seqStr, byte refBase) {
		char seqC = seqStr.charAt(pos);

		return ((seqC == 'T') && (refBase == 'C'));
	}
	
	public static boolean checkPattern(String motif, ReferenceContext ref, boolean negStrand){
		byte[] refBytes = new byte[motif.length()];
		byte[] motifSeq = motif.getBytes();
		if(negStrand)
			motifSeq = BaseUtilsMore.simpleReverseIupacCodeComplement(motifSeq);
		int start = negStrand? -(motif.length()-1) : 0;
		int end = negStrand? 0 : motif.length()-1;
		

		for(int i = start, index = 0; i <= end; i++, index++){
			GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(), ref.getLocus().getStart()+i );
			if( !ref.getWindow().containsP(loc) )
				return false;
			
			ReferenceContext tmpRef = new ReferenceContext(ref.getGenomeLocParser(),loc, ref.getWindow(),ref.getBases());
			refBytes[index] = tmpRef.getBase();
			if( !BaseUtils.isRegularBase(refBytes[index]) || !BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(motifSeq[index], refBytes[index]))
				return false;
		}

		return true;
	}
	
	public static String now(String dateFormat) {
		Calendar cal = Calendar.getInstance();
		SimpleDateFormat sdf = new SimpleDateFormat(dateFormat);
		return sdf.format(cal.getTime());

	}

	public class methyStatus {
		DiploidGenotype genotype;
		double ratio;
	}
	
	
	/*
	public static boolean goodPileupElement(PileupElement p, BisulfiteArgumentCollection BAC, ReferenceContext refContext) {
		short cPos=1;
		return goodPileupElement(p, BAC, refContext, "C", cPos);
	}

	//
	public static boolean goodPileupElement(PileupElement p, BisulfiteArgumentCollection BAC, ReferenceContext refContext, String patConv5, short posCinPatConv5) {
		BitSet bitSetRead;
		int offSet = p.getOffset();
		GATKSAMRecord read = p.getRead();
		if (read.containsTemporaryAttribute(BisulfiteSAMConstants.BIT_WISE_TAG)) {
			synchronized(read){
				bitSetRead = (BitSet)read.getTemporaryAttribute(BisulfiteSAMConstants.BIT_WISE_TAG);
			}
			//System.err.println(bitSetRead + " BIT_WISE_TAG: " + p.getRead().getTemporaryAttribute(BIT_WISE_TAG));
		}else{
			
			bitSetRead = new BitSet(read.getReadLength());
			 if(badReads(read,refContext, BAC) || badReadsByMismatchesInWindow(read,refContext, BAC, p) || badReadsByBisulfiteConversion(read,refContext, BAC,p)){
				 //System.err.println("no good reads " + read.getReadName());
			 }else{
				 bitSetRead = notTrimmedBase(offSet, read, BAC, bitSetRead);
				 BitSet bitSetRead2 = goodBasePassBisulfiteConversionFivePrime(offSet, read, BAC,refContext, patConv5, posCinPatConv5);
				 bitSetRead.and(bitSetRead2);
				// System.err.println("BitWiseTag" + bitSetRead.toString() + "\t" + bitSetRead.cardinality());
			 }
			 synchronized(read){
				 read.setTemporaryAttribute(BisulfiteSAMConstants.BIT_WISE_TAG, bitSetRead);
			 }
			 
		}
		if(bitSetRead.get(offSet)){
			//System.err.println(goodBase(p, BAC, refContext, offSet));
			return goodBase(p, BAC, refContext, offSet);

		}else{
			return false;
		}
		
		
	}
	*/
	
	public static boolean goodBaseInPileupElement(PileupElement p, BisulfiteArgumentCollection BAC, ReferenceContext refContext) {
		byte[] quals = p.getRead().getBaseQualities();
		if(BAC.useOriginalQual){
			quals = p.getRead().getOriginalBaseQualities();
		}
		//System.err.println(!notTrimmedBase(p.getOffset(), p.getRead(), BAC) + "\t" + !goodConvBase(p.getOffset(), p.getRead(), BAC, refContext));
		if (p.isDeletion() || p.getQual() == 0 || BaseUtils.simpleBaseToBaseIndex(p.getBase()) == -1 || quals[p.getOffset()] < BAC.MIN_BASE_QUALTY_SCORE || !notTrimmedBase(p.getOffset(), p.getRead(), BAC) || !goodConvBase(p.getOffset(), p.getRead(), BAC, refContext))
			return false;
		return true;
	}
	
	//filter out bad base
	public static boolean goodBase(PileupElement p, BisulfiteArgumentCollection BAC, ReferenceContext refContext, int offSet) {
		byte[] quals = p.getRead().getBaseQualities();
		if(BAC.useOriginalQual){
			quals = p.getRead().getOriginalBaseQualities();
		}
		if (p.isDeletion() || p.getQual() == 0 || BaseUtils.simpleBaseToBaseIndex(p.getBase()) == -1 || quals[offSet] < BAC.MIN_BASE_QUALTY_SCORE )
			return false;
		return true;
	}
	
	//TO DO: check in 2nd end, what is the offset?
	/*
	public static BitSet notTrimmedBase(int offset, SAMRecord samRecord, BisulfiteArgumentCollection BAC, BitSet bitSetRead){

		if(BAC.trim2nd){
			if(samRecord.getReadPairedFlag() && samRecord.getSecondOfPairFlag()){
				if(samRecord.getReadNegativeStrandFlag()){
					bitSetRead.set(BAC.trim3, samRecord.getReadLength() - BAC.trim5);	
				}
				else{
					bitSetRead.set(BAC.trim5, samRecord.getReadLength() - BAC.trim3);	
				}
			}
			
		}
		else{
			if(samRecord.getReadNegativeStrandFlag()){
				bitSetRead.set(BAC.trim3, samRecord.getReadLength() - BAC.trim5);		
			}
			else{
				bitSetRead.set(BAC.trim5, samRecord.getReadLength() - BAC.trim3);	
					
			}
		}
		return bitSetRead;
		
	}
	*/
	public static boolean notTrimmedBase(int offset, SAMRecord samRecord, BisulfiteArgumentCollection BAC){

		if(BAC.trim2nd){
			if(samRecord.getReadPairedFlag() && samRecord.getSecondOfPairFlag()){
				if(samRecord.getReadNegativeStrandFlag()){
					return offset >= BAC.trim3 && offset <= samRecord.getReadLength() - BAC.trim5;
					
				}
				else{
					return offset >= BAC.trim5 && offset <= samRecord.getReadLength() - BAC.trim3;
						
				}
			}
			else{
				return true;
			}
		}
		else{
			if(samRecord.getReadNegativeStrandFlag()){
				return offset >= BAC.trim3 && offset <= samRecord.getReadLength() - BAC.trim5;
					
			}
			else{
				return offset >= BAC.trim5 && offset <= samRecord.getReadLength() - BAC.trim3;

			}
		}
		
	}
	
	public static boolean goodConvBase(int offset, GATKSAMRecord samRecord, BisulfiteArgumentCollection BAC, ReferenceContext ref){
		if(samRecord.containsTemporaryAttribute(BisulfiteSAMConstants.BIT_WISE_TAG)){
			int convStart = (Integer) samRecord.getTemporaryAttribute(BisulfiteSAMConstants.BIT_WISE_TAG);
			if(ref == null || convStart == -1)
				return true;
			byte refBase=ref.getBase();
			if(samRecord.getReadNegativeStrandFlag()){
				if(offset > samRecord.getReadLength()-convStart){ //unconverted region in the read, the cytosine strand are filter out for methylation calling
					if(samRecord.getReadPairedFlag() && samRecord.getSecondOfPairFlag()){
						if(BaseUtils.basesAreEqual(refBase, BaseUtils.C)){
							return false;
						}
					}else{
						if(BaseUtils.basesAreEqual(refBase, BaseUtils.G)){
							return false;
						}
					}
				}
			}
			else{
				if(offset < convStart){//unconverted region in the read, the cytosine strand are filter out for methylation calling
					if(samRecord.getReadPairedFlag() && samRecord.getSecondOfPairFlag()){
						if(BaseUtils.basesAreEqual(refBase, BaseUtils.G)){
							return false;
						}
					}else{
						if(BaseUtils.basesAreEqual(refBase, BaseUtils.C)){
							return false;
						}
					}
				}
			}
			
		}
			
			return true;
	}
	
	/*
	//filter out bad reads
	public static boolean badReads(SAMRecord samRecord,ReferenceContext ref,BisulfiteArgumentCollection BAC) {
		
		//check bad single end reads
		if ( samRecord.getMappingQuality() < BAC.MIN_MAPPING_QUALTY_SCORE || samRecord.getReadUnmappedFlag() || samRecord.getDuplicateReadFlag()
			|| samRecord.getReadFailsVendorQualityCheckFlag() || samRecord.getNotPrimaryAlignmentFlag()) {
			return true;
		}
		
		//check bad paired end reads
		boolean paired = samRecord.getReadPairedFlag();
		if (paired) {
			//check Inverted duplicated reads
			if (samRecord.getAlignmentStart() == samRecord.getMateAlignmentStart() && samRecord.getReadNegativeStrandFlag() == samRecord.getMateNegativeStrandFlag()) {
				if(BAC.invDups == INVERT_DUPS.USE_ONLY_1ST_END){
					if (samRecord.getSecondOfPairFlag()) 
						return true;
				}
				else if(BAC.invDups == INVERT_DUPS.NOT_TO_USE){
					return true;
				}

			}
			//check improper paired or bad mate reads
			if (paired && !BAC.USE_BADLY_MATED_READS && BadMateFilter.hasBadMate(samRecord) && !samRecord.getProperPairFlag()) {
				return true;
			}
			
		}
		return false;
		
	}
	*/
	
	/*
	//filter out reads have many nonbisulfite mismatches.
	public static boolean badReadsByMismatchesInWindow(SAMRecord read, ReferenceContext ref, BisulfiteArgumentCollection BAC, PileupElement p) {
		int readLength = read.getReadLength();
		boolean negativeStrand = read.getReadNegativeStrandFlag();
		boolean secondEnd = read.getReadPairedFlag() && read.getSecondOfPairFlag();

		byte[] bases = read.getReadBases();
		
		//TODO: need to check window length in boundary...
		byte[] refBases =getRefBaseInRefWindow(ref, readLength, p.getOffset());
		//System.err.println(new String(refBasesInWidnow) + "\t" + refBasesInWidnow.length);
		//System.err.println(new String(refBases) + "\t" + len);
		//System.err.println(new String(bases) + "\t" + readLength + "\t" + negativeStrand + "\t" + secondEnd);
		int numberOfMismatches = 0;
		for(int i = 0; i < readLength; i++){
			if( !BaseUtils.basesAreEqual(refBases[i], bases[i]) && BaseUtilsMore.isBisulfiteMismatch(refBases[i], bases[i],negativeStrand, secondEnd))
				numberOfMismatches++;
			if(numberOfMismatches > BAC.MAX_MISMATCHES)
				return true;
		}
		//System.err.println(numberOfMismatches);
		return false;
	}
	*/
	
	/*
	//mask five prime pattern before converted pattern as 0 in bitset
	public static BitSet goodBasePassBisulfiteConversionFivePrime(int offset, SAMRecord samRecord, BisulfiteArgumentCollection BAC,ReferenceContext ref, String patConv5, short posCinPatConv5) {	
		int readLength = samRecord.getReadLength();
		BitSet bitSetRead = new BitSet(readLength);
		bitSetRead.flip(0, readLength);
		byte[] bases = samRecord.getReadBases();
		byte[] baseQs = samRecord.getBaseQualities();
		if (ref == null)
			return bitSetRead;
		
		//TODO: need to check window length in boundary...
		byte[] refBases =getRefBaseInRefWindow(ref, readLength,offset);
		
		
		boolean negStrand = samRecord.getReadNegativeStrandFlag();
		if(samRecord.getReadPairedFlag() && samRecord.getSecondOfPairFlag())
			negStrand = !negStrand;
		
		if(negStrand){
			bases = BaseUtils.simpleReverseComplement(bases);
			baseQs = BaseUtilsMore.simpleReverse(baseQs);
			refBases = BaseUtils.simpleReverseComplement(refBases);	
		}
		byte[] patterns = patConv5.getBytes();
		
		if(readLength < patterns.length )
			return bitSetRead;
		
		short convertedCount = 0;
		for(int i = 0; i <= readLength-patterns.length; i++){
			short numMatchesInRef=0;
			short numMatchesInReads=0;
			boolean conv = false;
			for(int j = i, index = 0; index < patterns.length; index++, j++){
				if(baseQs[j] < BAC.MIN_BASE_QUALTY_SCORE){	
					break;
				}
					
				if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(patterns[index], refBases[j])){
					numMatchesInRef++;
						
				}
				if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(patterns[index], bases[j])){
					numMatchesInReads++;
				}
				if(index == (posCinPatConv5-1) && BaseUtils.basesAreEqual(bases[j], BaseUtils.T)){
					conv=true;
				}
				
			}
			if(numMatchesInReads == patterns.length){
				if( numMatchesInRef == (patterns.length-1)){
					if(conv)
						convertedCount++;
					
					if(convertedCount < BAC.minConv){
						bitSetRead.flip(i);
					}
				}
			}
			
		}
		
		return bitSetRead;
	}
	*/
	
	/*
	//filter out reads not completed converted (HCH by default)
	public static boolean badReadsByBisulfiteConversion(SAMRecord read, ReferenceContext ref,BisulfiteArgumentCollection BAC, PileupElement p) {
		int readLength = read.getReadLength();
		boolean negativeStrand = read.getReadNegativeStrandFlag();
		if(read.getReadPairedFlag() && read.getSecondOfPairFlag())
			negativeStrand = !negativeStrand;
		
		String pattern = BAC.patConv;
		byte[] bases = read.getReadBases();
		byte[] baseQs = read.getBaseQualities();
		short numberOfPatternInRef = 0; 

		
		//TODO: need to check window length in boundary...
		
		byte[] refBases =getRefBaseInRefWindow(ref, readLength, p.getOffset());
		

		if(negativeStrand){
			bases = BaseUtils.simpleReverseComplement(bases);
			baseQs = BaseUtilsMore.simpleReverse(baseQs);
			if(numberOfPatternInRef != -1){
				refBases = BaseUtils.simpleReverseComplement(refBases);
			}
			
		}
		
		if(readLength < pattern.length())
			return true;
		int numberOfPattern = 0;
		
		byte[] patterns = pattern.getBytes();
		
		
		for(int i = 0; i <= readLength-patterns.length; i++){
			short numMatchesInRef=0;
			short numMatchesInReads=0;
			for(int j = i, index = 0; index < patterns.length; index++, j++){
				if(baseQs[j] < BAC.MIN_BASE_QUALTY_SCORE){	
					break;
				}
				if(numberOfPatternInRef != -1){ //mean reads is mapped
					
					if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(patterns[index], refBases[j])){
						numMatchesInRef++;
						
					}
					
				}
				if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(patterns[index], bases[j])){
					numMatchesInReads++;
				}
				
			}
			if(numMatchesInReads == patterns.length){
				numberOfPattern++;
			}
			if(numberOfPatternInRef != -1 && numMatchesInRef == patterns.length){
				numberOfPatternInRef++;
			}
		}
		if(((numberOfPatternInRef == -1 || numberOfPatternInRef == 0) && numberOfPattern == 0) || (numberOfPatternInRef != -1 && (double)numberOfPattern/(double)numberOfPatternInRef < BAC.minPatConv)){
			return false;
		}else{
			//System.err.println(new String(refBasesInWidnow) + "\t" + refBasesInWidnow.length);
			//		System.err.println(new String(refBases) + "\t" + len);
			//		System.err.println(new String(bases) + "\t" + readLength + "\t" + negativeStrand + "\t" + read.getSecondOfPairFlag() + "\t" + numberOfPattern + "\t" + numberOfPatternInRef);
			//		System.err.println(new String(baseQs) + "\t" + read.getReadName());
			return true;
		}
		
	}
	*/
	
	
	public static byte[] getRefBaseInRefWindow(ReferenceContext ref, int readLen, int offsetInRead){
		byte[]  refBasesInWidnow = ref.getBases();	
		int a = BisulfiteSAMConstants.REF_WINDOW_LEN - offsetInRead;
		byte[] refBases = new byte[readLen];
		for (int i = 0; i < readLen; i++) {
			refBases[i]=refBasesInWidnow[a+i];
		}
		return refBases;
	}
	
	public static byte[] refStrFromMd(SAMRecord read) throws Exception{
		return refStrFromMd(read.getReadString(),read.getStringAttribute(BisulfiteSAMConstants.MD_TAG), read.getCigar(), read.getReadName());
	}
	
	
	//From ben's picardUtils, need to test it in paired end space
	public static byte[] refStrFromMd(String seq, String md, Cigar cigar, String readName)
			throws Exception
			{
				if (seq == null) throw new Exception("Can not run refStrFromMd with a null seq variable");
				if (md == null) throw new Exception("Can not run refStrFromMd with a null MD variable. Reads is: " + readName);
				
				// Use sb as the reference output string
				StringBuilder sb = new StringBuilder(500*2+1);

				Matcher match = mdPat.matcher(md);
				int curSeqPos = 0;
				//int curMdPos = 0; // Not the same as seq pos when you have indels

				int savedBases = 0;
				for (final CigarElement cigEl : cigar.getCigarElements()) 
				{
					int cigElLen = cigEl.getLength();
					CigarOperator cigElOp = cigEl.getOperator();
//					System.err.printf("\tCigar El: len=%d, op=%s, consumesRead=%b, consumesRef=%b\n",
//							cigElLen,cigElOp,cigElOp.consumesReadBases(), cigElOp.consumesReferenceBases());
					
					
					// If it consumes reference bases, it's either a match or a deletion in the sequence
					// read.  Either way, we're going to need to parse throught the MD.
					if (cigElOp.consumesReferenceBases())
					{
						// We have a match region, go through the MD
						int basesMatched = 0;
						
						// Do we have any saved matched bases?
						while ((savedBases>0) && (basesMatched < cigElLen))
						{
							sb.append(seq.charAt(curSeqPos++));
							savedBases--;
							basesMatched++;
//							System.err.printf("\t\tDepleting saved bases, saved=%d, curSeqPos=%d, basesMatched=%d\n",savedBases,curSeqPos,basesMatched); 
						}

						while (basesMatched < cigElLen)
						{
							boolean matched = match.find();
							if (matched)
							{
//								System.err.println("Matched , basesMatched=" + basesMatched + ", match=" + match.group() + "," + match.group(1) + "," + match.group(2) + ", start=" + match.start());
								String mg;
								if ( ((mg = match.group(1)) !=null) && (mg.length() > 0) )
								{
									// It's a number , meaning a series of matches
									int num = Integer.parseInt(mg);
									for (int i = 0; i < num; i++)
									{
										if (basesMatched<cigElLen)
										{
											sb.append(seq.charAt(curSeqPos));
											curSeqPos++;
										}
										else
										{
											savedBases++;
										}
										basesMatched++;
									}
								}

								else if ( ((mg = match.group(2)) !=null) && (mg.length() > 0) )
								{
									// It's a single nucleotide, meaning a mismatch
									if (basesMatched<cigElLen)
									{
										sb.append(mg.charAt(0));
										curSeqPos++;
									}
									else
									{
										savedBases++;
									}
									basesMatched++;
								}
								else if ( ((mg = match.group(3)) !=null) && (mg.length() > 0) )
								{
									// It's a deletion, starting with a caret
									// don't include caret
									for (int i = 1; i < mg.length(); i++)
									{
										// Since this function is actually just meant to make a reference that lines up nucleotide 
										//  for nucleotide with the sequence read, we don't actually add the insertion to the reference.
										//sb.append(mg.charAt(i));
										basesMatched++;
									}
									
									// Check just to make sure.
									if (basesMatched != cigElLen)
									{
										throw new Exception("Got a deletion in CIGAR (" + cigar + ", deletion " + cigElLen + 
												" length) with an unequal ref insertion in MD (" + md + ", md " + basesMatched + " length" + "in Read: " + readName);
									}
									if (cigElOp != CigarOperator.DELETION)
									{
										throw new Exception ("Got an insertion in MD ("+md+") without a corresponding deletion in cigar ("+cigar+")" + "in Read: " + readName);
									}
									
								}
								else
								{
									matched = false;
								}
							}

							if (!matched)
							{
								throw new Exception("Illegal MD pattern: " + md + "in Read: " + readName);
							}

//							System.err.println("SavedBasesMatched=" + savedBases);
						}

					}
					else if (cigElOp.consumesReadBases())
					{
						// We have an insertion in read
						for (int i = 0; i < cigElLen; i++)
						{
							char c = (cigElOp == CigarOperator.SOFT_CLIP) ? '0' : '-';
							sb.append( c );
							curSeqPos++;
						}
					}
					else
					{
						// It's an op that consumes neither read nor reference bases.  Do we just ignore??
					}

				}
				
				return sb.toString().getBytes();
			}
	
	public static THashMap<String, AlignmentContext> splitContextBySampleName(AlignmentContext context) {
        GenomeLoc loc = context.getLocation();
        THashMap<String, AlignmentContext> contexts = new THashMap<String, AlignmentContext>();

        for(String sample: context.getPileup().getSamples()) {
            ReadBackedPileup pileupBySample = context.getPileup().getPileupForSample(sample);

            // Don't add empty pileups to the split context.
            if(pileupBySample.getNumberOfElements() == 0)
                continue;

            if(sample != null)
                contexts.put(sample, new AlignmentContext(loc, pileupBySample));
            else {
                    throw new UserException.ReadMissingReadGroup(pileupBySample.iterator().next().getRead());
                
            }
        }

        return contexts;
    }
	
}
