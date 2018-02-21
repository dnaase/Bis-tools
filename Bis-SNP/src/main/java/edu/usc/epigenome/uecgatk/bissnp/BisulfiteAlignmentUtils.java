package main.java.edu.usc.epigenome.uecgatk.bissnp;

import java.util.BitSet;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.sam.AlignmentUtils.MismatchCount;

import main.java.edu.usc.epigenome.uecgatk.bissnp.BisulfiteEnums.MethylSNPModel;

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

public class BisulfiteAlignmentUtils {

	public BisulfiteAlignmentUtils() { }

	public static MismatchCount getMismatchCount(SAMRecord r, byte[] refSeq, int refIndex) {
		return getMismatchCount(r, refSeq, refIndex, 0, r.getReadLength());
	}

	public static MismatchCount getMismatchCount(SAMRecord r, byte[] refSeq, int refIndex, int startOnRead, int nReadBases) {
		MismatchCount mc = new MismatchCount();

		int readIdx = 0;
		// index of the last base on read we want to count
		int endOnRead = startOnRead + nReadBases - 1;
		byte[] readSeq = r.getReadBases();
		Cigar c = r.getCigar();
		for (int i = 0; i < c.numCigarElements(); i++) {

			if (readIdx > endOnRead)
				break;

			CigarElement ce = c.getCigarElement(i);
			switch (ce.getOperator()) {

			case M:
				for (int j = 0; j < ce.getLength(); j++, refIndex++, readIdx++) {
					if (refIndex >= refSeq.length)
						continue;
					if (readIdx < startOnRead)
						continue;
					if (readIdx > endOnRead)
						break;
					byte refChr = refSeq[refIndex];
					byte readChr = readSeq[readIdx];
					// Note: we need to count X/N's as mismatches because that's
					// what SAM requires
					if (readChr != refChr) {
						
						if (BaseUtilsMore.isBisulfiteMismatch(refChr, readChr, r.getReadNegativeStrandFlag(), r.getReadPairedFlag() && r.getSecondOfPairFlag())) {
							mc.numMismatches++;
							mc.mismatchQualities += r.getBaseQualities()[readIdx];
						}

					}
				}
				break;
			case I:
			case S:
				readIdx += ce.getLength();
				break;
			case D:
			case N:
				refIndex += ce.getLength();
				break;
			case H:
			case P:
				break;
			default:
				throw new ReviewedGATKException("The " + ce.getOperator() + " cigar element is not currently supported");
			}

		}
		//System.err.println(mc.numMismatches + "\t" + mc.mismatchQualities);
		return mc;
	}

	/**
	 * Returns the number of mismatches in the pileup element within the given reference context in
	 * bisulfite-seq space.
	 * 
	 * @param read
	 *            the SAMRecord
	 * @param ref
	 *            the reference context
	 * @param maxMismatches
	 *            the maximum number of surrounding mismatches we tolerate to consider a base good
	 * @param windowSize
	 *            window size (on each side) to test
	 * @param sequencingMode
	 *            in Bisulfite mode, GNOMe-seq mode or Normal-seq mode
	 * @param pairedend
	 *            paired-end reads or not
	 * @return a bitset representing which bases are good
	 */
	public static BitSet mismatchesInRefWindow(SAMRecord read, ReferenceContext ref, int maxMismatches, int windowSize, MethylSNPModel sequencingMode, boolean pairedend) {
		// first determine the positions with mismatches
		int readLength = read.getReadLength();
		BitSet mismatches = new BitSet(readLength);

		int readStartPos = Math.max(read.getAlignmentStart(), ref.getLocus().getStart() - windowSize);
		int currentReadPos = read.getAlignmentStart();

		byte[] refBases = ref.getBases();
		int refIndex = readStartPos - ref.getWindow().getStart();
		if (refIndex < 0) {
			throw new IllegalStateException("When calculating mismatches, we somehow don't have enough previous reference context for read " + read.getReadName() + " at position " + ref.getLocus());
		}

		byte[] readBases = read.getReadBases();
		int readIndex = 0;

		Cigar c = read.getCigar();
		boolean negStrand = read.getReadNegativeStrandFlag();
		boolean secondPair = false;
		if (pairedend) {
			secondPair = read.getSecondOfPairFlag();
		}

		for (int i = 0; i < c.numCigarElements(); i++) {
			CigarElement ce = c.getCigarElement(i);
			int cigarElementLength = ce.getLength();
			switch (ce.getOperator()) {
			case M:
				for (int j = 0; j < cigarElementLength; j++, readIndex++) {
					// skip over unwanted bases
					if (currentReadPos++ < readStartPos)
						continue;

					// if reads extend beyond the contig end
					if (refIndex >= refBases.length)
						break;

					byte refChr = refBases[refIndex];
					byte readChr = readBases[readIndex];
					if (readChr != refChr) {
						if (sequencingMode == MethylSNPModel.BM || sequencingMode == MethylSNPModel.GM) {
							if (!negStrand) {
								if (secondPair) {
									if (BaseUtils.basesAreEqual(readChr, BaseUtils.Base.G.base) && BaseUtils.basesAreEqual(readChr, BaseUtils.Base.A.base)) {

									} else {
										mismatches.set(readIndex);
									}
								} else {
									if (BaseUtils.basesAreEqual(readChr, BaseUtils.Base.C.base) && BaseUtils.basesAreEqual(readChr, BaseUtils.Base.T.base)) {

									} else {
										mismatches.set(readIndex);
									}
								}

							} else {
								if (secondPair) {
									if (BaseUtils.basesAreEqual(readChr, BaseUtils.Base.C.base) && BaseUtils.basesAreEqual(readChr, BaseUtils.Base.T.base)) {

									} else {
										mismatches.set(readIndex);
									}
								} else {
									if (BaseUtils.basesAreEqual(readChr, BaseUtils.Base.G.base) && BaseUtils.basesAreEqual(readChr, BaseUtils.Base.A.base)) {

									} else {
										mismatches.set(readIndex);
									}
								}

							}
						} else {
							mismatches.set(readIndex);
						}

					}

					refIndex++;
				}
				break;
			case I:
			case S:
				readIndex += cigarElementLength;
				break;
			case D:
			case N:
				if (currentReadPos >= readStartPos)
					refIndex += cigarElementLength;
				currentReadPos += cigarElementLength;
				break;
			case H:
			case P:
				break;
			}
		}

		// all bits are set to false by default
		BitSet result = new BitSet(readLength);

		int currentPos = 0, leftPos = 0, rightPos;
		int mismatchCount = 0;

		// calculate how many mismatches exist in the windows to the left/right
		for (rightPos = 1; rightPos <= windowSize && rightPos < readLength; rightPos++) {
			if (mismatches.get(rightPos))
				mismatchCount++;
		}
		if (mismatchCount <= maxMismatches)
			result.set(currentPos);

		// now, traverse over the read positions
		while (currentPos < readLength) {
			// add a new rightmost position
			if (rightPos < readLength && mismatches.get(rightPos++))
				mismatchCount++;
			// re-penalize the previous position
			if (mismatches.get(currentPos++))
				mismatchCount++;
			// don't penalize the current position
			if (mismatches.get(currentPos))
				mismatchCount--;
			// subtract the leftmost position
			if (leftPos < currentPos - windowSize && mismatches.get(leftPos++))
				mismatchCount--;

			if (mismatchCount <= maxMismatches)
				result.set(currentPos);

		}

		return result;
	}

	public static long mismatchingQualities(SAMRecord r, byte[] refSeq, int refIndex) {
		return getMismatchCount(r, refSeq, refIndex).mismatchQualities;
	}
}
