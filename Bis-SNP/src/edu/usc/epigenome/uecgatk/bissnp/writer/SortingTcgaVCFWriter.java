package edu.usc.epigenome.uecgatk.bissnp.writer;

import org.broadinstitute.sting.utils.variantcontext.VariantContext;

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

public class SortingTcgaVCFWriter extends SortingVCFWriterOwnBase {
	protected TcgaVCFWriter tcgaInnerWriter = null;
	// the maximum START distance between records that we'll cache
	private int maxCachingStartDistance;

	public SortingTcgaVCFWriter(TcgaVCFWriter innerWriter, int maxCachingStartDistance) {
		this(innerWriter, maxCachingStartDistance, false);
		tcgaInnerWriter = innerWriter;
		this.maxCachingStartDistance = maxCachingStartDistance;
		// TODO Auto-generated constructor stub
	}

	/**
	 * create a local-sorting VCF writer, given an inner VCF writer to write to
	 * 
	 * @param innerWriter
	 *            the VCFWriter to write to
	 * @param maxCachingStartDistance
	 *            the maximum start distance between records that we'll cache
	 * @param takeOwnershipOfInner
	 *            Should this Writer close innerWriter when it's done with it
	 */
	public SortingTcgaVCFWriter(TcgaVCFWriter innerWriter, int maxCachingStartDistance, boolean takeOwnershipOfInner) {
		super(innerWriter, takeOwnershipOfInner);
		tcgaInnerWriter = innerWriter;
		this.maxCachingStartDistance = maxCachingStartDistance;
		// TODO Auto-generated constructor stub
	}

	public void enableDiscreteLoci(boolean enableDiscreteLoci) {
		this.enableDiscreteLoci = enableDiscreteLoci;
	}

	public TcgaVCFWriter getInnerWriter() {
		return this.tcgaInnerWriter;
	}

	public void writerFlush() {
		this.tcgaInnerWriter.writeFlush();
	}

	@Override
	protected void noteCurrentRecord(VariantContext vc) {
		// super.noteCurrentRecord(vc); // first, check for errors

		// then, update mostUpstreamWritableLoc:
		
		int mostUpstreamWritableIndex = vc.getStart() - maxCachingStartDistance;
		this.mostUpstreamWritableLoc = Math.max(BEFORE_MOST_UPSTREAM_LOC, mostUpstreamWritableIndex);
		//System.err.println( mostUpstreamWritableIndex + "\t" + mostUpstreamWritableLoc + "\t" + vc.getStart());
	}

}
