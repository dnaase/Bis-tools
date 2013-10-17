package edu.usc.epigenome.uecgatk.bissnp.writer;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import org.broad.tribble.TribbleException;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.codecs.vcf.StandardVCFWriter;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFilterHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderVersion;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;

import edu.usc.epigenome.uecgatk.bissnp.BisSNP;
import edu.usc.epigenome.uecgatk.bissnp.BisSNPUtils;
import edu.usc.epigenome.uecgatk.bissnp.BisulfiteVCFConstants;

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
 * VCF Writer to generate TCGA specific VCF file, it is only for sorted coordinate now, so in
 * multi-thread mode, this could not be used
 */
public class TcgaVCFWriter extends StandardVCFWriter {

	private SAMSequenceDictionary refDict = null;
	private boolean useOriginalHeader = false;
	private GenomeAnalysisEngine toolkit = null;
	
	public TcgaVCFWriter(File location, OutputStream output, SAMSequenceDictionary refDict, boolean enableOnTheFlyIndexing, boolean doNotWriteGenotypes) {
		super(location, output, refDict, enableOnTheFlyIndexing, doNotWriteGenotypes);
		this.refDict = refDict;
		// TODO Auto-generated constructor stub
	}

	public TcgaVCFWriter(File location, SAMSequenceDictionary refDict) {
		super(location, refDict);
		this.refDict = refDict;
		// TODO Auto-generated constructor stub
	}

	public TcgaVCFWriter(File location, SAMSequenceDictionary refDict, boolean enableOnTheFlyIndexing) {
		super(location, refDict, enableOnTheFlyIndexing);
		this.refDict = refDict;
		// TODO Auto-generated constructor stub
	}
	
	public TcgaVCFWriter(File location,SAMSequenceDictionary refDict, GenomeAnalysisEngine genomeAnalysisEngine, boolean useOriginalHeader) {
		
		super(location, refDict);
		this.refDict = refDict;
		this.useOriginalHeader = useOriginalHeader;
		this.toolkit = genomeAnalysisEngine;
		// TODO Auto-generated constructor stub
	}

	public TcgaVCFWriter(OutputStream output, SAMSequenceDictionary refDict) {
		super(output, refDict, false);
		this.refDict = refDict;
		// TODO Auto-generated constructor stub
	}

	public TcgaVCFWriter(OutputStream output, SAMSequenceDictionary refDict, boolean doNotWriteGenotypes) {
		super(output, refDict, doNotWriteGenotypes);
		this.refDict = refDict;
		// TODO Auto-generated constructor stub
	}

	// get the system time that this VCF file generated



	public void writeFlush() {

		try {
			mWriter.flush();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			throw new TribbleException("IOException writing the VCF flush to " + e);
		}

	}

	 @Override
	public void writeHeader(VCFHeader header) {
		 
		mHeader = doNotWriteGenotypes ? new VCFHeader(header.getMetaData()) : header;

		try {
			// the file format field needs to be written first, specially for
			// TCGA VCF header
			if(!useOriginalHeader){
				mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_FORMAT, "VCFv4.1").toString() + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_DATE, BisSNPUtils.now("yyyyMMdd")).toString() + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_TCGA_VERSION, "1.1").toString() + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR
						+ new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_LOG, "<InputVCF=<>, InputVCFSource=<>, InputVCFVer=<1.1>, InputVCFParam=<> InputVCFgeneAnno=<>>").toString() + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_REF, "<ID=" + refDict.getSequence(0).getAssembly() + ",Source=" + refDict.getSequence(0).getAttribute("UR") + ">").toString() + "\n");
				for(SAMSequenceRecord record : refDict.getSequences()){
					mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_CONTIG, "<ID=" + record.getSequenceName() + ",length=" + record.getSequenceLength()
							+ ",assembly=" + record.getAssembly() + ",md5=" + record.getAttribute("M5") + ",species=" + record.getAttribute("SP") + ">").toString() + "\n");
				}
				mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_ASSEMBLY, "none") + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_CENTER, "USC Epigenome Center").toString() + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_PHASE, "none").toString() + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_GAF, "none").toString() + "\n");
				
				for (VCFHeaderLine line : mHeader.getMetaData()) {
					if (line.getKey().equals(VCFHeaderVersion.VCF4_0.getFormatString()) || line.getKey().equals(VCFHeaderVersion.VCF3_3.getFormatString())
							|| line.getKey().equals(VCFHeaderVersion.VCF3_2.getFormatString()))
						continue;

					// are the records filtered (so we know what to put in the
					// FILTER column of passing records) ?
					if (line instanceof VCFFilterHeaderLine)
						filtersWereAppliedToContext = true;

					mWriter.write(VCFHeader.METADATA_INDICATOR);
					mWriter.write(line.toString());
					mWriter.write("\n");
				}
			}
			else{
				mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_FORMAT, "VCFv4.1").toString() + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_DATE, BisSNPUtils.now("yyyyMMdd")).toString() + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + mHeader.getOtherHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_TCGA_VERSION) + "\n");
				//mWriter.write(VCFHeader.METADATA_INDICATOR + mHeader.getOtherHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_LOG) + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR
						+ new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_LOG, "<InputVCF=<" + toolkit.getRodDataSources().get(0).getFile().getAbsolutePath() + ">, InputVCFSource=<>, InputVCFVer=<" + mHeader.getOtherHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_TCGA_VERSION).getValue() 
								+">, InputVCFParam=<> InputVCFgeneAnno=<" + mHeader.getOtherHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_GAF).getValue() + ">>").toString() + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + mHeader.getOtherHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_REF) + "\n");
				for(SAMSequenceRecord record : refDict.getSequences()){
					mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_CONTIG, "<ID=" + record.getSequenceName() + ",length=" + record.getSequenceLength()
							+ ",assembly=" + record.getAssembly() + ",md5=" + record.getAttribute("M5") + ",species=" + record.getAttribute("SP") + ">").toString() + "\n");
				}
				mWriter.write(VCFHeader.METADATA_INDICATOR + mHeader.getOtherHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_ASSEMBLY) + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + mHeader.getOtherHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_CENTER) + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + mHeader.getOtherHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_PHASE) + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + mHeader.getOtherHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_GAF) + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + mHeader.getOtherHeaderLine(BisulfiteVCFConstants.PROGRAM_ARGS) + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.PROGRAM_ARGS, BisSNP.getBisSNPArgumantsInput()).toString() + "\n");
				
			
				
				
				for (VCFHeaderLine line : mHeader.getMetaData()) {
					
					if (line instanceof VCFFilterHeaderLine)
						filtersWereAppliedToContext = true;

					if (line instanceof VCFInfoHeaderLine || line instanceof VCFFormatHeaderLine || line instanceof VCFFilterHeaderLine){
						mWriter.write(VCFHeader.METADATA_INDICATOR);
						mWriter.write(line.toString());
						mWriter.write("\n");
					}
					
					
				}
			}
			
			

			// write out the column line
			mWriter.write(VCFHeader.HEADER_INDICATOR);
			for (VCFHeader.HEADER_FIELDS field : mHeader.getHeaderFields()) {
				mWriter.write(field.toString());
				mWriter.write(VCFConstants.FIELD_SEPARATOR);
			}

			if (mHeader.hasGenotypingData()) {
				mWriter.write("FORMAT");
				for (String sample : mHeader.getGenotypeSamples()) {
					mWriter.write(VCFConstants.FIELD_SEPARATOR);
					mWriter.write(sample);
				}
			}

			mWriter.write("\n");
			mWriter.flush(); // necessary so that writing to an output stream
								// will work
		} catch (IOException e) {
			throw new TribbleException("IOException writing the VCF header to " + e);
		}
	}
}
