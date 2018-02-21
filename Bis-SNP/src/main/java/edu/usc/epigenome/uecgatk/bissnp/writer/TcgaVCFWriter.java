package main.java.edu.usc.epigenome.uecgatk.bissnp.writer;

import java.io.*;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import htsjdk.samtools.util.LocationAware;
import htsjdk.samtools.util.PositionalOutputStream;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.index.DynamicIndexCreator;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexCreator;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import htsjdk.variant.variantcontext.writer.SortingVariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;

import main.java.edu.usc.epigenome.uecgatk.bissnp.BisSNP;
import main.java.edu.usc.epigenome.uecgatk.bissnp.BisSNPUtils;
import main.java.edu.usc.epigenome.uecgatk.bissnp.BisulfiteVCFConstants;

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
public class TcgaVCFWriter implements VariantContextWriter {

	private SAMSequenceDictionary refDict = null;
	private boolean useOriginalHeader = false;
	private GenomeAnalysisEngine toolkit = null;
	private Writer mWriter;
	protected VCFHeader mHeader = null;
	private OutputStream outputStream;
	private boolean doNotWriteGenotypes = false;
	private boolean allowMissingFieldsInHeader;
	private boolean writeFullFormatField;
	protected boolean filtersWereAppliedToContext = false;
	private IndexCreator indexer;
	private File location;
	private String name;
	private LocationAware locationSource;
	private static int INITIAL_BUFFER_SIZE = 16384;
	private final ByteArrayOutputStream lineBuffer = new ByteArrayOutputStream(16384);
	private VCFEncoder vcfEncoder = null;
	
	public TcgaVCFWriter(File location, OutputStream output, SAMSequenceDictionary refDict,
						 boolean enableOnTheFlyIndexing, boolean doNotWriteGenotypes,
						 boolean allowMissingFieldsInHeader, boolean writeFullFormatField) {
		//super(location, output, refDict, enableOnTheFlyIndexing, doNotWriteGenotypes);
		this.name = writerName(location, output);
		this.locationSource = null;
		this.indexer = null;
		this.location = location;
		this.outputStream = output;
		this.refDict = refDict;
		if (enableOnTheFlyIndexing) {
			this.initIndexingWriter(new DynamicIndexCreator(location, IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME));
		}
		this.mWriter = new BufferedWriter(new OutputStreamWriter(this.lineBuffer, VCFEncoder.VCF_CHARSET));
		this.doNotWriteGenotypes = doNotWriteGenotypes;
		this.allowMissingFieldsInHeader = allowMissingFieldsInHeader;
		this.writeFullFormatField = writeFullFormatField;
	}

	public TcgaVCFWriter(File location, OutputStream output, SAMSequenceDictionary refDict, IndexCreator indexCreator,
						 boolean enableOnTheFlyIndexing, boolean doNotWriteGenotypes, boolean allowMissingFieldsInHeader,
						 boolean writeFullFormatField) {
		this.name = writerName(location, output);
		this.locationSource = null;
		this.indexer = null;
		this.location = location;
		this.outputStream = output;
		this.refDict = refDict;
		if (enableOnTheFlyIndexing) {
			this.initIndexingWriter(indexCreator);
		}

		this.mWriter = new BufferedWriter(new OutputStreamWriter(this.lineBuffer, VCFEncoder.VCF_CHARSET));
		this.doNotWriteGenotypes = doNotWriteGenotypes;
		this.allowMissingFieldsInHeader = allowMissingFieldsInHeader;
		this.writeFullFormatField = writeFullFormatField;
	}

	// get the system time that this VCF file generated

	private void initIndexingWriter(IndexCreator idxCreator) {
		this.indexer = idxCreator;
		if (this.outputStream instanceof LocationAware) {
			this.locationSource = (LocationAware)this.outputStream;
		} else {
			PositionalOutputStream positionalOutputStream = new PositionalOutputStream(this.outputStream);
			this.locationSource = positionalOutputStream;
			this.outputStream = positionalOutputStream;
		}

	}

	public void writeFlush() {

		try {
			mWriter.flush();
			this.getOutputStream().write(this.lineBuffer.toByteArray());
			this.lineBuffer.reset();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			throw new TribbleException("IOException writing the VCF flush to " + e);
		}

	}

	@Override
	public boolean checkError() {
		return this.getOutputStream() instanceof PrintStream && ((PrintStream)PrintStream.class.cast(this.getOutputStream())).checkError();
	}

	public OutputStream getOutputStream() {
		return this.outputStream;
	}

	 @Override
	public void writeHeader(VCFHeader header) {
		 
		mHeader = doNotWriteGenotypes ? new VCFHeader(header.getMetaDataInInputOrder()) : header;

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
				mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_CENTER, "none").toString() + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_PHASE, "none").toString() + "\n");
				mWriter.write(VCFHeader.METADATA_INDICATOR + new VCFHeaderLine(BisulfiteVCFConstants.VCF_HEADER_VERSION_GAF, "none").toString() + "\n");
				
				for (VCFHeaderLine line : mHeader.getMetaDataInInputOrder()) {
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
				
			
				
				
				for (VCFHeaderLine line : mHeader.getMetaDataInInputOrder()) {
					
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
			writeFlush(); // necessary so that writing to an output stream
								// will work
			this.vcfEncoder = new VCFEncoder(this.mHeader, this.allowMissingFieldsInHeader, this.writeFullFormatField);

		} catch (IOException e) {
			throw new TribbleException("IOException writing the VCF header to " + e);
		}
	}

	public void close() {
		try {
			this.mWriter.close();
			this.outputStream.close();
			if (this.indexer != null) {
				this.indexer.setIndexSequenceDictionary(this.refDict);
				Index index = this.indexer.finalizeIndex(this.locationSource.getPosition());
				index.writeBasedOnFeatureFile(this.location);
			}

		} catch (IOException var2) {
			throw new RuntimeIOException("Unable to close index for " + this.getStreamName(), var2);
		}
	}

	public void add(VariantContext vc) {
		if (this.indexer != null) {
			this.indexer.addFeature(vc, this.locationSource.getPosition());
		}
		try {
			if (this.doNotWriteGenotypes) {
				mWriter.write(this.vcfEncoder.encode((new VariantContextBuilder(vc)).noGenotypes().make()));
			} else {
				mWriter.write(this.vcfEncoder.encode(vc));
			}

			mWriter.write("\n");
			this.writeFlush();
		} catch (IOException var3) {
			throw new RuntimeIOException("Unable to write the VCF object to " + this.getStreamName(), var3);
		}

	}

	public String getStreamName() {
		return this.name;
	}

	protected static final String writerName(File location, OutputStream stream) {
		return location == null ? stream.toString() : location.getAbsolutePath();
	}
}
