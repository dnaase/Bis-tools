package main.java.edu.usc.epigenome.uecgatk.bissnp.writer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;

import htsjdk.tribble.TribbleException;
import htsjdk.tribble.index.DynamicIndexCreator;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.samtools.util.PositionalOutputStream;

public abstract class FormatWriterBase {
	protected DynamicIndexCreator indexer = null;
	protected File indexFile = null;
	protected BufferedWriter mWriter;
	protected PositionalOutputStream positionalStream = null;
	LittleEndianOutputStream idxStream = null;
	File location = null;

	public FormatWriterBase(File location) {
		this(location, openOutputStream(location));
	}

	public FormatWriterBase(File location, OutputStream output) {
		this.location = location;

		positionalStream = new PositionalOutputStream(output);
		output = positionalStream;

		mWriter = new BufferedWriter(new OutputStreamWriter(output));

	}

	public FormatWriterBase(OutputStream output) {
		mWriter = new BufferedWriter(new OutputStreamWriter(output));
	}

	protected static OutputStream openOutputStream(File location) {
		try {
			return new FileOutputStream(location);
		} catch (FileNotFoundException e) {
			throw new TribbleException("Unable to create bed file at location: " + location);
		}
	}

	public abstract void add(genomeObject obj);

	public abstract void addHeader(Object o);

	public void close() {
		// try to close the vcf stream
		try {
			mWriter.flush();
			mWriter.close();
		} catch (IOException e) {
			throw new TribbleException("Unable to close " + locationString() + " because of " + e.getMessage());
		}

	}

	public void writerFlush() {
		try {
			mWriter.flush();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private String locationString() {
		return location == null ? mWriter.toString() : location.getAbsolutePath();
	}

}
