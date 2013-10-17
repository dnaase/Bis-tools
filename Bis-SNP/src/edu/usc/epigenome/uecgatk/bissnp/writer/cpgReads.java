package edu.usc.epigenome.uecgatk.bissnp.writer;

import jonelo.jacksum.algorithm.Crc64;
import jonelo.jacksum.util.Service;

public class cpgReads implements genomeObject {

	private byte baseQ;
	private String chr;
	// private CRC32 encrypt;
	private Crc64 encrypt64;
	private int genomeLoc;
	private char methyStatus;
	private String readID;
	private char strand;
	private String refStrand;
	
	public cpgReads(String chr, int genomeLoc, char methyStatus, byte baseQ, char strand, String readID) {
		this.chr = chr;
		this.genomeLoc = genomeLoc;
		this.methyStatus = methyStatus;
		this.baseQ = baseQ;
		this.strand = strand;
		this.readID = readID;
		// this.encrypt = new CRC32();
		// encrypt.update(readID.getBytes());
		this.encrypt64 = new Crc64();
		encrypt64.update(readID.getBytes());
	}
	
	public cpgReads(String chr, int genomeLoc, char methyStatus, byte baseQ, char strand, String readID, String refStrand) {
		this.chr = chr;
		this.genomeLoc = genomeLoc;
		this.methyStatus = methyStatus;
		this.baseQ = baseQ;
		this.strand = strand;
		this.readID = readID;
		// this.encrypt = new CRC32();
		// encrypt.update(readID.getBytes());
		this.encrypt64 = new Crc64();
		encrypt64.update(readID.getBytes());
		this.refStrand=refStrand;
	}

	public byte getbaseQ() {
		return this.baseQ;
	}

	@Override
	public String getChr() {
		// TODO Auto-generated method stub
		return this.chr;
	}

	public long getEncryptID64() {
		return this.encrypt64.getValue();
	}

	public String getEncryptID64Ascii() {

		return Service.format(this.encrypt64.getByteArray());
	}

	public char getMethyStatus() {
		return this.methyStatus;
	}

	public String getReadID() {
		return this.readID;
	}

	@Override
	public int getStart() {
		// TODO Auto-generated method stub
		return this.genomeLoc;
	}

	public char getstrand() {
		return this.strand;
	}
	
	public String getrefStrand() {
		return this.refStrand;
	}

}
