package edu.usc.epigenome.uecgatk.bissnp.writer;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

public class cpgReadsWriterImp extends FormatWriterBase {

	private boolean notEncrypt = false;
	private boolean addRefStrand = false;
	
	public cpgReadsWriterImp(File location) {
		super(location);
		// TODO Auto-generated constructor stub
	}

	public cpgReadsWriterImp(File location, OutputStream output) {
		super(location, output);
		// TODO Auto-generated constructor stub
	}

	public cpgReadsWriterImp(OutputStream output) {
		super(output);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void add(genomeObject obj) {

		String readsLine = String.format("%s\t%d\t%c\t%d\t%c\t%s\n", obj.getChr(), obj.getStart(), ((cpgReads) obj).getMethyStatus(), ((cpgReads) obj).getbaseQ(), ((cpgReads) obj).getstrand(),
				((cpgReads) obj).getEncryptID64Ascii());
		if(addRefStrand){
			readsLine = String.format("%s\t%d\t%c\t%d\t%c\t%s\t%s\n", obj.getChr(), obj.getStart(), ((cpgReads) obj).getMethyStatus(), ((cpgReads) obj).getbaseQ(), ((cpgReads) obj).getstrand(),
					((cpgReads) obj).getEncryptID64Ascii(), ((cpgReads) obj).getrefStrand());
		}
		if(this.notEncrypt){
			readsLine = String.format("%s\t%d\t%c\t%d\t%c\t%s\n", obj.getChr(), obj.getStart(), ((cpgReads) obj).getMethyStatus(), ((cpgReads) obj).getbaseQ(), ((cpgReads) obj).getstrand(),
					((cpgReads) obj).getReadID());
			if(addRefStrand){
				readsLine = String.format("%s\t%d\t%c\t%d\t%c\t%s\t%s\n", obj.getChr(), obj.getStart(), ((cpgReads) obj).getMethyStatus(), ((cpgReads) obj).getbaseQ(), ((cpgReads) obj).getstrand(),
						((cpgReads) obj).getReadID(), ((cpgReads) obj).getrefStrand());
			}
		}
		
		
		try {
			mWriter.write(readsLine);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	@Override
	public void addHeader(Object o) {
		// TODO Auto-generated method stub
		String header = "#chr\tpos\tmethyStatus\tbaseQ\tstrand\treadID\n";
		if(addRefStrand){
			header = "#chr\tpos\tmethyStatus\tbaseQ\tstrand\treadID\trefCStrand\n";
		}
		
		try {
			mWriter.write(header);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void setEncrypt(boolean notEncrypt){
		this.notEncrypt = notEncrypt ;
	}
	
	public void addRefStrand(boolean addRefStrand){
		this.addRefStrand = addRefStrand ;
	}

}
