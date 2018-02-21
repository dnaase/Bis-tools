package main.java.edu.usc.epigenome.uecgatk.bissnp.writer;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Iterator;

public class bedObjectWriterImp extends FormatWriterBase {

	public bedObjectWriterImp(File location) {
		super(location);
		// TODO Auto-generated constructor stub
	}

	public bedObjectWriterImp(File location, OutputStream output) {
		super(location, output);
		// TODO Auto-generated constructor stub
	}

	public bedObjectWriterImp(OutputStream output) {
		super(output);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void add(genomeObject obj) {
		// TODO Auto-generated method stub
		String readsLine = String.format("%s\t%d\t%d\t%s", obj.getChr(), obj.getStart(), ((bedObject) obj).getEnd(), ((bedObject) obj).getStrandAsString());

		try {
			mWriter.write(readsLine);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Iterator<Object> it = ((bedObject) obj).getValueObject().iterator();

		while (it.hasNext()) {

			try {
				mWriter.write("\t" + it.next());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		try {
			mWriter.write("\n");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	@Override
	public void addHeader(Object o) {
		// TODO Auto-generated method stub
		try {
			mWriter.write(o.toString());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
