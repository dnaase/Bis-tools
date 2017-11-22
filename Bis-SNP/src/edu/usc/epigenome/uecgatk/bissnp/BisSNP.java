package edu.usc.epigenome.uecgatk.bissnp;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.ResourceBundle;

import net.sf.picard.filter.SamRecordFilter;

import org.broadinstitute.sting.commandline.*;

import org.broadinstitute.sting.gatk.CommandLineExecutable;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.gatk.io.stubs.OutputStreamArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.SAMFileReaderArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.VCFWriterArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.ApplicationDetails;
import org.broadinstitute.sting.utils.text.ListFileUtils;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;

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
  
public class BisSNP extends CommandLineExecutable {

	private static String argCommandline = "";
  
	private static String BisVersion = "BisSNP-0.91";
 
	public Walker<?, ?> walker = null;

	@Argument(fullName = "analysis_type", shortName = "T", doc = "Type of analysis to run")
	private String analysisName = null;

	// to record it is in second iteration or not
	// private static boolean secondIteration = false;

	// argument collection, the collection of command line args we accept. copy
	// from GATK, since they are private class in GATK
	@ArgumentCollection
	private GATKArgumentCollection argCollection = new GATKArgumentCollection();

	private final Collection<Object> bisulfiteArgumentSources = new ArrayList<Object>();

	public static List<String> createApplicationHeader() {

		List<String> header = new ArrayList<String>();
		header.add(String.format("The %s, Compiled %s", getBisSNPVersionNumber(), getBuildTime()));
		header.add(String
				.format("Based on The Genome Analysis Toolkit (GATK) v%s (prebuild GATK package could be download here: ftp://ftp.broadinstitute.org/pub/gsa/GenomeAnalysisTK/GenomeAnalysisTK-1.5-3-gbb2c10b.tar.bz2)",
						getVersionNumber()));
		header.add("Copyright (c) 2011 USC Epigenome Center");
		header.add("Please view our documentation at http://epigenome.usc.edu/publicationdata/bissnp2011/");
		header.add("For support, please send email to lyping1986@gmail.com or benbfly@gmail.com");
		return header;
	}

	public static String getBisSNPArgumantsInput() {
		return argCommandline;
	}

	public static String getBisSNPVersionNumber() {
		return BisVersion;
	}

	public static String getBuildTime() {
		ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText");
		return headerInfo.containsKey("build.timestamp") ? headerInfo.getString("build.timestamp") : "<unknown>";
	}

	public static String getVersionNumber() {
		ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText");

		return headerInfo.containsKey("org.broadinstitute.sting.gatk.version") ? headerInfo.getString("org.broadinstitute.sting.gatk.version") : "<unknown>";
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		try {
			BisSNP instance = new BisSNP();
			for (String arg : args) {
				argCommandline = argCommandline + " " + arg;
			}
			start(instance, args);

			System.exit(CommandLineProgram.result);
		} catch (UserException e) {
			exitSystemWithUserError(e);
		} catch (Exception e) {
			exitSystemWithError(e);
		}

	}

	@Override
	public String getAnalysisName() {
		// TODO Auto-generated method stub
		return analysisName;
	}

	@Override
	protected int execute() throws Exception {

		try {

			engine.setParser(parser);
			bisulfiteArgumentSources.add(this);

			walker = engine.getWalkerByName(getAnalysisName());

			engine.setArguments(getArgumentCollection());

			engine.setSAMFileIDs(ListFileUtils.unpackBAMFileList(getArgumentCollection().samFiles, parser));

			engine.setWalker(walker);
			walker.setToolkit(engine);

			Collection<ReadFilter> filters = engine.createFilters();
			engine.setFilters(filters);

			loadArgumentsIntoObject(walker);
			bisulfiteArgumentSources.add(walker);

			Collection<RMDTriplet> rodBindings = ListFileUtils.unpackRODBindings(parser.getRodBindings(), parser);

			

			engine.setReferenceMetaDataFiles(rodBindings);

			for (SamRecordFilter filter : filters) {
				loadArgumentsIntoObject(filter);
				bisulfiteArgumentSources.add(filter);
			}

			engine.execute();

		} catch (Exception e) {
			throw e;
		}
		return 0;
	}

	@Override
	protected ApplicationDetails getApplicationDetails() {
		return new ApplicationDetails(createApplicationHeader(), Collections.<String> emptyList(), ApplicationDetails.createDefaultRunningInstructions(getClass()), null);
	}

	@Override
	protected GATKArgumentCollection getArgumentCollection() {
		// TODO Auto-generated method stub
		return argCollection;
	}

	/**
	 * Subclasses of CommandLinePrograms can provide their own types of command-line arguments.
	 * 
	 * @return A collection of type descriptors generating implementation-dependent placeholders.
	 */
	@Override
	protected Collection<ArgumentTypeDescriptor> getArgumentTypeDescriptors() {
		return Arrays.asList(new VCFWriterArgumentTypeDescriptor(engine, System.out, bisulfiteArgumentSources), new SAMFileReaderArgumentTypeDescriptor(engine),
				new SAMFileWriterArgumentTypeDescriptor(engine, System.out), new OutputStreamArgumentTypeDescriptor(engine, System.out));
	}

}
