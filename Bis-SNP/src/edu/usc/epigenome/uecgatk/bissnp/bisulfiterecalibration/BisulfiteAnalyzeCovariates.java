/**
 * 
 */
package edu.usc.epigenome.uecgatk.bissnp.bisulfiterecalibration;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.regex.Pattern;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.analyzecovariates.AnalyzeCovariates;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.gatk.walkers.recalibration.Covariate;
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum;
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.R.RScriptExecutor;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.io.Resource;
import org.broadinstitute.sting.utils.text.XReadLines;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Apr 13, 2012 7:52:10 PM
 * 
 */
@DocumentedGATKFeature(groupName = "BisulfiteAnalyzeCovariates", summary = "Package to plot residual accuracy versus error covariates for the base quality score recalibrator")
public class BisulfiteAnalyzeCovariates extends CommandLineProgram {
	protected static final String EOF_MARKER = "EOF";

	final private static Logger logger = Logger.getLogger(AnalyzeCovariates.class);
	private static final String PLOT_INDEL_QUALITY_RSCRIPT = "plot_indelQuality.R";
	private static final String PLOT_RESDIUAL_ERROR_OTHER_COVARIATE = "plot_residualError_OtherCovariate.R";

	private static final String PLOT_RESDIUAL_ERROR_QUALITY_SCORE_COVARIATE = "plot_residualError_QualityScoreCovariate.R";
	private final Pattern COMMENT_PATTERN = Pattern.compile("^#.*");
	private final Pattern COVARIATE_PATTERN = Pattern.compile("^ReadGroup,QualityScore,.*");
	// ///////////////////////////
	// Private Member Variables
	// ///////////////////////////
	private BisulfiteAnalysisDataManager dataManager;

	@Hidden
	@Argument(fullName = "do_indel_quality", shortName = "indels", required = false, doc = "If supplied, do indel quality plotting")
	private boolean DO_INDEL_QUALITY = false;

	@Argument(fullName = "ignoreQ", shortName = "ignoreQ", doc = "Ignore bases with reported quality less than this number.", required = false)
	private int IGNORE_QSCORES_LESS_THAN = 5;

	/**
	 * This argument is useful for comparing before/after plots and you want the axes to match each
	 * other.
	 */
	@Argument(fullName = "max_histogram_value", shortName = "maxHist", required = false, doc = "If supplied, this value will be the max value of the histogram plots")
	private int MAX_HISTOGRAM_VALUE = 0;

	/**
	 * Combinations of covariates in which there are zero mismatches technically have infinite
	 * quality. We get around this situation by capping at the specified value. We've found that Q40
	 * is too low when using a more completely database of known variation like dbSNP build 132 or
	 * later.
	 */
	@Argument(fullName = "max_quality_score", shortName = "maxQ", required = false, doc = "The integer value at which to cap the quality scores, default is 50")
	private int MAX_QUALITY_SCORE = 50;
	@Argument(fullName = "numRG", shortName = "numRG", doc = "Only process N read groups. Default value: -1 (process all read groups)", required = false)
	private int NUM_READ_GROUPS_TO_PROCESS = -1; // -1 means process all read
													// groups
	private final Pattern OLD_RECALIBRATOR_HEADER = Pattern.compile("^rg,.*");
	@Argument(fullName = "output_dir", shortName = "outputDir", doc = "The directory in which to output all the plots and intermediate data files", required = false)
	private File OUTPUT_DIR = new File("analyzeCovariates");
	// ///////////////////////////
	// Command Line Arguments
	// ///////////////////////////
	/**
	 * After the header, data records occur one per line until the end of the file. The first
	 * several items on a line are the values of the individual covariates and will change depending
	 * on which covariates were specified at runtime. The last three items are the data- that is,
	 * number of observations for this combination of covariates, number of reference mismatches,
	 * and the raw empirical quality score calculated by phred-scaling the mismatch rate.
	 */
	@Input(fullName = "recal_file", shortName = "recalFile", doc = "The input recal csv file to analyze", required = false)
	private String RECAL_FILE = "output.recal_data.csv";
	private ArrayList<Covariate> requestedCovariates;

	public static void main(String args[]) {
		try {
			BisulfiteAnalyzeCovariates clp = new BisulfiteAnalyzeCovariates();
			start(clp, args);
			System.exit(CommandLineProgram.result);
		} catch (Exception e) {
			exitSystemWithError(e);
		}
	}

	@Override
	protected int execute() {

		// create the output directory where all the data tables and plots will
		// go
		if (!OUTPUT_DIR.exists() && !OUTPUT_DIR.mkdirs())
			throw new UserException.BadArgumentValue("--output_dir/-outDir", "Unable to create output directory: " + OUTPUT_DIR);

		if (!RScriptExecutor.RSCRIPT_EXISTS)
			Utils.warnUser(logger, "Rscript not found in environment path. Plots will not be generated.");

		// initialize all the data from the csv file and allocate the list of
		// covariates
		logger.info("Reading in input csv file...");
		initializeData();
		logger.info("...Done!");

		// output data tables for Rscript to read in
		logger.info("Writing out intermediate tables for R...");
		writeDataTables();
		logger.info("...Done!");

		// perform the analysis using Rscript and output the plots
		logger.info("Calling analysis R scripts and writing out figures...");
		callRScripts();
		logger.info("...Done!");

		return 0;
	}

	private void addCSVData(String line) {
		String[] vals = line.split(",");

		// Check if the data line is malformed, for example if the read group
		// string contains a comma then it won't be parsed correctly
		if (vals.length != requestedCovariates.size() + 3) {
			throw new RuntimeException("Malformed input recalibration file. Found data line with too many fields: " + line
					+ " --Perhaps the read group string contains a comma and isn't being parsed correctly.");
		}

		Object[] key = new Object[requestedCovariates.size()];
		Covariate cov;
		int iii;
		for (iii = 0; iii < requestedCovariates.size(); iii++) {
			cov = requestedCovariates.get(iii);
			// key[iii] = cov.getValue( vals[iii] );
			if (cov instanceof BisulfiteDinucCovariate) {
				key[iii] = ((BisulfiteDinucCovariate) cov).getValueBisulfite(vals[iii]);
			} else {
				key[iii] = cov.getValue(vals[iii]);
			}
		}
		// Create a new datum using the number of observations, number of
		// mismatches, and reported quality score
		RecalDatum datum = new RecalDatum(Long.parseLong(vals[iii]), Long.parseLong(vals[iii + 1]), Double.parseDouble(vals[1]), 0.0);
		// Add that datum to all the collapsed tables which will be used in the
		// sequential calculation
		dataManager.addToAllTables(key, datum, IGNORE_QSCORES_LESS_THAN);
	}

	private void callRScripts() {
		int numReadGroups = 0;

		// for each read group
		for (Object readGroupKey : dataManager.getCollapsedTable(0).data.keySet()) {
			if (++numReadGroups <= NUM_READ_GROUPS_TO_PROCESS || NUM_READ_GROUPS_TO_PROCESS == -1) {

				String readGroup = readGroupKey.toString();
				logger.info("Analyzing read group: " + readGroup);

				// for each covariate
				for (int iii = 1; iii < requestedCovariates.size(); iii++) {
					final Covariate cov = requestedCovariates.get(iii);
					final File outputFile = new File(OUTPUT_DIR, readGroup + "." + cov.getClass().getSimpleName() + ".dat");
					if (DO_INDEL_QUALITY) {
						RScriptExecutor executor = new RScriptExecutor();
						executor.addScript(new Resource(PLOT_INDEL_QUALITY_RSCRIPT, AnalyzeCovariates.class));
						// The second argument is the name of the covariate in
						// order to make the plots look nice
						executor.addArgs(outputFile, cov.getClass().getSimpleName().split("Covariate")[0]);
						executor.exec();
					} else {
						if (iii == 1) {
							// Analyze reported quality
							RScriptExecutor executor = new RScriptExecutor();
							executor.addScript(new Resource(PLOT_RESDIUAL_ERROR_QUALITY_SCORE_COVARIATE, AnalyzeCovariates.class));
							// The second argument is the Q scores that should
							// be turned pink in the plot because they were
							// ignored
							executor.addArgs(outputFile, IGNORE_QSCORES_LESS_THAN, MAX_QUALITY_SCORE, MAX_HISTOGRAM_VALUE);
							executor.exec();
						} else { // Analyze all other covariates
							RScriptExecutor executor = new RScriptExecutor();
							executor.addScript(new Resource(PLOT_RESDIUAL_ERROR_OTHER_COVARIATE, AnalyzeCovariates.class));
							// The second argument is the name of the covariate
							// in order to make the plots look nice
							executor.addArgs(outputFile, cov.getClass().getSimpleName().split("Covariate")[0]);
							executor.exec();
						}
					}
				}
			} else { // at the maximum number of read groups so break out
				break;
			}
		}
	}

	private void initializeData() {

		// Get a list of all available covariates
		Collection<Class<? extends Covariate>> classes = new PluginManager<Covariate>(Covariate.class).getPlugins();

		int lineNumber = 0;
		boolean foundAllCovariates = false;

		// Read in the covariates that were used from the input file
		requestedCovariates = new ArrayList<Covariate>();

		try {
			for (String line : new XReadLines(new File(RECAL_FILE))) {
				lineNumber++;
				if (COMMENT_PATTERN.matcher(line).matches() || OLD_RECALIBRATOR_HEADER.matcher(line).matches() || line.equals(EOF_MARKER)) {
					; // Skip over the comment lines, (which start with '#')
				} else if (COVARIATE_PATTERN.matcher(line).matches()) {
					if (foundAllCovariates) {
						throw new RuntimeException("Malformed input recalibration file. Found covariate names intermingled with data in file: " + RECAL_FILE);
					} else { // Found the covariate list in input file, loop
								// through all of them and instantiate them
						String[] vals = line.split(",");
						for (int iii = 0; iii < vals.length - 3; iii++) {
							boolean foundClass = false;
							for (Class<?> covClass : classes) {
								if ((vals[iii] + "Covariate").equalsIgnoreCase(covClass.getSimpleName())) {
									foundClass = true;
									try {
										Covariate covariate = (Covariate) covClass.newInstance();
										requestedCovariates.add(covariate);
									} catch (Exception e) {
										throw new DynamicClassResolutionException(covClass, e);
									}
								}
							}

							if (!foundClass) {
								throw new RuntimeException("Malformed input recalibration file. The requested covariate type (" + (vals[iii] + "Covariate") + ") isn't a valid covariate option.");
							}
						}

					}

				} else { // Found a line of data
					if (!foundAllCovariates) {

						foundAllCovariates = true;

						// At this point all the covariates should have been
						// found and initialized
						if (requestedCovariates.size() < 2) {
							throw new RuntimeException("Malformed input recalibration file. Covariate names can't be found in file: " + RECAL_FILE);
						}

						// Initialize any covariate member variables using the
						// shared argument collection
						for (Covariate cov : requestedCovariates) {
							cov.initialize(new RecalibrationArgumentCollection());
						}

						// Initialize the data hashMaps
						dataManager = new BisulfiteAnalysisDataManager(requestedCovariates.size());

					}
					addCSVData(line); // Parse the line and add the data to the
										// HashMap
				}
			}

		} catch (FileNotFoundException e) {
			throw new RuntimeException("Can not find input file: " + RECAL_FILE);
		} catch (NumberFormatException e) {
			throw new RuntimeException("Error parsing recalibration data at line " + lineNumber + ". Perhaps your table was generated by an older version of CovariateCounterWalker.");
		}
	}

	private void writeDataTables() {

		int numReadGroups = 0;

		// for each read group
		for (Object readGroupKey : dataManager.getCollapsedTable(0).data.keySet()) {

			if (NUM_READ_GROUPS_TO_PROCESS == -1 || ++numReadGroups <= NUM_READ_GROUPS_TO_PROCESS) {
				String readGroup = readGroupKey.toString();
				RecalDatum readGroupDatum = (RecalDatum) dataManager.getCollapsedTable(0).data.get(readGroupKey);
				logger.info(String.format("Writing out data tables for read group: %s\twith %s observations\tand aggregate residual error = %.3f", readGroup, readGroupDatum.getNumObservations(),
						readGroupDatum.empiricalQualDouble(0, MAX_QUALITY_SCORE) - readGroupDatum.getEstimatedQReported()));

				// for each covariate
				for (int iii = 1; iii < requestedCovariates.size(); iii++) {
					Covariate cov = requestedCovariates.get(iii);

					// Create a PrintStream
					File outputFile = new File(OUTPUT_DIR, readGroup + "." + cov.getClass().getSimpleName() + ".dat");
					PrintStream output;
					try {
						output = new PrintStream(FileUtils.openOutputStream(outputFile));
					} catch (IOException e) {
						throw new UserException.CouldNotCreateOutputFile(outputFile, e);
					}

					try {
						// Output the header
						output.println("Covariate\tQreported\tQempirical\tnMismatches\tnBases");

						for (Object covariateKey : ((Map) dataManager.getCollapsedTable(iii).data.get(readGroupKey)).keySet()) {
							output.print(covariateKey.toString() + "\t"); // Covariate
							RecalDatum thisDatum = (RecalDatum) ((Map) dataManager.getCollapsedTable(iii).data.get(readGroupKey)).get(covariateKey);
							output.print(String.format("%.3f", thisDatum.getEstimatedQReported()) + "\t"); // Qreported
							output.print(String.format("%.3f", thisDatum.empiricalQualDouble(0, MAX_QUALITY_SCORE)) + "\t"); // Qempirical
							output.print(thisDatum.getNumMismatches() + "\t"); // nMismatches
							output.println(thisDatum.getNumObservations()); // nBases
						}
					} finally {
						// Close the PrintStream
						IOUtils.closeQuietly(output);
					}
				}
			} else {
				break;
			}

		}
	}
}
