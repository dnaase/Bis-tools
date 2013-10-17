package edu.usc.epigenome.uecgatk.bissnp.bisulfiterecalibration;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Gather;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.recalibration.CountCovariatesGatherer;
import org.broadinstitute.sting.gatk.walkers.recalibration.Covariate;
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalDataManager;
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.recalibration.TableRecalibrationWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.collections.NestedHashMap;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import edu.usc.epigenome.uecgatk.bissnp.BaseUtilsMore;
import edu.usc.epigenome.uecgatk.bissnp.BisSNPUtils;

@BAQMode(ApplicationTime = BAQ.ApplicationTime.FORBIDDEN)
@By(DataSource.READS)
// Only look at covered loci, not every loci of the reference file
// @ReadFilters({MappingQualityZeroFilter.class,
// MappingQualityUnavailableFilter.class})
@ReadFilters({ MappingQualityZeroFilter.class, UnmappedReadFilter.class, BadMateFilter.class, NotPrimaryAlignmentFilter.class, DuplicateReadFilter.class })
// Filter out all reads with zero or unavailable mapping quality
@Requires({ DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES })
// This walker requires both -I input.bam and -R reference.fasta
@PartitionBy(PartitionType.LOCUS)
@Reference(window = @Window(start = -500, stop = 500))
public class BisulfiteCountCovariatesWalker extends LocusWalker<BisulfiteCountCovariatesWalker.CountedData, BisulfiteCountCovariatesWalker.CountedData> implements
		TreeReducible<BisulfiteCountCovariatesWalker.CountedData> {

	private static final String COVARS_ATTRIBUTE = "COVARS";
	private static int DBSNP_VALIDATION_CHECK_FREQUENCY = 1000000;
	private static final double DBSNP_VS_NOVEL_MISMATCH_RATE = 2.0;

	private static final String SEEN_ATTRIBUTE = "SEEN";

	// ///////////////////////////
	// Constants
	// ///////////////////////////
	private static final String SKIP_RECORD_ATTRIBUTE = "SKIP";

	// ///////////////////////////
	// Command Line Arguments
	// ///////////////////////////
	/**
	 * This algorithm treats every reference mismatch as an indication of error. However, real
	 * genetic variation is expected to mismatch the reference, so it is critical that a database of
	 * known polymorphic sites is given to the tool in order to skip over those sites. This tool
	 * accepts any number of RodBindings (VCF, Bed, etc.) for use as this database. For users
	 * wishing to exclude an interval list of known variation simply use -XL my.interval.list to
	 * skip over processing those sites. Please note however that the statistics reported by the
	 * tool will not accurately reflected those sites skipped by the -XL argument.
	 */
	@Input(fullName = "knownSites", shortName = "knownSites", doc = "A database of known polymorphic sites to skip over in the recalibration algorithm", required = false)
	public List<RodBinding<Feature>> knownSites = Collections.emptyList();

	/**
	 * After the header, data records occur one per line until the end of the file. The first
	 * several items on a line are the values of the individual covariates and will change depending
	 * on which covariates were specified at runtime. The last three items are the data- that is,
	 * number of observations for this combination of covariates, number of reference mismatches,
	 * and the raw empirical quality score calculated by phred-scaling the mismatch rate.
	 */
	@Output(fullName = "recal_file", shortName = "recalFile", required = true, doc = "Filename for the output covariates table recalibration file")
	@Gather(CountCovariatesGatherer.class)
	public PrintStream RECAL_FILE;

	// ///////////////////////////
	// Shared Arguments
	// ///////////////////////////
	@ArgumentCollection
	protected RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();
	/**
	 * See the -list argument to view available covariates.
	 */
	@Argument(fullName = "covariate", shortName = "cov", doc = "Covariates to be used in the recalibration. Each covariate is given as a separate cov parameter. ReadGroup and ReportedQuality are required covariates and are already added for you.", required = false)
	private String[] COVARIATES = null;

	// ///////////////////////////
	// Private Member Variables
	// ///////////////////////////
	private final BisulfiteRecalDataManager dataManager = new BisulfiteRecalDataManager();

	// ///////////////////////////
	// Debugging-only Arguments
	// ///////////////////////////
	@Argument(fullName = "dont_sort_output", shortName = "unsorted", required = false, doc = "If specified, the output table recalibration csv file will be in an unsorted, arbitrary order to save some run time.")
	private boolean DONT_SORT_OUTPUT = false;

	@Argument(fullName = "list", shortName = "ls", doc = "List the available covariates and exit", required = false)
	private boolean LIST_ONLY = false;
	private final ArrayList<Covariate> requestedCovariates = new ArrayList<Covariate>();
	/**
	 * This calculation is critically dependent on being able to skip over known polymorphic sites.
	 * Please be sure that you know what you are doing if you use this option.
	 */
	@Argument(fullName = "run_without_dbsnp_potentially_ruining_quality", shortName = "run_without_dbsnp_potentially_ruining_quality", required = false, doc = "If specified, allows the recalibrator to be used without a dbsnp rod. Very unsafe and for expert users only.")
	private boolean RUN_WITHOUT_DBSNP = false;
	@Argument(fullName = "standard_covs", shortName = "standard", doc = "Use the standard set of covariates in addition to the ones listed using the -cov argument", required = false)
	private boolean USE_STANDARD_COVARIATES = false;

	/**
	 * Update the mismatch / total_base counts for a given class of loci.
	 * 
	 * @param counter
	 *            The CountedData to be updated
	 * @param context
	 *            The AlignmentContext which holds the reads covered by this locus
	 * @param refBase
	 *            The reference base
	 */
	private static void updateMismatchCounts(CountedData counter, final AlignmentContext context, final byte refBase) {
		for (PileupElement p : context.getBasePileup()) {
			byte readBase = p.getBase();
			int readBaseIndex = BaseUtils.simpleBaseToBaseIndex(readBase);
			int refBaseIndex = BaseUtils.simpleBaseToBaseIndex(refBase);
			boolean negStrand = p.getRead().getReadNegativeStrandFlag();
			if (readBaseIndex != -1 && refBaseIndex != -1) {
				
				if (p.getRead().getReadPairedFlag()) {

					if (!BaseUtils.basesAreEqual(refBase, readBase) && BaseUtilsMore.isBisulfiteMismatch(refBase,readBase, negStrand,p.getRead().getSecondOfPairFlag())){
						counter.novelCountsMM++;
						//System.err.println(counter.novelCountsMM);
					}
					/*
					if (p.getRead().getSecondOfPairFlag()) {
						//readBase = BaseUtils.simpleComplement(readBase);
						//readBaseIndex = BaseUtils.simpleBaseToBaseIndex(readBase);
						negStrand = !negStrand;
						
						if (readBaseIndex != refBaseIndex) {
							if ((BisSNPUtils.isCytosine(refBase, false) && BisSNPUtils.isCytosine(readBase, true) )
									|| (BisSNPUtils.isCytosine(BaseUtils.simpleComplement(refBase), false) && BisSNPUtils.isCytosine(BaseUtils.simpleComplement(readBase), true))) {
								System.err.println((char)refBase + "\t" + (char)readBase + "\t" + negStrand);
							} else {
								//System.err.println((char)refBase + "\t" + (char)readBase + "\t" + negStrand);
								counter.novelCountsMM++;
							}
						}

					} else {
						if (readBaseIndex != refBaseIndex) {
							if ((BisSNPUtils.isCytosine(refBase, false) && BisSNPUtils.isCytosine(readBase, true) && !negStrand)
									|| (negStrand && BisSNPUtils.isCytosine(BaseUtils.simpleComplement(refBase), false) && BisSNPUtils.isCytosine(BaseUtils.simpleComplement(readBase), true))) {

							} else {
								counter.novelCountsMM++;
							}
						}

					}
					*/
				}

			} else {
				if (readBaseIndex != refBaseIndex) {
					if ((BisSNPUtils.isCytosine(refBase, false) && BisSNPUtils.isCytosine(readBase, true) && !negStrand)
							|| (negStrand && BisSNPUtils.isCytosine(BaseUtils.simpleComplement(refBase), false) && BisSNPUtils.isCytosine(BaseUtils.simpleComplement(readBase), true))) {

					} else {
						counter.novelCountsMM++;
					}
				}

			}

			counter.novelCountsBases++;
		}
	}

	// ---------------------------------------------------------------------------------------------------------------
	//
	// initialize
	//
	// ---------------------------------------------------------------------------------------------------------------

	/**
	 * Parse the -cov arguments and create a list of covariates to be used here Based on the
	 * covariates' estimates for initial capacity allocate the data hashmap
	 */
	@Override
	public void initialize() {

		if (RAC.FORCE_PLATFORM != null) {
			RAC.DEFAULT_PLATFORM = RAC.FORCE_PLATFORM;
		}

		// Get a list of all available covariates
		final List<Class<? extends Covariate>> covariateClasses = new PluginManager<Covariate>(Covariate.class).getPlugins();
		final List<Class<? extends RequiredCovariate>> requiredClasses = new PluginManager<RequiredCovariate>(RequiredCovariate.class).getPlugins();
		final List<Class<? extends StandardCovariate>> standardClasses = new PluginManager<StandardCovariate>(StandardCovariate.class).getPlugins();

		// Print and exit if that's what was requested
		if (LIST_ONLY) {
			logger.info("Available covariates:");
			for (Class<?> covClass : covariateClasses) {
				logger.info(covClass.getSimpleName());
			}
			logger.info("");

			System.exit(0); // Early exit here because user requested it
		}

		// Warn the user if no dbSNP file or other variant mask was specified
		if (knownSites.isEmpty() && !RUN_WITHOUT_DBSNP) {
			throw new UserException.CommandLineException(
					"This calculation is critically dependent on being able to skip over known variant sites. Please provide a VCF file containing known sites of genetic variation.");
		}

		if (USE_STANDARD_COVARIATES) {
			// We want the standard covariates to appear in a consistent order
			// but the packageUtils method gives a random order
			// A list of Classes can't be sorted, but a list of Class names can
			// be
			final List<String> standardClassNames = new ArrayList<String>();
			for (Class<?> covClass : standardClasses) {
				standardClassNames.add(covClass.getName());
			}
			Collections.sort(standardClassNames); // Sort the list of class
													// names
			for (String className : standardClassNames) {
				for (Class<?> covClass : standardClasses) { // Find the class
															// that matches this
															// class name
					if (covClass.getName().equals(className)) {
						try {
							final Covariate covariate = (Covariate) covClass.newInstance();
							requestedCovariates.add(covariate);
						} catch (Exception e) {
							throw new DynamicClassResolutionException(covClass, e);
						}
					}
				}
			}
		}
		// Finally parse the -cov arguments that were provided, skipping over
		// the ones already specified
		if (COVARIATES != null) {
			for (String requestedCovariateString : COVARIATES) {
				boolean foundClass = false;
				for (Class<?> covClass : covariateClasses) {
					if (requestedCovariateString.equalsIgnoreCase(covClass.getSimpleName())) { // -cov
																								// argument
																								// matches
																								// the
																								// class
																								// name
																								// for
																								// an
																								// implementing
																								// class
						foundClass = true;
						if (!requiredClasses.contains(covClass) && (!USE_STANDARD_COVARIATES || !standardClasses.contains(covClass))) {
							try {
								// Now that we've found a matching class, try to
								// instantiate it
								final Covariate covariate = (Covariate) covClass.newInstance();
								requestedCovariates.add(covariate);
							} catch (Exception e) {
								throw new DynamicClassResolutionException(covClass, e);
							}
						}
					}
				}

				if (!foundClass) {
					throw new UserException.CommandLineException("The requested covariate type (" + requestedCovariateString
							+ ") isn't a valid covariate option. Use --list to see possible covariates.");
				}
			}
		}

		logger.info("The covariates being used here: ");
		for (Covariate cov : requestedCovariates) {
			logger.info("\t" + cov.getClass().getSimpleName());
			cov.initialize(RAC); // Initialize any covariate member variables
									// using the shared argument collection
		}

	}

	// ---------------------------------------------------------------------------------------------------------------
	//
	// map
	//
	// ---------------------------------------------------------------------------------------------------------------

	/**
	 * For each read at this locus get the various covariate values and increment that location in
	 * the map based on whether or not the base matches the reference at this particular location
	 * 
	 * @param tracker
	 *            The reference metadata tracker
	 * @param ref
	 *            The reference context
	 * @param context
	 *            The alignment context
	 * @return Returns 1, but this value isn't used in the reduce step
	 */
	@Override
	public CountedData map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
		// Only use data from non-dbsnp sites
		// Assume every mismatch at a non-dbsnp site is indicative of poor
		// quality
		CountedData counter = new CountedData();
		// System.err.println("haha");
		if (tracker.getValues(knownSites).size() == 0) {
			// For each read at this locus
			for (final PileupElement p : context.getBasePileup()) {
				final GATKSAMRecord gatkRead = p.getRead();
				int offset = p.getOffset();

				if (gatkRead.containsTemporaryAttribute(SKIP_RECORD_ATTRIBUTE)) {
					continue;
				}

				if (!gatkRead.containsTemporaryAttribute(SEEN_ATTRIBUTE)) {
					gatkRead.setTemporaryAttribute(SEEN_ATTRIBUTE, true);
					RecalDataManager.parseSAMRecord(gatkRead, RAC);

					// Skip over reads with no calls in the color space if the
					// user requested it
					if (!(RAC.SOLID_NOCALL_STRATEGY == BisulfiteRecalDataManager.SOLID_NOCALL_STRATEGY.THROW_EXCEPTION) && RecalDataManager.checkNoCallColorSpace(gatkRead)) {
						gatkRead.setTemporaryAttribute(SKIP_RECORD_ATTRIBUTE, true);
						continue;
					}

					RecalDataManager.parseColorSpace(gatkRead);
					gatkRead.setTemporaryAttribute(COVARS_ATTRIBUTE, BisulfiteRecalDataManager.computeCovariates(gatkRead, requestedCovariates, ref));
				}

				// Skip this position if base quality is zero
				if (gatkRead.getBaseQualities()[offset] > 0) {

					byte[] bases = gatkRead.getReadBases();
					byte refBase = ref.getBase();

					// Skip if this base is an 'N' or etc.
					if (BaseUtils.isRegularBase(bases[offset])) {
						// This base finally passed all the checks for a good
						// base, so add it to the big data hashmap
						updateDataFromRead(counter, gatkRead, offset, refBase);
					}
				}
			}
			counter.countedSites++;
		} else { // We skipped over the dbSNP site, and we are only processing
					// every Nth locus
			counter.skippedSites++;
			
			updateMismatchCounts(counter, context, ref.getBase());
		}

		return counter;
	}

	/**
	 * Write out the full data hashmap to disk in CSV format
	 * 
	 * @param sum
	 *            The CountedData to write out to RECAL_FILE
	 */
	@Override
	public void onTraversalDone(CountedData sum) {
		logger.info("Writing raw recalibration data...");
		if (sum.countedBases == 0L) {
			throw new UserException.BadInput("Could not find any usable data in the input BAM file(s).");
		}
		outputToCSV(sum, RECAL_FILE);
		logger.info("...done!");
	}

	/**
	 * The Reduce method doesn't do anything for this walker.
	 * 
	 * @param mapped
	 *            Result of the map. This value is immediately ignored.
	 * @param sum
	 *            The summing CountedData used to output the CSV data
	 * @return returns The sum used to output the CSV data
	 */
	@Override
	public CountedData reduce(CountedData mapped, CountedData sum) {
		// Do a dbSNP sanity check every so often
		return validatingDbsnpMismatchRate(sum.add(mapped));
	}

	// ---------------------------------------------------------------------------------------------------------------
	//
	// reduce
	//
	// ---------------------------------------------------------------------------------------------------------------

	/**
	 * Initialize the reduce step by creating a PrintStream from the filename specified as an
	 * argument to the walker.
	 * 
	 * @return returns A PrintStream created from the -recalFile filename argument specified to the
	 *         walker
	 */
	@Override
	public CountedData reduceInit() {
		return new CountedData();
	}

	@Override
	public CountedData treeReduce(CountedData sum1, CountedData sum2) {
		return validatingDbsnpMismatchRate(sum1.add(sum2));
	}

	/**
	 * For each entry (key-value pair) in the data hashmap output the Covariate's values as well as
	 * the RecalDatum's data in CSV format
	 * 
	 * @param recalTableStream
	 *            The PrintStream to write out to
	 */
	private void outputToCSV(CountedData sum, final PrintStream recalTableStream) {
		recalTableStream.printf("# Counted Sites    %d%n", sum.countedSites);
		recalTableStream.printf("# Counted Bases    %d%n", sum.countedBases);
		recalTableStream.printf("# Skipped Sites    %d%n", sum.skippedSites);
		recalTableStream.printf("# Fraction Skipped 1 / %.0f bp%n", (double) sum.countedSites / sum.skippedSites);

		if (sum.solidInsertedReferenceBases != 0) {
			recalTableStream.printf("# Fraction SOLiD inserted reference 1 / %.0f bases%n", (double) sum.countedBases / sum.solidInsertedReferenceBases);
			recalTableStream.printf("# Fraction other color space inconsistencies 1 / %.0f bases%n", (double) sum.countedBases / sum.otherColorSpaceInconsistency);
		}

		// Output header saying which covariates were used and in what order
		for (Covariate cov : requestedCovariates) {
			recalTableStream.print(cov.getClass().getSimpleName().split("Covariate")[0] + ",");
		}
		recalTableStream.println("nObservations,nMismatches,Qempirical");

		if (DONT_SORT_OUTPUT) {
			printMappings(recalTableStream, 0, new Object[requestedCovariates.size()], dataManager.data.data);
		} else {
			printMappingsSorted(recalTableStream, 0, new Object[requestedCovariates.size()], dataManager.data.data);
		}

		// print out an EOF marker
		recalTableStream.println(TableRecalibrationWalker.EOF_MARKER);
	}

	private void printMappings(final PrintStream recalTableStream, final int curPos, final Object[] key, final Map data) {
		for (Object comp : data.keySet()) {
			key[curPos] = comp;
			final Object val = data.get(comp);
			if (val instanceof BisulfiteRecalDatumOptimized) {
				// For each Covariate in the key
				for (Object compToPrint : key) {
					// Output the Covariate's value
					recalTableStream.print(compToPrint + ",");
				}
				// Output the RecalDatum entry
				recalTableStream.println(((BisulfiteRecalDatumOptimized) val).outputToCSV());
			} else { // Another layer in the nested hash map
				printMappings(recalTableStream, curPos + 1, key, (Map) val);
			}
		}
	}

	private void printMappingsSorted(final PrintStream recalTableStream, final int curPos, final Object[] key, final Map data) {
		final ArrayList<Comparable> keyList = new ArrayList<Comparable>();
		for (Object comp : data.keySet()) {
			keyList.add((Comparable) comp);
		}

		Collections.sort(keyList);

		for (Comparable comp : keyList) {
			key[curPos] = comp;
			final Object val = data.get(comp);
			if (val instanceof BisulfiteRecalDatumOptimized) {
				// For each Covariate in the key
				for (Object compToPrint : key) {
					// Output the Covariate's value
					recalTableStream.print(compToPrint + ",");
				}
				// Output the RecalDatum entry
				recalTableStream.println(((BisulfiteRecalDatumOptimized) val).outputToCSV());
			} else { // Another layer in the nested hash map
				printMappingsSorted(recalTableStream, curPos + 1, key, (Map) val);
			}
		}
	}

	/**
	 * Major workhorse routine for this walker. Loop through the list of requested covariates and
	 * pick out the value from the read, offset, and reference Using the list of covariate values as
	 * a key, pick out the RecalDatum and increment, adding one to the number of observations and
	 * potentially one to the number of mismatches Lots of things are passed as parameters to this
	 * method as a strategy for optimizing the covariate.getValue calls because pulling things out
	 * of the SAMRecord is an expensive operation.
	 * 
	 * @param counter
	 *            Data structure which holds the counted bases
	 * @param gatkRead
	 *            The SAMRecord holding all the data for this read
	 * @param offset
	 *            The offset in the read for this locus
	 * @param refBase
	 *            The reference base at this locus
	 */
	private void updateDataFromRead(CountedData counter, final GATKSAMRecord gatkRead, final int offset, final byte refBase) {
		final Object[][] covars = (Comparable[][]) gatkRead.getTemporaryAttribute(COVARS_ATTRIBUTE);
		final Object[] key = covars[offset];

		// Using the list of covariate values as a key, pick out the RecalDatum
		// from the data HashMap
		final NestedHashMap data = dataManager.data; // optimization - create
														// local reference
		BisulfiteRecalDatumOptimized datum = (BisulfiteRecalDatumOptimized) data.get(key);
		if (datum == null) { // key doesn't exist yet in the map so make a new
								// bucket and add it
			// initialized with zeros, will be incremented at end of method
			datum = (BisulfiteRecalDatumOptimized) data.put(new BisulfiteRecalDatumOptimized(), true, key);
		}

		// Need the bases to determine whether or not we have a mismatch
		final byte base = gatkRead.getReadBases()[offset];
		final long curMismatches = datum.getNumMismatches();

		// Add one to the number of observations and potentially one to the
		// number of mismatches
		datum.incrementBaseCountsBisulfite(base, refBase,gatkRead.getReadNegativeStrandFlag(),(gatkRead.getReadPairedFlag() && gatkRead.getSecondOfPairFlag()));
	//	if (gatkRead.getReadNegativeStrandFlag()) {
	//		datum.incrementBaseCountsBisulfite(BaseUtils.simpleComplement(base), BaseUtils.simpleComplement(refBase));
	//	} else {
	//		datum.incrementBaseCountsBisulfite(base, refBase);
	//	}

		counter.countedBases++;
		counter.novelCountsBases++;
		counter.novelCountsMM += datum.getNumMismatches() - curMismatches;
	}

	/**
	 * Validate the dbSNP reference mismatch rates.
	 */
	private CountedData validatingDbsnpMismatchRate(CountedData counter) {
		if (++counter.lociSinceLastDbsnpCheck >= DBSNP_VALIDATION_CHECK_FREQUENCY) {
			counter.lociSinceLastDbsnpCheck = 0;

			if (counter.novelCountsBases != 0L && counter.dbSNPCountsBases != 0L) {
				final double fractionMM_novel = (double) counter.novelCountsMM / (double) counter.novelCountsBases;
				final double fractionMM_dbsnp = (double) counter.dbSNPCountsMM / (double) counter.dbSNPCountsBases;

				if (fractionMM_dbsnp < DBSNP_VS_NOVEL_MISMATCH_RATE * fractionMM_novel) {
					Utils.warnUser("The variation rate at the supplied list of known variant sites seems suspiciously low. Please double-check that the correct ROD is being used. "
							+ String.format("[dbSNP variation rate = %.4f, novel variation rate = %.4f]", fractionMM_dbsnp, fractionMM_novel));
					DBSNP_VALIDATION_CHECK_FREQUENCY *= 2;
				}
			}
		}

		return counter;
	}

	public static class CountedData {
		public long countedBases = 0;
		public long countedSites = 0;
		public long dbSNPCountsMM = 0, dbSNPCountsBases = 0;
		public int lociSinceLastDbsnpCheck = 0;
		public long novelCountsMM = 0, novelCountsBases = 0;

		public long otherColorSpaceInconsistency = 0;
		public long skippedSites = 0;
		public long solidInsertedReferenceBases = 0;

		/**
		 * Adds the values of other to this, returning this
		 * 
		 * @param other
		 * @return this object
		 */
		public CountedData add(CountedData other) {
			countedSites += other.countedSites;
			countedBases += other.countedBases;
			skippedSites += other.skippedSites;
			solidInsertedReferenceBases += other.solidInsertedReferenceBases;
			otherColorSpaceInconsistency += other.otherColorSpaceInconsistency;
			dbSNPCountsMM += other.dbSNPCountsMM;
			dbSNPCountsBases += other.dbSNPCountsBases;
			novelCountsMM += other.novelCountsMM;
			novelCountsBases += other.novelCountsBases;
			lociSinceLastDbsnpCheck += other.lociSinceLastDbsnpCheck;
			return this;
		}
	}
}
