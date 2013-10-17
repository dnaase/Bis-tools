/**
 * 
 */
package edu.usc.epigenome.dmntools.utils;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Oct 16, 2013 9:57:56 PM
 * 
 */
public class BenjaminiHochbergFDR {
	/*--------------------------------------------------------------
    FIELDS.
    --------------------------------------------------------------*/

    private HashEntry[] hash ;
    /**
     * the GO labels that have been tested (constructor input).
     */
    private static String [] goLabels;
    /**
     * the raw p-values that were given as input for the constructor, order corresponds to String [] goLabels.
     */
    private static String [] pvalues;
    /**
     * the goLabels ordened according to the ordened pvalues.
     */
    private static String [] ordenedGOLabels;
    /**
     * the raw p-values ordened in ascending order.
     */
    private static String [] ordenedPvalues;
    /**
     * the adjusted p-values ordened in ascending order.
     */
    private static String [] adjustedPvalues;

    /**
     * hashmap with the results (adjusted p-values) as values and the GO labels as keys.
     */
    private static HashMap correctionMap;

    /**
     * the significance level.
     */
    private static BigDecimal alpha;
    /**
     * the number of tests.
     */
    private static int m;
    /**
     * scale for the division in de method 'runFDR'.
     */
    private static final int RESULT_SCALE = 100;

    // Keep track of progress for monitoring:

    private int maxValue;
    private TaskMonitor taskMonitor = null;
    private boolean interrupted = false;
   

    /*--------------------------------------------------------------
    CONSTRUCTOR.
    --------------------------------------------------------------*/

    /**
     * Constructor.
     *
     * @param golabelstopvalues Hashmap of Strings with the goLabels and their pvalues.
     * @param alpha             String with the desired significance level.
     */

    public BenjaminiHochbergFDR(double[] pvalue) {

        this.adjustedPvalues = new String[m];

    }
 
    
    /*--------------------------------------------------------------
    * METHODS. adapted from cytoscape software: 
    * https://code.google.com/p/clusterviz-cytoscape/source/browse/trunk/src/clusterviz/BiNGO/BiNGO/BenjaminiHochbergFDR.java?r=33
    --------------------------------------------------------------*/

    /**
     * method that calculates the Benjamini and Hochberg correction of
     * the false discovery rate
     * NOTE : convert array indexes [0..m-1] to ranks [1..m].
     * orden raw p-values low .. high
     * test p<(i/m)*alpha from high to low (for i=m..1)
     * i* (istar) first i such that the inequality is correct.
     * reject hypothesis for i=1..i* : labels 1..i* are overrepresented
     * <p/>
     * adjusted p-value for i-th ranked p-value p_i^adj = min(k=i..m)[min(1,m/k p_k)]
     */

    public void calculate() {

        // ordening the pvalues.

        java.util.Arrays.sort(hash, new HashComparator()) ; 
        this.ordenedPvalues = parse(hash);
        // calculating adjusted p-values.
        BigDecimal min = new BigDecimal("" + 1);
        BigDecimal mkprk;
        for (int i = m; i > 0; i--) {
            mkprk = (new BigDecimal("" + m).multiply(new BigDecimal(ordenedPvalues[i - 1]))).divide(new BigDecimal("" + i), RESULT_SCALE, BigDecimal.ROUND_HALF_UP);
            if (mkprk.compareTo(min) < 0) {
                min = mkprk;
            }
            adjustedPvalues[i - 1] = min.toString();

        }
        correctionMap = new HashMap();
        for (int i = 0; i < adjustedPvalues.length && i < ordenedGOLabels.length; i++) {
            correctionMap.put(ordenedGOLabels[i], adjustedPvalues[i]);
        }
    }

    
    public String [] parse(HashEntry [] data) {
        String[] keys = new String[data.length];
        String[] values = new String[data.length];
        for(int i = 0; i < data.length; i++){
            keys[i] = data[i].key;
            values[i] = data[i].value;
        }
        ordenedGOLabels = keys;
        return values;
    }

    /*--------------------------------------------------------------
      GETTERS.
    --------------------------------------------------------------*/

    /**
     * getter for the map of corrected p-values.
     *
     * @return HashMap correctionMap.
     */
    public HashMap getCorrectionMap() {
        return correctionMap;
    }

    /**
     * getter for the ordened p-values.
     *
     * @return String[] with the ordened p-values.
     */
    public String[] getOrdenedPvalues() {
        return ordenedPvalues;
    }

    /**
     * getter for the adjusted p-values.
     *
     * @return String[] with the adjusted p-values.
     */
    public String[] getAdjustedPvalues() {
        return adjustedPvalues;
    }

    /**
     * getter for the ordened GOLabels.
     *
     * @return String[] with the ordened GOLabels.
     */
    public String[] getOrdenedGOLabels() {
        return ordenedGOLabels;
    }


    /**
     * Run the Task.
     */
    public void run() {
        calculate();
    }

    /**
     * Non-blocking call to interrupt the task.
     */
    public void halt() {
        this.interrupted = true;
    }

    /**
     * Sets the Task Monitor.
     *
     * @param taskMonitor TaskMonitor Object.
     */
    public void setTaskMonitor(TaskMonitor taskMonitor) {
        if (this.taskMonitor != null) {
            throw new IllegalStateException("Task Monitor is already set.");
        }
        this.taskMonitor = taskMonitor;
    }

    /**
     * Gets the Task Title.
     *
     * @return human readable task title.
     */
    public String getTitle() {
        return new String("Calculating FDR correction");
    }


}