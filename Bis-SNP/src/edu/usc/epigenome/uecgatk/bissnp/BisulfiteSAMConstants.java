/**
 * 
 */
package edu.usc.epigenome.uecgatk.bissnp;

import java.io.Serializable;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Sep 15, 2013 11:28:15 PM
 * 
 */
public class BisulfiteSAMConstants implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 8517920216814454852L;
	
	public static final  String BIT_WISE_TAG = "Xn"; //use BIT_WISE_TAG tag to mark reads that have been thought duplicated, it will greatly increase the speed... but it cause some position not consistent in methy value...
	public static final  String MD_TAG = "MD";
	public static final  int REF_WINDOW_LEN = 500;
	public static final int MAX_COV = 999999999;
	public static final String CONT_CHRM_UCSC = "chrM";
	public static final String CONT_CHRM_NCBI = "MT";
}
