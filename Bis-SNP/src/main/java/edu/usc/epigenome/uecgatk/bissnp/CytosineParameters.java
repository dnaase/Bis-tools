/**
 * 
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 20, 2012 11:07:32 AM
 * 
 */
public class CytosineParameters {

	public double cytosineMethylation;
	/**
	 * 
	 */
	// 1-based coordinate
	public int cytosinePosition;
	public char cytosineStrand;

	// CT heterozygous loci would appear in SNP.vcf file,
	// but not appear in CpG.vcf or C.vcf file, it will
	// not appear in summary statitics of heterozygous
	// cytosine pattern, since C/T heterozygous would
	// lead methylation level not estimatable..
	public boolean isCTHeterozygousLoci = false;
	public boolean isCytosinePattern = false;

	// implement heterozygous, need to do it in the future, useful for Allele
	// specific analysis.. mean position at cytosine are heterozygous,
	public boolean isHeterozygousCytosinePattern = false;
	// mean position at context position rather than cytosine position are heterozygous,
	public boolean isHeterozygousInContextPosition = false;

	// is this reference cytosine pattern in the reference genome?
	public boolean isReferenceCytosinePattern = false;
	public boolean methylationAutoestimated;
	public int numOfCReadsInCytosinePosInBisulfiteCStrandAlleleA = 0;
	public int numOfCReadsInCytosinePosInBisulfiteCStrandAlleleB = 0;
	public int numOfOtherReadsInCytosinePosInBisulfiteCStrandAlleleA = 0;
	public int numOfOtherReadsInCytosinePosInBisulfiteCStrandAlleleB = 0;
	public int numOfTReadsInCytosinePosInBisulfiteCStrandAlleleA = 0;
	public int numOfTReadsInCytosinePosInBisulfiteCStrandAlleleB = 0;

	public String patternOfAlleleA = null;

	public String patternOfAlleleB = null;

	public CytosineParameters() {

	}

}
