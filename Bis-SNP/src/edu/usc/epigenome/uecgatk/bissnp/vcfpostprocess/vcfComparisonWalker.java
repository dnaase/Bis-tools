/**
 * 
 */
package edu.usc.epigenome.uecgatk.bissnp.vcfpostprocess;

import java.util.Iterator;
import java.util.List;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time May 5, 2012 5:09:11 PM
 * 
 */
public class vcfComparisonWalker extends RodWalker<Integer, Integer> {

	/**
	 * 
	 */

	/**
	 * The variant file(s) to compare against.
	 */
	@Input(fullName = "comp", shortName = "comp", doc = "Input comparison file(s)", required = false)
	public RodBinding<VariantContext> comp;

	/**
	 * The variant file(s) to evaluate.
	 */
	@Input(fullName = "eval", shortName = "eval", doc = "Input evaluation file(s)", required = true)
	public RodBinding<VariantContext> eval;

	@Argument(fullName = "homCytosineComp", shortName = "homC", doc = "compare homozygous Cytosines, otherwise compare heterozygous SNP", required = false)
	protected boolean homC = false;

	@Argument(fullName = "minQualityEval", shortName = "mq", doc = "Minimum genotyping quality for evaluation vcf file", required = false)
	protected double mq = 20.0;

	@Argument(fullName = "minQualityStandard", shortName = "mqComp", doc = "Minimum genotyping quality for standard vcf file", required = false)
	protected double mqComp = 50.0;

	private int comp_num = 0;
	private int eval_num = 0;
	private int eval_num_filtered = 0;
	private int eval_num_no_call = 0;
	private int fn = 0;
	private int fp = 0;
	private int tn = 0;
	private int tp = 0;

	@Override
	public void initialize() {

	}

	public boolean isHetSnp(VariantContext vc) {
		if (vc.hasGenotypes()) {
			Iterator<Genotype> it = vc.getGenotypes().iterator();
			while (it.hasNext()) {
				if (it.next().isHet())
					return true;
			}
		}
		return false;
	}

	public boolean isHomoC(VariantContext vc) {
		if (vc.hasGenotypes()) {
			Iterator<Genotype> it = vc.getGenotypes().iterator();
			while (it.hasNext()) {
				Genotype tmp = it.next();
				// System.err.println(tmp.getGenotypeString());
				if (tmp.getGenotypeString().equalsIgnoreCase("C/C") || tmp.getGenotypeString().equalsIgnoreCase("G/G")) {
					return true;
				}
			}
		}
		return false;
	}

	@Override
	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

		if (tracker == null)
			return null;
		List<VariantContext> comp_bindings = tracker.getValues(comp);
		if (!comp_bindings.isEmpty()) {

			VariantContext vc_comp = comp_bindings.get(0);
			comparison(vc_comp, tracker);

			comp_num++;
		}

		return null;
	}

	@Override
	public void onTraversalDone(Integer result) {
		System.out.println("Evaluation loci: " + eval_num + "\tStandard loci: " + comp_num);
		System.out.println("Evaluation loci not called: " + eval_num_no_call + "\tEvaluation loci filtered out in false negative loci: " + eval_num_filtered);
		System.out.println("True Positive loci: " + tp);
		System.out.println("False Positive loci: " + fp);
		System.out.println("True Negative loci: " + tn);
		System.out.println("False Negative loci: " + fn);
		System.out.println("False Discoverary rate: " + (double) fp / (double) (fp + tp));
		System.out.println("Sensitivity rate: " + (double) tp / (double) (fn + tp));
	}

	@Override
	public Integer reduce(Integer value, Integer sum) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Integer reduceInit() {
		// TODO Auto-generated method stub
		return null;
	}

	private void comparison(VariantContext vc_comp, RefMetaDataTracker tracker) {
		List<VariantContext> eval_bindings = tracker.getValues(eval);
		if (!eval_bindings.isEmpty()) {
			VariantContext vc_eval = eval_bindings.get(0);

			eval_num++;
			if (homC) {
				if (isHomoC(vc_comp) && vc_comp.getPhredScaledQual() >= mq) {

					if (isHomoC(vc_eval) && vc_eval.getPhredScaledQual() >= mq) {
						tp++;
					} else {
						fn++;
						if (vc_eval.getPhredScaledQual() < mq) {
							eval_num_filtered++;
						}
					}

				} else {
					if (isHomoC(vc_eval) && vc_eval.getPhredScaledQual() >= mq) {
						fp++;
					} else {
						tn++;
					}
				}
			} else {
				if (isHetSnp(vc_comp) && vc_comp.getPhredScaledQual() >= mq) {
					if (isHetSnp(vc_eval) && vc_eval.getPhredScaledQual() >= mq) {
						tp++;
					} else {
						fn++;
						if (vc_eval.getPhredScaledQual() < mq) {
							eval_num_filtered++;
						}
					}

				} else {
					if (isHetSnp(vc_eval) && vc_eval.getPhredScaledQual() >= mq) {
						fp++;
					} else {
						tn++;
					}
				}
			}
		} else {
			eval_num_no_call++;
			if (homC) {
				if (isHomoC(vc_comp) && vc_comp.getPhredScaledQual() >= mq) {
					fn++;

				} else {
					tn++;
				}
			} else {
				if (isHetSnp(vc_comp) && vc_comp.getPhredScaledQual() >= mq) {
					fn++;
				} else {
					tn++;
				}
			}
		}
	}
}
