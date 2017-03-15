/**
 * 
 */
package edu.usc.epigenome.uecgatk.NOMeSeqWalker;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.math3.stat.inference.TTest;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteVCFConstants;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteVariantCallContext;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObject;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObjectWriterImp;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Feb 28, 2013 8:51:33 PM
 * 
 */
public class DmrInBisulfiteVcfWalker extends RodWalker<BisulfiteVariantCallContext, BisulfiteVcfSlidingWindow> {

	/**
     * The variant file(s) to evaluate.
     */
    @Input(fullName="case_file", shortName = "case", doc="Input case VCF file(s)", required=true)
    public RodBinding<VariantContext> cases;

    /**
     * The variant file(s) to compare against.
     */
    @Input(fullName="control_file", shortName = "control", doc="Input control VCF file(s)", required=true)
    public RodBinding<VariantContext> controls;

    
	@Output(fullName = "summary_dmr_file", shortName = "dmr", doc = "write result matrix file ", required = false)
	public String dmrFile = null;
	
	@Output(fullName = "sliding_matrix", shortName = "matrix", doc = "write sliding window matrix file ", required = false)
	public String matrixFile = null;
	
	 @Argument(fullName="pValue", shortName="pValue", doc="p value criteria as the differentiate methylation. default: 0.01", required=false)
	 public double pValue = 0.01;
	
	 @Argument(fullName="minCT", shortName="minCT", doc="Minimum number of CT reads for each cytosine to call methylation level. default: 2", required=false)
	 public int minCT = 2;
	 
	 @Argument(fullName="minNumC", shortName="minNumC", doc="Minimum number of cytosines in the window to call methylation level. default: 5", required=false)
	 public int minNumC = 10;
	 
	 @Argument(fullName="windowLen", shortName="windowLen", doc="Minimum number of cytosines in the window to call methylation level. default: 200", required=false)
	 public int windowLen = 200;
	 
	 @Argument(fullName="cytosinePattern", shortName="cytosinePattern", doc="cytosinePattern to look at inside the window. ", required=false)
	 public List<String> cytosinePattern = null;
	 
	 @Argument(fullName="mainCytosinePattern", shortName="mainCytosinePattern", doc=" main cytosinePattern to define the window. ", required=false)
	 public String mainCytosinePattern = "HCG";

	private BisulfiteVcfSlidingWindow caseWindow = null;
	private BisulfiteVcfSlidingWindow controlWindow = null;

	private bedObjectWriterImp dmrWriter = null;
	private bedObjectWriterImp matrixWriter = null;
	
	public void initialize() {
		if(cytosinePattern == null){
			cytosinePattern = new ArrayList<String>();
			cytosinePattern.add(mainCytosinePattern);
		}
		caseWindow = new BisulfiteVcfSlidingWindow(cytosinePattern, minCT, minNumC, windowLen, mainCytosinePattern);
		controlWindow = new BisulfiteVcfSlidingWindow(cytosinePattern, minCT, minNumC, windowLen, mainCytosinePattern);
		if(dmrFile != null){
			File fn1 = new File(dmrFile);
			dmrWriter = new bedObjectWriterImp(fn1);
		}
		if(matrixFile != null){
			File fn2 = new File(matrixFile);
			matrixWriter = new bedObjectWriterImp(fn2);
		}
		
        // maintain the full list of comps
    //    comps.addAll(compsProvided);
        
        
    }
	
	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.LocusWalker#map(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */
	@Override
	public BisulfiteVariantCallContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
		if(tracker==null)
			return null;
		List<VariantContext> case_bindings = tracker.getValues(cases);
		List<VariantContext> control_bindings = tracker.getValues(controls);
		int casePreAddLen = caseWindow.getWindowLen();
		int controlPreAddLen = controlWindow.getWindowLen();
		int casePostAddLen = 0;
		int controlPostAddLen = 0;
		boolean dataAddCase = false;
		boolean dataAddControl = false;
		if(!case_bindings.isEmpty() || !control_bindings.isEmpty()){
			if(!case_bindings.isEmpty()){
				VariantContext vc_case = case_bindings.get(0);
				if(!caseWindow.getWindow().isEmpty()){
					casePostAddLen = caseWindow.getWindow().peekFirst().ref.getLocus().distance(ref.getLocus());
				}
				
				BisulfiteVariantCallContext bvcCase = new BisulfiteVariantCallContext(null, vc_case, context, ref);
				dataAddCase = caseWindow.add(bvcCase, true);
				//System.err.println(bvcCase.getVariantContext().getGenotype(0).getAttributeAsString("CP", "."));
			}
			if(!control_bindings.isEmpty()){
				VariantContext vc_control = control_bindings.get(0);
				if(!controlWindow.getWindow().isEmpty()){
					controlPostAddLen = controlWindow.getWindow().peekFirst().ref.getLocus().distance(ref.getLocus());
				}
				
				BisulfiteVariantCallContext bvcControl = new BisulfiteVariantCallContext(null, vc_control, context, ref);
				dataAddControl = controlWindow.add(bvcControl, true);
				//System.err.println(bvcControl.getVariantContext());
			}
			//if(!controlWindow.getWindow().isEmpty() && controlWindow.getWindow().peekFirst().ref.getLocus().getStart() == 7002987){
			//	System.err.println(controlWindow.getWindow().peekFirst().ref.getLocus() + "\t" + controlWindow.getWindow().peekLast().ref.getLocus() + "\t" + controlWindow.getWindowLen() + "\t" + controlPreAddLen + "\t" + controlPostAddLen );
			//}
			
			if(casePreAddLen <= windowLen && casePostAddLen > windowLen && controlPreAddLen <= windowLen && controlPostAddLen > windowLen && (dataAddCase || dataAddControl)){
				List<Object> content = new ArrayList<Object>();
				boolean haveData = false;
				for(String pattern : cytosinePattern){
					if(caseWindow.getCytosineNum(pattern) >= minNumC && controlWindow.getCytosineNum(pattern) >= minNumC ){
						haveData = true;
						
						double methy_case = caseWindow.getCytosineMethyWeightMean(pattern);
						double methy_control = controlWindow.getCytosineMethyWeightMean(pattern);
						if(matrixWriter != null){
							//content.add(caseWindow.getWindow().peekFirst().ref.getLocus());
							//content.add(caseWindow.getWindow().peekLast().ref.getLocus());
							content.add(methy_case);
							//content.add(caseWindow.getCytosineMethyMean(pattern));
						//	content.add(caseWindow.getCytosineNum(pattern));
							
						//	content.add(caseWindow.getCytosineCReads(pattern));
						//	content.add(caseWindow.getCytosineCTReads(pattern));
							content.add(methy_control);
						//	content.add(controlWindow.getCytosineMethyMean(pattern));
						//	content.add(controlWindow.getCytosineNum(pattern));
							
						//	content.add(controlWindow.getCytosineCReads(pattern));
						//	content.add(controlWindow.getCytosineCTReads(pattern));
							double[] methy_vector_case = caseWindow.getMethyVector(pattern);
							double[] methy_vector_control = controlWindow.getMethyVector(pattern);
							content.add(methyChangesPvalue(methy_vector_case, methy_vector_control));
						}
						if(dmrWriter != null){
							
						}
						
					}
					
				}
				if(haveData){
					bedObject tmpLine = new bedObject(controlWindow.getWindow().peekFirst().ref.getLocus().getContig(), controlWindow.getWindow().peekFirst().ref.getLocus().getStart(),controlWindow.getWindow().peekLast().ref.getLocus().getStop(), content );
					matrixWriter.add(tmpLine);
				}
				
			}
		}
		

        return null;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduceInit()
	 */
	@Override
	public BisulfiteVcfSlidingWindow reduceInit() {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public BisulfiteVcfSlidingWindow reduce(BisulfiteVariantCallContext value, BisulfiteVcfSlidingWindow sum) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public void onTraversalDone(BisulfiteVcfSlidingWindow result) {
		if(dmrWriter != null){
			dmrWriter.close();
		}
		if(matrixWriter != null){
			matrixWriter.close();
		}
		
		logger.info("Finished!");
	}
	
	public double methyChangesPvalue(double[] methy_vector_case, double[] methy_vector_control){
		TTest ttest = new TTest();
		return ttest.tTest(methy_vector_case, methy_vector_control);
	}
	
	/*
	public methyStatusChanges methyChanges(double methy_case, double methy_control, double[] methy_vector_case, double[] methy_vector_control){
		TTest ttest = new TTest();
		double p = ttest.tTest(methy_vector_case, methy_vector_control);
		if(p < pValue){
			if(){
				
			}
			else{
				
			}
		}
		else{
			if(){
				
			}
			else{
				
			}
		}
	}
	*/
	public enum methyStatusChanges{
		UNKNOWN, //vague..
		MH, //methylation high in both case
		MR, //methylation resistant
		MP, //methylation prone
		ML, //methylation loss
	}

}
