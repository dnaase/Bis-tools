/**
 * maybe it is better to write rodWlker rather than Locus walker for this kind of alignment..
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp.otherwalker;

import java.io.File;
import java.util.LinkedList;
import java.util.List;

import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.bed.SimpleBEDFeature;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;

import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.datasources.rmd.ReferenceOrderedDataSource;

import org.broadinstitute.gatk.engine.filters.DuplicateReadFilter;
import org.broadinstitute.gatk.engine.filters.MappingQualityFilter;
import org.broadinstitute.gatk.engine.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.gatk.engine.filters.UnmappedReadFilter;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import org.broadinstitute.gatk.utils.refdata.utils.LocationAwareSeekableRODIterator;

import org.broadinstitute.gatk.utils.refdata.utils.RODRecordList;

import org.broadinstitute.gatk.engine.walkers.By;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.Downsample;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.engine.walkers.ReadFilters;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.engine.walkers.Window;

import org.broadinstitute.gatk.utils.GenomeLoc;


import main.java.edu.usc.epigenome.uecgatk.bissnp.BisulfiteArgumentCollection;
import main.java.edu.usc.epigenome.uecgatk.bissnp.BisulfiteGenotyperEngine;
import main.java.edu.usc.epigenome.uecgatk.bissnp.BisulfiteVariantCallContext;

import org.broadinstitute.gatk.engine.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.gatk.engine.filters.MappingQualityZeroFilter;
import main.java.edu.usc.epigenome.uecgatk.bissnp.filters.*;
import main.java.edu.usc.epigenome.uecgatk.bissnp.writer.bedObject;
import main.java.edu.usc.epigenome.uecgatk.bissnp.writer.bedObjectWriterImp;
import java.util.HashSet;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 */

// require input bam file, reference sequence, dbSNP rod, -L interval list bed file, -B feature list bed file
@ReadFilters( {UnmappedReadFilter.class, DuplicateReadFilter.class, NotPrimaryAlignmentFilter.class, FailsVendorQualityCheckFilter.class, MappingQualityFilter.class, MappingQualityZeroFilter.class, NotProperPairedReadFilter.class, InvertedDupsReadFilter.class, BisulfiteMismatchReadsFilter.class, BisulfiteIncompleteConvReadsFilter.class, BisulfiteFivePrimeConvReadsFilter.class} ) // Filter out all reads with zero mapping quality
@Reference(window=@Window(start=-500,stop=500))
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.NONE)
public class MethyPatternFeatureByBamWalker extends LocusWalker<Boolean, Boolean>
		implements TreeReducible<Boolean> {

	@ArgumentCollection private static BisulfiteArgumentCollection BAC = new BisulfiteArgumentCollection();
	
	@Input(fullName="aligned_feature", shortName = "feature" , doc="Input feature location to align", required=false)
	public RodBinding<BEDFeature> feature;
	
	@Output(fullName = "methy_file_name", shortName = "methyFile", doc = "Output file name with methylation pattern around the feature", required = false)
    public String methyFile = null;
	
	
	@Output(fullName = "datapoint_file_name", shortName = "dataPoint", doc = "Output file name with Data point/Tag density pattern around the feature", required = false)
    public String dataPoint = null;
	
//	@Argument(fullName = "feature_name", shortName = "feature", doc = "Feature name provide in -B:<name>,<type> <filename> option", required = false)
 //   public String feature = null;
	
	@Argument(fullName = "search_distance_to_feature", shortName = "distance", doc = "define the distance before or after feature", required = false)
    public int distance = 2000;
	
	@Argument(fullName = "bin_size", shortName = "binSize", doc = "define the bin size when sliding window. default: 1", required = false)
    public int binSize = 1;
	
	@Argument(fullName = "minium_CT_reads_count", shortName = "minCTdepth", doc = "minium number of CT reads should contained to calculate methylation value", required = false)
    public int minCTdepth = 1;
	
	@Argument(fullName = "enable_orientation", shortName = "orientated", doc = "orientated by strand or not", required = false)
    public boolean orientated = false;
	
	@Argument(fullName = "motif_alignment_type", shortName = "alignmentType", doc = "motif aligned at FiveEnd, ThreeEnd or Center", required = false)
    public MotifAlignmentType alignmentType = MotifAlignmentType.Center;
	
	private bedObjectWriterImp methyWriter = null;
	
	private bedObjectWriterImp dataPointWriter = null;
	
	private boolean inFeature = false;
	private boolean writtenObject = false;
	
	private LinkedList<Double> tmpMethyValueList = null;
	
	private LinkedList<Integer> tmpDataPointList = null;
	
	private ReferenceOrderedDataSource rodIt = null;
	
	private BisulfiteGenotyperEngine BG_engine = null;
	
	private String chr = null;
	private int bedStart = 0;
	private int bedEnd = 0;
	private SimpleBEDFeature bed = null;
	private Strand strand = Strand.NONE;
	private FeatureCondition summary= null;
	
	public void initialize(){
		 File fn1 = new File(methyFile);
		 File fn2= null;
		 if(dataPoint != null)
			fn2 = new File(dataPoint);
		 methyWriter = new bedObjectWriterImp(fn1);

		 if(dataPoint != null)
			 dataPointWriter = new bedObjectWriterImp(fn2);
		 
		 tmpMethyValueList = new LinkedList<Double>();

		 tmpDataPointList = new LinkedList<Integer>();
		 rodIt = getToolkit().getRodDataSources().get(1);
		 if(BAC.cytosineContextsAcquired.isEmpty()){
			 BAC.cytosineContextsAcquired.add("CG,1");	
		}
		 BAC.makeCytosine();
		 summary = new FeatureCondition();
	}
	   
	@Override
	public Boolean map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {

	     GenomeLoc loc = ref.getLocus().getLocation();
    	 GenomeLoc searchLoc = getToolkit().getGenomeLocParser().createGenomeLoc(ref.getLocus().getLocation().getContig(), ref.getLocus().getLocation().getStart()-distance-10, ref.getLocus().getLocation().getStart()+distance+10);
    	 
    	 LocationAwareSeekableRODIterator locRodIt = rodIt.seek(searchLoc);
    	 if(!locRodIt.hasNext()){
    		 inFeature = false; 
    		// locRodIt.close();
    		 rodIt.close(locRodIt);
    		
    		 return null;
    	 }
    	 else{
    		   
    		 if(!inFeature){
    			 RODRecordList rodList = locRodIt.seekForward(searchLoc);
    			// System.err.println(rodList.get(0).getUnderlyingObject());
    			 //rodList.get(0).getUnderlyingObject()
    			// for(Object it : rodList){
    				 if(rodList.get(0).getUnderlyingObject() instanceof SimpleBEDFeature){
    					 SimpleBEDFeature bedTmp = (SimpleBEDFeature)rodList.get(0).getUnderlyingObject();
    					// System.err.println(rodList.get(0).getUnderlyingObject());
    					 int featureAlignStart = (bedTmp.getStart() + bedTmp.getEnd())/2;
    					 
    					 if(alignmentType == MotifAlignmentType.FiveEnd){
    						 featureAlignStart = bedTmp.getStart();
    						 if(bedTmp.getStrand() == Strand.NEGATIVE){
    							 featureAlignStart = bedTmp.getEnd();
    						 }
    					 }
    					 else if(alignmentType == MotifAlignmentType.ThreeEnd){
    						 featureAlignStart = bedTmp.getEnd();
    						 if(bedTmp.getStrand() == Strand.NEGATIVE){
    							 featureAlignStart = bedTmp.getStart();
    						 }
    					 }
    					 
        				 if(loc.distance(getToolkit().getGenomeLocParser().createGenomeLoc(bedTmp.getChr(), featureAlignStart, featureAlignStart)) <= distance){
                			 bed = bedTmp;
                			 strand = bedTmp.getStrand();
            	    		 chr = bedTmp.getChr();
            	    		 bedStart = bedTmp.getStart();
            	    		 bedEnd = bedTmp.getEnd();
            	    		 inFeature = true;
            	    		 writtenObject = false;
            	    		// System.err.println(bed.getStart());
            	    		// break;
            	    	 }
        			// }
        			 
        			 }
        			 
    			 //} 
    		}
    		else{
     			// System.err.println(loc.distance(getToolkit().getGenomeLocParser().createGenomeLoc(bed.getChr(), (bed.getStart() + bed.getEnd())/2, (bed.getStart() + bed.getEnd())/2)) + "\t" + bed.toString());
    			int featureAlignStart = (bed.getStart() + bed.getEnd())/2;
				 
				 if(alignmentType == MotifAlignmentType.FiveEnd){
					 featureAlignStart = bed.getStart();
					 if(bed.getStrand() == Strand.NEGATIVE){
						 featureAlignStart = bed.getEnd();
					 }
				 }
				 else if(alignmentType == MotifAlignmentType.ThreeEnd){
					 featureAlignStart = bed.getEnd();
					 if(bed.getStrand() == Strand.NEGATIVE){
						 featureAlignStart = bed.getStart();
					 }
				 }	 
				 
    			if(bed != null && loc.distance(getToolkit().getGenomeLocParser().createGenomeLoc(bed.getChr(), featureAlignStart, featureAlignStart)) > distance){

     	    		 inFeature = false;
     	    		 writtenObject = false;
     	    		// System.err.println(bed.getStart());
     	    		// break;
     			}
     		}
    			 
            
    	 }
    	 
    	
    	 rodIt.close(locRodIt);
    	
	     if(!inFeature){
	    	 if(!writtenObject && !tmpMethyValueList.isEmpty()){
	    		// System.err.println(tmpMethyValueListGch.size() + "\t" + tmpMethyValueListWcg.size());
	    		 if(tmpMethyValueList.size() == distance * 2 + 1){
	    			 
	    			 bedObject bedLineGch = new bedObject(chr, bedStart, bedEnd, strand, (List)tmpMethyValueList); //chr, start, end, strand, aveMethyNDR, gchNumNDR, gchDepthNDR, gchCTdepthNDR, aveMethyLinker, gchNumLinker, gchDepthLinker, gchCTdepthLinker
		    		 methyWriter.add(bedLineGch);
		    		 if(dataPoint != null){
		    			 bedObject bedLineGch2 = new bedObject(chr, bedStart, bedEnd, strand, (List)tmpDataPointList);
			    		 dataPointWriter.add(bedLineGch2);
		    		 }
		    		 
	    		 }
	    		
	    		

	    		 tmpMethyValueList.clear();

	    		 if(dataPoint != null)
	    			 tmpDataPointList.clear();
	    		 
	    		 writtenObject = true;
	    		 logger.info(chr + "\t" + bedStart + "\t" + bedEnd);
	    	 }
	    	
	    	 return null;
	     }
	     else{

	 		if(context == null){
	 			
	 				addContextToList(null, strand, tmpMethyValueList);

	 				if(dataPoint != null)
	 					addCoverageToList(null, strand, tmpDataPointList);
	 			
	 				
	 			return null;
	 		}
	 		boolean isCytosinePattern =false;

	 		BG_engine = new BisulfiteGenotyperEngine(tracker, ref, context, BAC, getToolkit());
	 		BisulfiteVariantCallContext bvc = BG_engine.getBisulfiteVariantCallContext();
	 		if(bvc == null || bvc.getSummaryAcrossRG().cytosinePatternConfirmedSet ==null){
	 			addContextToList(null, strand, tmpMethyValueList);

 				if(dataPoint != null)
 					addCoverageToList(null, strand, tmpDataPointList);
 				return null;
	 			
	 		}
	 		//	return null;
	 		HashSet<String> cytosinePatternConfirmedList = bvc.getSummaryAcrossRG().cytosinePatternConfirmedSet;
	 		
	 		
	 			if(!cytosinePatternConfirmedList.isEmpty()){
	 				isCytosinePattern=true;
	 			}

	 				
	 		if(isCytosinePattern){
	 			
	 			addContextToList(bvc, strand, tmpMethyValueList);
	 			if(dataPoint != null)
	 				addCoverageToList(bvc, strand, tmpDataPointList);
	 			//if positive strand, offerLast(), if negative strand, offerFirst()
	 			
	 		}
	 		else{
	 			addContextToList(null, strand, tmpMethyValueList);
	 			if(dataPoint != null)
	 				addCoverageToList(null, strand, tmpDataPointList);
	 		}
	 		
	 		
	 		
	 		return null;
	     }

		
	}

	@Override
	public Boolean reduceInit() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Boolean reduce(Boolean value, Boolean sum) {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public Boolean treeReduce(Boolean lhs, Boolean rhs) {
		return null;
	}

	public void onTraversalDone(Boolean result) {
		
		methyWriter.close();

		if(dataPoint != null)
			dataPointWriter.close();
		logger.info("Finished!");
	}
	
	private void addContextToList(BisulfiteVariantCallContext bvc, Strand strand, LinkedList<Double> list){
		if(orientated){
			if(strand == Strand.NEGATIVE){
				if(bvc == null){
					list.offerFirst(Double.NaN);
				}
				else{
					if(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT >= minCTdepth){
						list.offerFirst(100*(double)bvc.getSummaryAcrossRG().numC/(double)(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT));
					}
					else{
						list.offerFirst(Double.NaN);
					}
					
				}
			}
			else{
				if(bvc == null){
					list.offerLast(Double.NaN);
				}
				else{
					if(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT >= minCTdepth){
						list.offerLast(100*(double)bvc.getSummaryAcrossRG().numC/(double)(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT));
					}
					else{
						list.offerFirst(Double.NaN);
					}
					
				}
			}
			
		}
		else{
			if(bvc == null){
				list.offerLast(Double.NaN);
			}
			else{
				if(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT >= minCTdepth){
					list.offerLast(100*(double)bvc.getSummaryAcrossRG().numC/(double)(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT));
				}
				else{
					list.offerFirst(Double.NaN);
				}
				
			}
		}
	}
	
	private void addCoverageToList(BisulfiteVariantCallContext bvc, Strand strand, LinkedList<Integer> list){
		if(orientated){
			if(strand == Strand.NEGATIVE){
				if(bvc == null){
					list.offerFirst(0);
				}
				else{
					list.offerFirst(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT);
				}
			}
			else{
				if(bvc == null){
					list.offerLast(0);
				}
				else{
					list.offerLast(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT);
				}
			}
			
		}
		else{
			if(bvc == null){
				list.offerLast(0);
			}
			else{
				list.offerLast(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT);
			}
		}
	}
/*	
	private List<Object> smoothList(LinkedList<Object> list, int binSize){
		if(binSize == 1){
			return list;
		}
		else{
			
			LinkedList<Object> returnedList = new LinkedList<Object>();
			
			Iterator<Object> it = list.iterator();
			double sum = 0;
			int num = 0;
			while(it.hasNext()){
				num++;
				if(num )
				sum += (Double)it.next();
			}
		}
		return list;
		
	}
*/	
	private boolean listStatistics(LinkedList<Double> list, String type){
		boolean filter= true;
		int numNaElement = 0;
		if(type.equalsIgnoreCase("GCH")){
			
		}
		else if(type.equalsIgnoreCase("HCG")){
			
		}
		else if(type.equalsIgnoreCase("WCG")){
			
		}
		
		return filter;
	}
	
	private enum MotifAlignmentType{
		FiveEnd,
		ThreeEnd,
		Center
	}

	private class FeatureCondition{
		public long visitedFeatures = 0;
		public long nNoHcgValueFeatures = 0; //number of feature with not even one value of HCG (all are NA..)
		public long nNoWcgValueFeatures = 0; //number of feature with not even one value of WCG
		public long nNoGchValueFeatures = 0; //number of feature with not even one value of GCH
		
		public long nHcgConfidantlyCalled = 0;
		public long nWcgConfidantlyCalled = 0;
		public long nGchConfidantlyCalled = 0;
		
		public long sumReadsHcgConfidantlyCalled = 0;
		public long sumReadsWcgConfidantlyCalled = 0;
		public long sumReadsGchConfidantlyCalled = 0;
		
		public long sumCTReadsHcgConfidantlyCalled = 0;
		public long sumCTReadsWcgConfidantlyCalled = 0;
		public long sumCTReadsGchConfidantlyCalled = 0;
		
		public double sumMethyReadsHcgConfidantlyCalled = 0;
		public double sumMethyReadsWcgConfidantlyCalled = 0;
		public double sumMethyReadsGchConfidantlyCalled = 0;
		
		public double averageOfReadsCoveredHcg(){return (double)sumReadsHcgConfidantlyCalled/(double)nHcgConfidantlyCalled;}
		public double averageOfReadsCoveredWcg(){return (double)sumReadsWcgConfidantlyCalled/(double)nWcgConfidantlyCalled;}
		public double averageOfReadsCoveredGch(){return (double)sumReadsGchConfidantlyCalled/(double)nGchConfidantlyCalled;}
		
		public double averageOfCTReadsCoveredHcg(){return (double)sumCTReadsHcgConfidantlyCalled/(double)nHcgConfidantlyCalled;}
		public double averageOfCTReadsCoveredWcg(){return (double)sumCTReadsWcgConfidantlyCalled/(double)nWcgConfidantlyCalled;}
		public double averageOfCTReadsCoveredGch(){return (double)sumCTReadsGchConfidantlyCalled/(double)nGchConfidantlyCalled;}
		
		public double averageOfMethyCoveredHcg(){return sumMethyReadsHcgConfidantlyCalled/(double)nHcgConfidantlyCalled;}
		public double averageOfMethyCoveredWcg(){return sumMethyReadsWcgConfidantlyCalled/(double)nWcgConfidantlyCalled;}
		public double averageOfMethyCoveredGch(){return sumMethyReadsGchConfidantlyCalled/(double)nGchConfidantlyCalled;}
		
	}
	

	
}
