/**
 * 
 */
package edu.usc.epigenome.uecgatk.bissnp;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 16, 2012 3:48:16 PM
 * 
 */
public class GATKSAMRecordFilterStorage {

	private BisulfiteArgumentCollection BAC;

	private GATKSAMRecord GATKrecord = null;
	private boolean goodBase = false;

	public GATKSAMRecordFilterStorage(GATKSAMRecord GATKrecord, BisulfiteArgumentCollection BAC, ReferenceContext refContext, int offset) {

		this.GATKrecord = GATKrecord;
		this.BAC = BAC;
		//setGoodBases(offset, refContext);

	}

	public GATKSAMRecordFilterStorage(GATKSAMRecord GATKrecord, BisulfiteArgumentCollection BAC, ReferenceContext refContext, int offset, boolean filterByOriginalQual) {
		// TODO Auto-generated constructor stub
		this.GATKrecord = GATKrecord;
		this.BAC = BAC;
		//setGoodBases(offset, refContext, filterByOriginalQual);

	}

	public boolean isGoodBase() {

		return goodBase;
	}
	/*
	private void setGoodBases(int offset, ReferenceContext refContext) {
		byte[] quals = GATKrecord.getBaseQualities();

		if (GATKrecord.getMappingQuality() >= BAC.MIN_MAPPING_QUALTY_SCORE && quals[offset] >= BAC.MIN_BASE_QUALTY_SCORE
				&& (BAC.USE_BADLY_MATED_READS || (!BadMateFilter.hasBadMate(GATKrecord)) && !GATKrecord.getNotPrimaryAlignmentFlag())) {
			if(BAC.trim2nd){
				if(GATKrecord.getReadPairedFlag() && GATKrecord.getSecondOfPairFlag()){
					if(GATKrecord.getReadNegativeStrandFlag()){
						if(offset >= BAC.trim3 && offset <= GATKrecord.getReadLength() - BAC.trim5){
							goodBase = true;
						}
					}
					else{
						if(offset >= BAC.trim5 && offset <= GATKrecord.getReadLength() - BAC.trim3){
							goodBase = true;
						}
					}
				}
				else{
					goodBase = true;
				}
			}
			else{
				if(GATKrecord.getReadNegativeStrandFlag()){
					if(offset >= BAC.trim3 && offset <= GATKrecord.getReadLength() - BAC.trim5){
						goodBase = true;
					}
				}
				else{
					if(offset >= BAC.trim5 && offset <= GATKrecord.getReadLength() - BAC.trim3){
						goodBase = true;
					}
				}
			}
			
			//System.err.println(offset + "\t" + GATKrecord.getReadLength() + "\t" + GATKrecord.getReadBases()[offset] + "\t" + refContext.getBaseAsChar() + "\t" + refContext.getLocus().getStart() + "\t" + GATKrecord.getAlignmentStart()  + "\t" + GATKrecord.getReadString() + "\t" + GATKrecord.getReadNegativeStrandFlag());
			
		}

	}

	private void setGoodBases(int offset, ReferenceContext refContext, boolean filterByOriginalQual) {
		byte[] quals = GATKrecord.getOriginalBaseQualities();

		if (GATKrecord.getMappingQuality() >= BAC.MIN_MAPPING_QUALTY_SCORE && quals[offset] >= BAC.MIN_BASE_QUALTY_SCORE
				&& (BAC.USE_BADLY_MATED_READS || (!BadMateFilter.hasBadMate(GATKrecord)) && !GATKrecord.getNotPrimaryAlignmentFlag())) {
			if(BAC.trim2nd){
				if(GATKrecord.getReadPairedFlag() && GATKrecord.getSecondOfPairFlag()){
					if(GATKrecord.getReadNegativeStrandFlag()){
						if(offset >= BAC.trim3 && offset <= GATKrecord.getReadLength() - BAC.trim5){
							goodBase = true;
						}
					}
					else{
						if(offset >= BAC.trim5 && offset <= GATKrecord.getReadLength() - BAC.trim3){
							goodBase = true;
						}
					}
				}
				else{
					goodBase = true;
				}
			}
			else{
				if(GATKrecord.getReadNegativeStrandFlag()){
					if(offset >= BAC.trim3 && offset <= GATKrecord.getReadLength() - BAC.trim5){
						goodBase = true;
					}
				}
				else{
					if(offset >= BAC.trim5 && offset <= GATKrecord.getReadLength() - BAC.trim3){
						goodBase = true;
					}
				}
			}
		}
	

	}
	*/
}
