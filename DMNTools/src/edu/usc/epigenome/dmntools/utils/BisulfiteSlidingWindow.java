/**
 * 
 */
package edu.usc.epigenome.dmntools.utils;

import java.util.LinkedList;


/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Oct 16, 2013 6:48:12 PM
 * 
 */
public class BisulfiteSlidingWindow {

		public LinkedList<BisulfiteGenomicLocHmm> windowList;
		private int numC = 0;
		private int numCT = 0;
		private int minCT;
		private int minGch;
		private long minLen;
		private int segHmmState;
		
		
		public BisulfiteSlidingWindow(int minCT, int minGch, long minLen, int segHmmState){
			windowList = new LinkedList<BisulfiteGenomicLocHmm>();
			this.minCT = minCT;
			this.minGch = minGch;
			this.minLen = minLen;
			this.segHmmState = segHmmState;
			
		}
		
		public void addLast(BisulfiteGenomicLocHmm data){
			windowList.offerLast(data);
			if(data.hmmState != segHmmState){
				numCT += data.numCT;
				numC += data.numC;
			}
			
			
		}
		
		public void addFirst(BisulfiteGenomicLocHmm data){
			windowList.offerFirst(data);
			if(data.hmmState != segHmmState){
				numCT += data.numCT;
				numC += data.numC;
			}
			
			
		}
		
		public LinkedList<BisulfiteGenomicLocHmm> addLast(BisulfiteGenomicLocHmm data, boolean automateRemoveFirstBatch){
			LinkedList<BisulfiteGenomicLocHmm> dataList = new LinkedList<BisulfiteGenomicLocHmm>();
			addLast(data);
			while(getGchNum() > minGch && getLength() > minLen && getCtReadsNum() > minCT){
				dataList.offerLast(removeFirst());
			}
			if(!dataList.isEmpty()){
				addFirst(dataList.pollLast());
			}
			
			return dataList;
		}
		
		public BisulfiteGenomicLocHmm removeLast(){
			BisulfiteGenomicLocHmm data = windowList.pollLast();
			if(data.hmmState != segHmmState){
				numCT -= data.numCT;
				numC -= data.numC;
			}
			
			return data;
		}
		
		public BisulfiteGenomicLocHmm removeFirst(){
			BisulfiteGenomicLocHmm data = windowList.pollFirst();
			if(data.hmmState != segHmmState){
				numCT -= data.numCT;
				numC -= data.numC;
			}
			
			return data;
		}
		
		public double getMean(){
			
			return (double)numC/(double)numCT;
		}
		
		public int getGchNum(){
			return windowList.size();
		}
		
		public int getLength(){
			if(!windowList.isEmpty()){
				BisulfiteGenomicLocHmm last = windowList.peekLast();
				BisulfiteGenomicLocHmm first = windowList.peekFirst();
				if(last.contig.equalsIgnoreCase(first.contig))
					return last.end - first.start;
			}
			return 0;
		}
		
		public int getCtReadsNum(){

			return numCT;
		}
		
		public int getCReadsNum(){
			
			return numC;
		}
		

}
