/**
 * 
 */
package edu.usc.epigenome.dmntools.hmm;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.List;

import edu.usc.epigenome.dmntools.distribution.OpdfBetaFactory;

import be.ac.ulg.montefiore.run.jahmm.CentroidFactory;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.KMeansCalculator;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfFactory;
import be.ac.ulg.montefiore.run.jahmm.ViterbiCalculator;
import be.ac.ulg.montefiore.run.jahmm.learn.KMeansLearner;


/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Dec 11, 2012 11:59:36 PM
 * 
 */
public class NhmmKMeansLearner<O extends Observation & CentroidFactory<? super O>> {

	private Clusters<O> clusters;
	private int nbStates;
	private List<? extends List<? extends O>> obsSeqs;
	private OpdfFactory<? extends Opdf<O>> opdfFactory;
	private boolean terminated;

	
	/**
	 * Initializes a K-Means algorithm implementation.  This algorithm
	 * finds a HMM that models a set of observation sequences.
	 *
	 * @param nbStates  The number of states the resulting HMM will be made of.
	 * @param opdfFactory A class that builds the observation probability
	 *                    distributions associated to the states of the HMM.
	 * @param sequences A vector of observation sequences.  Each observation
	 *                sequences is a vector of
	 *                {@link be.ac.ulg.montefiore.run.jahmm.Observation
	 *                observations} compatible with the 
	 *                {@link be.ac.ulg.montefiore.run.jahmm.CentroidFactory
	 *                k-means algorithm}.
	 */
	public NhmmKMeansLearner(int nbStates,
			OpdfFactory<? extends Opdf<O>> opdfFactory,
			List<? extends List<? extends O>> sequences)
	{	
		this.obsSeqs = sequences;
		this.opdfFactory = opdfFactory;
		this.nbStates = nbStates;
		
		List<? extends O> observations = flat(sequences);
		clusters = new Clusters<O>(nbStates, observations);
		terminated = false;
	}
	



	/**
	 * @param sTATES
	 * @param opdfBetaFactory
	 * @param seqs
	 */
	public NhmmKMeansLearner(int nbStates, OpdfBetaFactory opdfBetaFactory, ArrayList<LinkedList<ObservationMethy>> sequences) {
		this.obsSeqs = (List<? extends List<? extends O>>) sequences;
		this.opdfFactory = (OpdfFactory<? extends Opdf<O>>) opdfBetaFactory;
		this.nbStates = nbStates;
		//System.err.println(sequences.get(0).get(0).value);
		List<? extends O> observations = flat((List<? extends List<? extends O>>)sequences);
		//for(O o : observations)
		//	System.err.println(o);
		clusters = new Clusters<O>(nbStates, observations);
		terminated = false;
	}




	/**
	 * Performs one iteration of the K-Means algorithm.
	 * In one iteration, a new HMM is computed using the current clusters, and
	 * the clusters are re-estimated using this HMM.
	 *
	 * @return A new, updated HMM.
	 */
	public Nhmm<O> iterate()
	{	
		Nhmm<O> nhmm = new Nhmm<O>(nbStates, opdfFactory);
	//	double[] pi = new double[] {0.33,0.33, 0.33};
	//	double[][] Aij = new double[][] {{0.5,0., 0.5},{0.5,0.5, 0.}, {0.33,0.33, 0.33}};
	//	ArrayList<Opdf<O>> opdfs = new ArrayList<Opdf<O>>(nbStates);
	//	for(int i = 0; i< nbStates; i++)
	//		opdfs.add(opdfFactory.factor());

	//	Hmm<O> hmm = new Hmm<O>(pi, Aij, opdfs);
	//	hmm.setPi(0, 0.8);
	//	hmm.setPi(1, 0.2);
	//	hmm.setAij(0, 1, 0.1);
	//	hmm.setAij(0, 0, 0.9);
	///	hmm.setAij(1, 0, 0.1);
	//	hmm.setAij(1, 1, 0.9);
		//hmm.setOpdf(0, new OpdfBeta(0.3,1.0));
		//System.out.println("Before K-means-hmm:" + hmm);
		learnPi(nhmm);
		System.out.println("learnPi:");
		learnAij(nhmm);
		System.out.println("learnAij:");
	//	System.out.println("K-means-hmm-before-learnOpdf:" + hmm);
		learnOpdf(nhmm);
		System.out.println("learnOpdf:");
		
		terminated = optimizeCluster(nhmm);
		System.out.println("K-means-nhmm:" + nhmm);
		return nhmm;
	}
	
	
	/**
	 * Returns <code>true</code> if the algorithm has reached a fix point,
	 * else returns <code>false</code>.
	 */
	public boolean isTerminated()
	{
		return terminated;
	}
	
	
	/**
	 * Does iterations of the K-Means algorithm until a fix point is reached.
	 * 
	 * @return The HMM that best matches the set of observation sequences given
	 *         (according to the K-Means algorithm).
	 */
	public Nhmm<O> learn()
	{	
		Nhmm<O> nhmm;
		int i = 0;
		do {
			System.out.println("iteration:" + i);
			nhmm = iterate();
		
		i++;}
		while(!isTerminated());
		
		return nhmm;
	}
	
	
	private void learnPi(Nhmm<?> hmm)
	{	
		double[] pi = new double[nbStates];
		
		for (int i = 0; i < nbStates; i++)
			pi[i] = 0.;
		
		for (List<? extends O> sequence : obsSeqs)
			pi[clusters.clusterNb(sequence.get(0))]++;
	//	for (List<? extends O> sequence : obsSeqs)
	//		for (int i = 1; i < sequence.size(); i++) 
	//			pi[clusters.clusterNb(sequence.get(i))]++;
	//	double sum = 0;
	//	for (int i = 0; i < hmm.nbStates(); i++)
	//			sum += pi[i];
	//	for (int i = 0; i < hmm.nbStates(); i++)
	//		hmm.setPi(i, pi[i] / sum);
		for (int i = 0; i < nbStates; i++)
			hmm.setPi(i, pi[i] / obsSeqs.size());
	}
	
	
	private void learnAij(Nhmm<O> hmm)
	{	
		for (int i = 0; i < hmm.nbStates(); i++)
			for (int j = 0; j < hmm.nbStates(); j++)
				hmm.setAij(i, j, 0.);
	//	int count = 0;
		for (List<? extends O> obsSeq : obsSeqs) {
			if (obsSeq.size() < 2)
				continue;
			
			int first_state;
			int second_state = clusters.clusterNb(obsSeq.get(0));
			for (int i = 1; i < obsSeq.size(); i++) {
				first_state = second_state;
				second_state =
					clusters.clusterNb(obsSeq.get(i));
			//	if (first_state == 0 || second_state == 0){
			//		count++;
			//	}
					//System.err.println( " i:" + i + "\tfirst_state:" + first_state + "\tsecond_state" + second_state + "\tseq" + obsSeq.get(i) + "\t" + (hmm.getAij(first_state, second_state)+1.));
				hmm.setAij(first_state, second_state,
						hmm.getAij(first_state, second_state)+1.);
			}
		}
		//System.err.println(hmm.getAij(0, 1) + "\t" + hmm.getAij(0, 2) + "\t" + hmm.getAij(0, 0) + "\tcount:" + count);
		/* Normalize Aij array */
		for (int i = 0; i < hmm.nbStates(); i++) {
			double sum = 0;
			
			for (int j = 0; j < hmm.nbStates(); j++){
			//	System.err.println(i + "," + j + ":" + hmm.getAij(i, j));
				sum += hmm.getAij(i, j);
			}
			//System.err.println( "i:" + i + "sum:" + sum);	
				
			if (sum == 0.)
				for (int j = 0; j < hmm.nbStates(); j++) 
					hmm.setAij(i, j, 1. / hmm.nbStates());     // Arbitrarily
			else
				for (int j = 0; j < hmm.nbStates(); j++)
					hmm.setAij(i, j, hmm.getAij(i, j) / sum);
		}
	}
	
	
	private void learnOpdf(Nhmm<O> hmm)
	{
		for (int i = 0; i < hmm.nbStates(); i++) {
			Collection<O> clusterObservations = clusters.cluster(i);
			//System.err.println("cluster " + i + "\t size:" + clusterObservations.size());
			if (clusterObservations.isEmpty())
				hmm.setOpdf(i, opdfFactory.factor());
			else
				hmm.getOpdf(i).fit(clusterObservations);
		}
	}
	
	
	/* Return true if no modification */
	private boolean optimizeCluster(Nhmm<O> hmm)
	{	
		boolean modif = false;
		
		for (List<? extends O> obsSeq : obsSeqs) {
			NhmmViterbiCalculator vc = new NhmmViterbiCalculator(obsSeq, hmm);
			//ViterbiCalculator vc = new ViterbiCalculator(obsSeq, hmm);
			int states[] = vc.stateSequence();
			
			for (int i = 0; i < states.length; i++) {
				O o = obsSeq.get(i);
				//if(states[i] == 2)
				//System.err.print(states[i]);
				if (clusters.clusterNb(o) != states[i]) {
					modif = true;
					System.err.print(clusters.clusterNb(o) + "\t" + states[i]);
					clusters.remove(o, clusters.clusterNb(o));
					clusters.put(o, states[i]);
				}
			}

			//	System.err.println();
		}
		
		return !modif;
	}
	
	
	static <T> List<T> flat(List<? extends List<? extends T>> lists)
	{	
		List<T> v = new ArrayList<T>();
		
		for (List<? extends T> list : lists)
			v.addAll(list);
		
		return v;
	}
}


/*
 * This class holds the matching between observations and clusters.
 */
class Clusters<O extends CentroidFactory<? super O>>
{	
	class Value
	{
		private int clusterNb;
		
		Value(int clusterNb)
		{
			this.clusterNb = clusterNb;
		}
		
		void setClusterNb(int clusterNb)
		{
			this.clusterNb = clusterNb;
		}
		
		int getClusterNb()
		{
			return clusterNb;
		}
	}
	
	
	private Hashtable<O,Value> clustersHash;
	private ArrayList<Collection<O>> clusters;
	
	
	public Clusters(int k, List<? extends O> observations)
	{
		
		clustersHash = new Hashtable<O,Value>();
		clusters = new ArrayList<Collection<O>>();
		
		KMeansCalculator<O> kmc = new KMeansCalculator<O>(k, observations);
		
		for (int i = 0; i < k; i++) {
			Collection<O> cluster = kmc.cluster(i);
			clusters.add(cluster);
			
			for (O element : cluster) 
				clustersHash.put(element, new Value(i));
		}
	}
	
	
	public boolean isInCluster(Observation o, int clusterNb)
	{
		return clusterNb(o) == clusterNb;
	}
	
	
	public int clusterNb(Observation o)
	{
		return clustersHash.get(o).getClusterNb();
	}
	
	
	public Collection<O> cluster(int clusterNb)
	{
		return clusters.get(clusterNb);
	}
	
	
	public void remove(Observation o, int clusterNb)
	{
		clustersHash.get(o).setClusterNb(-1);
		clusters.get(clusterNb).remove(o);
	}
	
	
	public void put(O o, int clusterNb)
	{
		clustersHash.get(o).setClusterNb(clusterNb);
		clusters.get(clusterNb).add(o);
	}
}
