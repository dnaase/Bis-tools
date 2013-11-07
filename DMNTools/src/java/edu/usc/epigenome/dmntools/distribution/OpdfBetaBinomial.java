/**
 * 
 */
package edu.usc.epigenome.dmntools.distribution;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Collection;

import org.apache.commons.math3.special.Gamma;

import edu.usc.epigenome.dmntools.hmm.ObservationMethy;

import be.ac.ulg.montefiore.run.distributions.GaussianDistribution;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 14, 2012 7:21:33 PM
 * 
 */
public class OpdfBetaBinomial implements Opdf<ObservationMethy> {

	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private BetaBinomialDistribution distribution;
	
	private double epi = 1e-5;
	/**
	 * 
	 */
	public OpdfBetaBinomial() {
		distribution = new BetaBinomialDistribution();
	}
	
	public OpdfBetaBinomial(int n, double a, double b) {
		distribution = new BetaBinomialDistribution(n,a,b);
	}
	
	public OpdfBetaBinomial(Collection<? extends ObservationMethy> co) {
		double[] weights = new double[co.size()];
		Arrays.fill(weights, 1. / co.size());
		fit(co, weights, true);
		
	}

	/**
	 * Returns this distribution's mean value.
	 *
	 * @return This distribution's mean value.
	 */
	public double mean() 
	{
		return distribution.mean();
	}
	
	
	/**
	 * Returns this distribution's variance.
	 *
	 * @return This distribution's variance.
	 */
	public double variance()
	{
		return distribution.variance();
	}
	
	
	public double probability(ObservationMethy o) 
	{	
		return distribution.probability(o.value);
	}
	
	
	public ObservationMethy generate()
	{
		return new ObservationMethy(distribution.generate());
	}

	

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.jahmm.Opdf#fit(O[])
	 */
	@Override
	public void fit(ObservationMethy... oa) {
		fit(Arrays.asList(oa));
		
	}

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.jahmm.Opdf#fit(java.util.Collection)
	 */
	@Override
	public void fit(Collection<? extends ObservationMethy> co) {
		double[] weights = new double[co.size()];
		Arrays.fill(weights, 1. / co.size());
		
		fit(co, weights);
		
	}

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.jahmm.Opdf#fit(O[], double[])
	 */
	@Override
	public void fit(ObservationMethy[] o, double[] weights) {
		fit(Arrays.asList(o), weights);
		
	}

	/* Molaro 2011 Cell paper!!
	 * @see be.ac.ulg.montefiore.run.jahmm.Opdf#fit(java.util.Collection, double[])
	 */
	@Override
	public void fit(Collection<? extends ObservationMethy> co, double[] weights) {
		if (co.isEmpty() || co.size() != weights.length)
			throw new IllegalArgumentException();
		
		// Compute trials
		double trials = 0.;
		int i = 0;
		for (ObservationMethy o : co){
			
			trials += o.coverage * weights[i++];
		//	System.err.println(weights[i-1] + "\t" + d  + "\t" + o.value + "\t" + mean);
			
		}
		
		// Compute alpha
		double alpha = 0.;
		double logP = 0.;
		i = 0;
		for (ObservationMethy o : co){
			
			//logP += Math.log(o.value);
			logP += Math.log(o.value) * weights[i++];
		//	System.err.println(logP + "\tvalue: " + o.value);
		//	System.err.println(weights[i-1] + "\t" + d  + "\t" + o.value + "\t" + mean);
			
		}
		//logP /= co.size();
		logP += digammaFunction(distribution.getAlpha() + distribution.getBeta());
		alpha = bisectionToFindFunctionRoot(logP);
		System.err.println(logP + "\talpha: " + alpha);
		// Compute beta
		double beta = 0.;
		double log1minusP = 0.;
		i = 0;
		for (ObservationMethy o : co) {
			//log1minusP += Math.log(1-o.value);
			log1minusP += Math.log(1-o.value) * weights[i++];
		//	System.err.println("log1minusP: " + log1minusP + "\tvalue: " + o.value);
			if(Double.isNaN(log1minusP))
				System.exit(1);
		}
		//log1minusP /= co.size();
		log1minusP += digammaFunction(distribution.getAlpha() + distribution.getBeta());
		beta = bisectionToFindFunctionRoot(log1minusP);
		System.err.println(co.size() + "\t" + log1minusP + "\tbeta: " + beta);
		distribution = new BetaBinomialDistribution((int)trials, alpha, beta);
		
	}
		
	/*
	 * initiate beta-binomial distribution
	 */
		public void fit(Collection<? extends ObservationMethy> co, double[] weights, boolean initiate) {
			if (co.isEmpty() || co.size() != weights.length)
				throw new IllegalArgumentException();
			
			// Compute trials
			double trials = 0.;
			int i = 0;
			for (ObservationMethy o : co){
				
				trials += o.coverage * weights[i++];
			//	System.err.println(weights[i-1] + "\t" + d  + "\t" + o.value + "\t" + mean);
				
			}
			
			// Compute alpha
			double alpha = 0.;
			double logP = 0.;
			i = 0;
			for (ObservationMethy o : co){
				
				//logP += Math.log(o.value);
				logP += Math.log(o.value)  * weights[i++];
			//	System.err.println(weights[i-1] + "\t" + d  + "\t" + o.value + "\t" + mean);
				
			}
			//logP /= co.size();
			alpha = bisectionToFindFunctionRoot(logP);
			
			// Compute beta
			double beta = 0.;
			double log1minusP = 0.;
			i = 0;
			for (ObservationMethy o : co) {
				//log1minusP += Math.log(1-o.value);
				log1minusP += Math.log(1-o.value) * weights[i++];
			}
			//log1minusP /= co.size();
			beta = bisectionToFindFunctionRoot(log1minusP);
			
			distribution = new BetaBinomialDistribution((int)trials, alpha, beta);
			
		}
	
		/*
		 * funtion is: exp(x)-theta = 0;
		 * 
		 */
		private double bisectionToFindFunctionRoot(double theta){
			double x = 0, a = 0.01, b = 1000000;
			
		    double dx = b-a;
		    int k = 0;
		    while (Math.abs(dx) > epi ) {
		      x = ((a+b)/2);
		      if (((digammaFunction(a)-theta)*(digammaFunction(x)-theta)) < 0) {
		      //if (((Gamma.digamma(a)-theta)*(Gamma.digamma(x)-theta)) < 0) {
		        b  = x;
		        dx = b-a;
		      }
		      else {
		        a = x;
		        dx = b-a;
		      }
		    //  System.err.println(digammaFunction(a) + "\t" + digammaFunction(b) + "\ttheta:" + theta + "\ta: " + a + "\tb: " + b + "\tx: " + x + "\tk: " + k);
		     k++;
		    }
		  //  System.err.println(digammaFunction(a) + "\t" + digammaFunction(b) + "\ttheta:" + theta + "\ta: " + a + "\tb: " + b + "\tx: " + x + "\tk: " + k);
		    return x;
		}
		
		
		
		/*
		 * digamma function is approximated by Bernardo, JosŽ M. (1976). "Algorithm AS 103 psi(digamma function) computation". Applied Statistics 25: 315Ð317.
		 * http://en.wikipedia.org/wiki/Digamma_function#cite_ref-2 
		 */
		private double digammaFunction(double x){
			return Math.log(x) - 1/(2*x) - 1/(12*Math.pow(x, 2)) + 1/(120*Math.pow(x, 4)) - 1/(252*Math.pow(x, 6));
		}
		
		private double inverseDigammaFunction(double x){
			return Math.exp(x) + 0.2;
		}
		
	public OpdfBetaBinomial clone()
	{
		try {
			return (OpdfBetaBinomial) super.clone();
		} catch(CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
	}

	/* should alpha and beta more useful..
	 * @see be.ac.ulg.montefiore.run.jahmm.Opdf#toString(java.text.NumberFormat)
	 */
	@Override
	public String toString(NumberFormat numberFormat) {
		return "BetaBinomial distribution --- " +
				"Mean: " + distribution.mean() +
						//numberFormat.format(distribution.mean()) +
				" Variance " + distribution.variance() + 
						//numberFormat.format(distribution.variance());
		" Alpha: " + distribution.getAlpha() +

		" Beta " + distribution.getBeta() + 
		" Trial " + distribution.getTrials();
		
	}

}
