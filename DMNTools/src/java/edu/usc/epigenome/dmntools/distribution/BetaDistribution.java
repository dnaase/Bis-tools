/**
 * 
 */
package edu.usc.epigenome.dmntools.distribution;

import cern.jet.random.Beta;
import cern.jet.random.Gamma;
import cern.jet.random.engine.MersenneTwister64;
import cern.jet.random.engine.RandomEngine;
import be.ac.ulg.montefiore.run.distributions.RandomDistribution;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 15, 2012 5:18:59 PM
 * 
 */
public class BetaDistribution implements RandomDistribution {

	 /**
     * @uml.property name="alpha"
     */
    private double alpha;

    /**
     * @uml.property name="beta"
     */
    private double beta;
    
    private double scale;
    
    
	public BetaDistribution() {
		this(0.5, 0.5);
	}
	
	public BetaDistribution(double a, double b) {

	        if (a <= 0) a = 1;
	        if (b <= 0) b = 1;
	        // if (a>=b) b=a+1;

	        alpha = a;
	        beta=b;
	        scale = logGamma(alpha + beta) - logGamma(alpha) - logGamma(beta);
	}

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.distributions.RandomDistribution#generate()
	 */
	@Override
	public double generate() {
		Beta randomBeta = new Beta(alpha, beta, new MersenneTwister64());
		return randomBeta.nextDouble();
	}

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.distributions.RandomDistribution#probability(double)
	 */
	@Override
	public double probability(double x) {
		if ((x < 0) | (x > 1)) return 0;
        else if ((x == 0) & (alpha == 1)) return beta;
        else if ((x == 0) & (alpha < 1)) return Double.POSITIVE_INFINITY;
        else if ((x == 0) & (alpha > 1)) return 0;
        else if ((x == 1) & (beta == 1)) return alpha;
        else if ((x == 1) & (beta < 1)) return Double.POSITIVE_INFINITY;
        else if ((x == 1) & (beta > 1)) return 0;
        else return Math.exp(scale + (alpha - 1) * Math.log(x) + (beta - 1)
                * Math.log(1 - x));

	}
	
	public double mean() {
        return alpha/(alpha+beta);
    }

    /** Specify the variance in close form */
    public double variance() {
        //return alpha*((1+alpha)+beta)/((alpha+beta)*(1+alpha+beta)) - mean()*mean();
        return alpha * beta
                / ((alpha + beta) * (alpha + beta) * (alpha + beta + 1));
    }

    /** Specify the SD in close form */
    public double getSD() {
        return Math.sqrt(variance());
    }
    
    

    /**
     * Get alpha parameter
     * @uml.property name="alpha"
     */
    public double getAlpha() {
        return alpha;
    }


    /**
     * Get beta parameter
     * @uml.property name="beta"
     */
    public double getBeta() {
        return beta;
    }
	
	public double logGamma(double x) {
		double coef[] = { 76.18009173, -86.50532033, 24.01409822, -1.231739516,
				0.00120858003, -0.00000536382 };
		double step = 2.50662827465, fpf = 5.5, t, tmp, ser;
		t = x - 1;
		tmp = t + fpf;
		tmp = (t + 0.5) * Math.log(tmp) - tmp;
		ser = 1;
		for (int i = 1; i <= 6; i++) {
			t = t + 1;
			ser = ser + coef[i - 1] / t;
		}
		return tmp + Math.log(step * ser);
	}

	/** This method computes the gamma function. */
	public double gamma(double x) {
		return Math.exp(logGamma(x));
	}
	
	 /** overwrites the method in distribution for estimating parameters */
    public void paramEstimate(double[] distData) {

        double sMean = sampleMean(distData);
        double sVar = sampleVar(distData, sMean);
        // estimate alpha and beta a,b from sample statistics
        double A = minDouble(distData);
        double B = maxDouble(distData);
        double q = (sMean - A) / (B - A);
        //alpha = ( (B-A)*(B-A)*(1-q)*q*q/sVar) - 1;
        //beta = alpha*(1-q)/q;
        beta = (1 - sMean) * (1 - sMean) * sMean / sVar - (1 - sMean);
        alpha = sMean * beta / (1 - sMean);

    }
    
 // estimate sample mean from data
 	public double sampleMean(double[] array) {
 		double mean = 0;
 		double sumX = array[0];
 		for (int i = 1; i < array.length; i++)
 			sumX += array[i];
 		mean = sumX / array.length;
 		return mean;
 	}

 	public double sampleVar(double[] array, double mean) {
 		double variance = 1;
 		double sumX2 = array[0] * array[0];
 		for (int i = 1; i < array.length; i++)
 			sumX2 += array[i] * array[i];
 		//variance=mean*mean-sumX2/(array.length*array.length);
 		variance = -mean * mean + sumX2 / array.length;
 		return variance;
 	}
 	
 	public double maxDouble(double[] X) {
        double iMax = X[0];
        for (int i = 0; i < X.length; i++) {
            if (X[i] > iMax) iMax = X[i];
        }
        return iMax;
    }

    public double minDouble(double[] X) {
        double iMin = X[0];
        for (int i = 0; i < X.length; i++) {
            if (X[i] < iMin) iMin = X[i];
        }
        return iMin;
    }
    
    /** Compute the maximum getDensity */
    public double getMaxDensity() {
        double mode;
        if (alpha < 1) mode = 0.01;
        else if (beta <= 1) mode = 0.99;
        else mode = (alpha - 1) / (alpha + beta - 2);
       // System.err.println("mode: " + mode);
        return probability(mode);
    }

}
