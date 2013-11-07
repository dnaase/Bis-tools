/**
 * 
 */
package edu.usc.epigenome.dmntools.distribution;

import cern.jet.random.Beta;
import cern.jet.random.Binomial;
import cern.jet.random.engine.MersenneTwister64;
import be.ac.ulg.montefiore.run.distributions.RandomDistribution;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 14, 2012 1:10:41 PM
 * 
 */
public class BetaBinomialDistribution implements RandomDistribution {

	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
     * @uml.property name="trials"
     */
    private int trials;

    /**
     * @uml.property name="alpha"
     */
    private double alpha;

    /**
     * @uml.property name="beta"
     */
    private double beta;
    
    public BetaBinomialDistribution() {
    	this(10, 0.2, 0.8);
	}
	
	public BetaBinomialDistribution(int n, double a, double b) {
		 if (n < 1) n = 1;
	        if (a <= 0) a = 1;
	        if (b <= 0) b = 1;
	        // if (a>=b) b=a+1;
	        trials = n;
	        alpha = a;
	        beta = b;
	}

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.distributions.RandomDistribution#generate()
	 * see Random number generation and Monte Carlo methods / James E. Gentle -2nd ed. 2003
	 * page 187, chapter 5.2.7 
	 * from cern.jet.random.Binomial class
	 * validate the random betaBinomial generater by rbetabin.ab() function's source code in R package: "VGAM"
	 */
	@Override
	public double generate() {
		Beta randomBeta =  new Beta(alpha, beta, new MersenneTwister64());
		Binomial randomBinomial = new Binomial(trials, randomBeta.nextDouble(), new MersenneTwister64());
		return randomBinomial.nextInt();
	//	return randomBinomial.nextDouble();
	}

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.distributions.RandomDistribution#probability(double)
	 */
	@Override
	public double probability(double x) {
		int k = (int) Math.rint(x);
		
        if (k < 0 | k > trials) return 0;
        return (gamma(k+alpha)*gamma(trials-k+beta)*gamma(alpha+beta)*gamma(trials+2))/
        	((trials+1)*gamma(alpha+beta+trials)*gamma(alpha)*gamma(beta)*gamma(k+1)*gamma(trials-k+1));

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
	
	 public double mean() {
	        return trials*alpha/(alpha+beta);
	    }

	    /** Specify the variance in close form */
	    public double variance() {
	        return trials*alpha*(trials*(1+alpha)+beta)/((alpha+beta)*(1+alpha+beta)) - mean()*mean();
	    }

	    /** Specify the SD in close form */
	    public double getSD() {
	        return Math.sqrt(variance());
	    }
	    
	    

	    /**
	     * Get the number of trials
	     * @uml.property name="trials"
	     */
	    public int getTrials() {
	        return trials;
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

}
