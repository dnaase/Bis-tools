/**
 * 
 */
package edu.usc.epigenome.dmntools.hmm;

import java.text.NumberFormat;

import be.ac.ulg.montefiore.run.jahmm.Centroid;
import be.ac.ulg.montefiore.run.jahmm.CentroidFactory;
import be.ac.ulg.montefiore.run.jahmm.CentroidObservationReal;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 15, 2012 12:01:39 PM
 * 
 */
public class ObservationMethy extends Observation
implements CentroidFactory<ObservationMethy>{

	public int coverage;
	public double distance;
	final public double value;
	/**
	 * @param value
	 */
	public ObservationMethy(double value) {
		this.value = value;

	}
	
	/**
	 * @param value
	 */
	public void setCoverage(int coverage) {
		this.coverage = coverage;

	}

	/**
	 * @param normalized distance to the next element
	 */
	public void setDistance(double distance) {
		this.distance = distance;

	}
	
	/**
	 * Returns the centroid matching this observation.
	 *
	 * @return The corresponding observation.
	 */
	public Centroid<ObservationMethy> factor()
	{
		return new CentroidObservationMethy(this);
	}

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.jahmm.Observation#toString(java.text.NumberFormat)
	 */
	@Override
	public String toString(NumberFormat numberFormat) {
		// TODO Auto-generated method stub
		return numberFormat.format(value);
	}
}
