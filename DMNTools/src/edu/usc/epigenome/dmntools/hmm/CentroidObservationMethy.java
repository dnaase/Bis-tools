/**
 * 
 */
package edu.usc.epigenome.dmntools.hmm;

import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.Centroid;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Dec 12, 2012 3:25:52 PM
 * 
 */
public class CentroidObservationMethy implements Centroid<ObservationMethy> {

	
	private double value;
	
	
	public CentroidObservationMethy(ObservationMethy o)
	{
		this.value = o.value;
	}
	
	
	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.jahmm.Centroid#reevaluateAdd(java.lang.Object, java.util.List)
	 */
	@Override
	public void reevaluateAdd(ObservationMethy e, List<? extends ObservationMethy> v) {
		value =	(value * (double) v.size() + e.value) / (v.size()+1.);

	}

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.jahmm.Centroid#reevaluateRemove(java.lang.Object, java.util.List)
	 */
	@Override
	public void reevaluateRemove(ObservationMethy e, List<? extends ObservationMethy> v) {
		value = ((value * (double) v.size()) - e.value) / (v.size()-1.);

	}

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.jahmm.Centroid#distance(java.lang.Object)
	 */
	@Override
	public double distance(ObservationMethy e) {
		// TODO Auto-generated method stub
		return Math.abs(e.value - value);
	}

}
