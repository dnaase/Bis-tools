/**
 * 
 */
package main.java.edu.usc.epigenome.uecgatk.bissnp;


import java.util.HashMap;

import java.util.List;
import java.util.Map.Entry;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 19, 2012 9:12:14 PM
 * 
 */
public class CytosinePatternsUserDefined {

	private HashMap<String, CytosineParameters> contexts;
	/**
	 * 
	 */
	private List<String> cytosineContextsAcquired;
	private int maxCytosinePatternLen = 0;
	private BisulfiteEnums.MethylSNPModel sequencingMode;

	public CytosinePatternsUserDefined(List<String> cytosineContextsAcquired, BisulfiteEnums.MethylSNPModel sequencingMode) {
		this.cytosineContextsAcquired = cytosineContextsAcquired;
		this.sequencingMode = sequencingMode;
		contexts = new HashMap<String, CytosineParameters>();
		setContextsDefinedByUsers();
		setMaxCytosinePatternLen();

	}

	@Override
	public CytosinePatternsUserDefined clone() {
		CytosinePatternsUserDefined returnType = new CytosinePatternsUserDefined(cytosineContextsAcquired, sequencingMode);
		returnType.cytosineContextsAcquired = cytosineContextsAcquired;
		returnType.contexts = contexts;
		returnType.maxCytosinePatternLen = maxCytosinePatternLen;
		returnType.sequencingMode = sequencingMode;
		return returnType;
	}

	/**
	 * @return
	 */
	public HashMap<String, CytosineParameters> getContextDefined() {
		return contexts;
	}

	public int getMaxCytosinePatternLen() {
		return maxCytosinePatternLen;
	}

	private void setContextsDefinedByDefault() {
		if (sequencingMode == BisulfiteEnums.MethylSNPModel.BM) {
			CytosineParameters cp1 = new CytosineParameters();
			cp1.cytosinePosition = 1;
			cp1.cytosineMethylation = 0.5;
			cp1.methylationAutoestimated = true;
			contexts.put("CG", cp1);
			CytosineParameters cp2 = new CytosineParameters();
			cp2.cytosinePosition = 1;
			cp2.cytosineMethylation = 0.5;
			cp2.methylationAutoestimated = true;
			contexts.put("CH", cp2);
			CytosineParameters cp3 = new CytosineParameters();
			cp3.cytosinePosition = 1;
			cp3.cytosineMethylation = 0.5;
			cp3.methylationAutoestimated = true;
			contexts.put("C", cp3);
			/*
			 * CytosineParameters cp2 = new CytosineParameters(); cp2.cytosinePosition = 1;
			 * cp2.cytosineMethylation = 0.5; cp2.methylationAutoestimated = true;
			 * contexts.put("CHG", cp2); CytosineParameters cp3 = new CytosineParameters();
			 * cp3.cytosinePosition = 1; cp3.cytosineMethylation = 0.5; cp3.methylationAutoestimated
			 * = true; contexts.put("CHH", cp3); CytosineParameters cp4 = new CytosineParameters();
			 * cp4.cytosinePosition = 1; cp4.cytosineMethylation = 0.5; cp4.methylationAutoestimated
			 * = true; contexts.put("C", cp4);
			 */
		} else if (sequencingMode == BisulfiteEnums.MethylSNPModel.GM) {
			CytosineParameters cp1 = new CytosineParameters();
			cp1.cytosinePosition = 2;
			cp1.cytosineMethylation = 0.5;
			cp1.methylationAutoestimated = true;
			contexts.put("HCG", cp1);
			// CytosineParameters cp2 = new CytosineParameters();
			// cp2.cytosinePosition = 2;
			// cp2.cytosineMethylation = 0.5;
			// cp2.methylationAutoestimated = true;
			// contexts.put("WCG", cp2);
			CytosineParameters cp3 = new CytosineParameters();
			cp3.cytosinePosition = 2;
			cp3.cytosineMethylation = 0.5;
			cp3.methylationAutoestimated = true;
			contexts.put("GCG", cp3);
			CytosineParameters cp4 = new CytosineParameters();
			cp4.cytosinePosition = 2;
			cp4.cytosineMethylation = 0.5;
			cp4.methylationAutoestimated = true;
			contexts.put("GCH", cp4);
			CytosineParameters cp5 = new CytosineParameters();
			cp5.cytosinePosition = 2;
			cp5.cytosineMethylation = 0.5;
			cp5.methylationAutoestimated = true;
			contexts.put("HCH", cp5);
			// CytosineParameters cp6 = new CytosineParameters();
			// cp6.cytosinePosition = 2;
			// cp6.cytosineMethylation = 0.5;
			// cp6.methylationAutoestimated = true;
			// contexts.put("WCH", cp6);
			// CytosineParameters cp7 = new CytosineParameters();
			// cp7.cytosinePosition = 1;
			// cp7.cytosineMethylation = 0.5;
			// cp7.methylationAutoestimated = true;
			// contexts.put("C", cp7);
		}

	}

	private void setContextsDefinedByUsers() {

		if (cytosineContextsAcquired.isEmpty()) {
			setContextsDefinedByDefault();
		}

		for (String content : cytosineContextsAcquired) {

			String[] contentsArray = content.split(","); // -c CH,1,0.01
			CytosineParameters cp = new CytosineParameters();
			cp.cytosinePosition = Integer.parseInt(contentsArray[1]);
			if (contentsArray.length > 2) {
				cp.cytosineMethylation = Double.parseDouble(contentsArray[2]);
				cp.methylationAutoestimated = false;
			} else {
				cp.cytosineMethylation = 0.5;
				cp.methylationAutoestimated = true;
			}
			String cytosineContext = contentsArray[0].toUpperCase();
			contexts.put(cytosineContext, cp);

		}

	}

	private void setMaxCytosinePatternLen() {

		//for (String cytosinePattern : contexts.keySet()) {
		for ( Entry<String, CytosineParameters> contextsEntrySet: contexts.entrySet()) {
			String cytosinePattern = contextsEntrySet.getKey();
			int tmpLength = cytosinePattern.length();
			int cytosinePos = contexts.get(cytosinePattern).cytosinePosition;
			if (Math.max(tmpLength - cytosinePos, cytosinePos - 1) > maxCytosinePatternLen) {
				maxCytosinePatternLen = Math.max(tmpLength - cytosinePos, cytosinePos - 1);
			}
		}

	}

}
