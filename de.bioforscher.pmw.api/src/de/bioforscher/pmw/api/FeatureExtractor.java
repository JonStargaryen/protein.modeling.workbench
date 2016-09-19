package de.bioforscher.pmw.api;

import de.bioforscher.pmw.model.FeatureType;
import de.bioforscher.pmw.model.Protein;

/**
 * A service to compute features (such as accessible surface area, secondary structure, ...) for a given protein structure.<br />
 * They can either be:<br />
 * <ul>
 * <li><b>Annotations</b>: Require 3D coordinates to be present. Quite reliable, but results from different algorithms may still differ.</li>
 * <li><b>Predictions</b>: Mainly sequence-based, thus, with little requirements. Potentially, erroneous values.</li>
 * </ul>
 * Available features are defined by {@link FeatureType} - thus, feature types are only loosely coupled to algorithms. This mapping is defined in the implementation of {@link FeatureExtractor}.<br 7>
 * Results will be written to the feature list of each {@link Residue}. The computed features will be presented via the availableFeatures field of {@link Protein}.
 * @author S
 *
 */
public interface FeatureExtractor {	
	/**
	 * Provides a way to compute several Features for the given protein.
	 * @param protein the protein to be processed - computed features will be added and the corresponding feature flags will be updated
	 * @param featuresToCompute the set of values to be computed
	 */
	void computeFeatures(Protein protein, FeatureType... featuresToCompute);
}
