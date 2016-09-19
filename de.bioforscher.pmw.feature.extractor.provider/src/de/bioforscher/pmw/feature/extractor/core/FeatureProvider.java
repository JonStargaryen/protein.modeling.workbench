package de.bioforscher.pmw.feature.extractor.core;

import de.bioforscher.pmw.model.FeatureType;
import de.bioforscher.pmw.model.Protein;

public interface FeatureProvider {
	String ALPHA_CARBON_NAME = "CA";
    String BETA_CARBON_NAME = "CB";
    /**
     * the distance threshold which helix-helix interactions should not exceed
     */
    double HELIX_HELIX_INTERACTION_CB_CUTOFF = 6.0;
    /**
     * the 'grace' distance additional to van-der-Waals-radii which must not be
     * exceeded in order for helix-helix interactions to occur in agreement to
     * the 2nd criterium
     */
    double HELIX_HELIX_INTERACTION_VDW_CUTOFF = 0.6;
    /**
     * the minimal length a helix must feature in order to span the cell
     * membrane - see {@link http://www.ncbi.nlm.nih.gov/books/NBK21570/}
     */
    int MINIMAL_LENGTH_OF_A_MEMBRANE_SPANNING_HELIX = 22;
	
	/**
	 * reports all missing {@link FeatureType} which are needed for the computation of the actual feature
	 * @param protein the protein for which the feature shall be computed and whose prerequirements will be checked
	 * @return all requirements which need to be computed beforehand
	 */
	FeatureType[] checkForMissingRequirements(Protein protein);
	
	/**
	 * @return all FeatureTypes which need to be provided by any {@link Protein} to be processed, the returned array can be empty, when this computation has no prerequirements
	 */
	FeatureType[] getRequiredFeatures();
	
	/**
	 * @return an array of at least 1 entry, containing the Features provided by this class
	 */
	FeatureType[] getProvidedFeatures();
	
	/**
	 * computes the feature provided by this FeatureProvider (yes, truly!)
	 * @param protein the protein to be processed
	 */
	void computeFeature(Protein protein);
}
