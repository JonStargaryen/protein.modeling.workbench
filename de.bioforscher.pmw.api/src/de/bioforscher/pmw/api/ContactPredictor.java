package de.bioforscher.pmw.api;

import de.bioforscher.pmw.model.Protein;

/**
 * Provides capabilities to predict contacts within protein structures. These
 * contact maps can subsequently be transformed to distance maps. Based on a
 * filled distance map, the {@link ReconstructionService} can finally create a
 * {@link Protein} structure which fulfills these predictions as well as
 * possible.
 * 
 * @author S
 *
 */
public interface ContactPredictor {
	/*
	 * values present in the contact matrices
	 */
	/**
	 * describes a unsafe contact (e.g. contradicting the predicted membrane
	 * topology)
	 */
	int UNSAFE_RESIDUE_CONTACT_VALUE = -1;
	/**
	 * describes no predicted contact for this position
	 */
	int NO_RESIDUE_CONTACT_VALUE = 0;
	/**
	 * a predicted residue-residue contact in agreement with other predictions
	 */
	int SAFE_RESIDUE_CONTACT_VALUE = 1;

	/*
	 * values present in the distance matrices
	 */
	/**
	 * the distance between the same atom
	 */
	double ZERO_DISTANCE = 0.0;
	/**
	 * the ca-ca-distance of 2 consecutive amino acids
	 */
	double CA_CA_DISTANCE = 3.8;
	/**
	 * the distance is not known (or determined yet)
	 */
	double UNKNOWN_DISTANCE = Double.MAX_VALUE;
}
