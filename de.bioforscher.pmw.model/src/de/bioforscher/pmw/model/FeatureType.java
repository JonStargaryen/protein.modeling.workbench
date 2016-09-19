package de.bioforscher.pmw.model;

/**
 * States which features ({@code double} values describing each {@link Residue}) were previously computed for a {@link Protein} structure and can now be visualized, analyzed or be the basis for further computations.<br />
 * It can also be used to request the computation of certain features.
 *
 * @author S
 *
 */
public enum FeatureType {
	MOTIF_ANNOTATION,
	SECONDARY_STRUCTURE,
	ACCESSIBLE_SURFACE_AREA,
	MEMBRANE_TOPOLOGY,
	INTERACTIONS
}
