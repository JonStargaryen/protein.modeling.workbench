package de.bioforscher.pmw.model;

/**
 * Used to state the current reconstruction level of {@link Protein} objects.<br />
 * Furthermore, this {@code enum} is used to request reconstruction operations for {@link Protein} objects by specifying the target reconstruction level.
 * 
 * @author S
 *
 */
public enum ReconstructionLevel {
	/**
	 * no reconstruction efforts performed - raw sequence and maybe some predicted features
	 */
	NONE,
	/**
	 * main step in the reconstruction process - translating the composed contact/distance map to spatial coordinates
	 */
	CALPHA,
	/**
	 * backbone atoms are placed
	 */
	BACKBONE,
	/**
	 * side chain atoms are placed
	 */
	SIDECHAIN,
	/**
	 * the coarsely reconstructed protein was furthermore refined/optimized/minimized
	 */
	REFINED,
	/**
	 * the structure is a PDB structure or was successfully validated by the PMW pipeline
	 */
	VALIDATED;
}
