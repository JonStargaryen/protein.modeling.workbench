package de.bioforscher.pmw.model;

/**
 * Represents secondary structure elements as annotated by DSSP. When compared, the types are sorted in the declaration order of the enum,
 * which is the DSSP preference of type assignment.
 * @author S
 *
 */
public enum SecondaryStructure {
	COIL,
	BEND,
	TURN,
	PIHELIX,
	THREE10HELIX,
	BRIDGE,
	EXTENDED,
	ALPHA_HELIX;
	
	public boolean isHelixType() {
		return this.name().contains("HELIX");
	}
}
