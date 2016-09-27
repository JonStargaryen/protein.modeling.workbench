package de.bioforscher.pmw.feature.extractor.sse;

public enum DSSPSecondaryStructureElement {
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

	public boolean isStrandType() {
		return this.equals(BRIDGE) || this.equals(EXTENDED);
	}
}
