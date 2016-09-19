package de.bioforscher.pmw.model;

import org.osgi.dto.DTO;

/**
 * Defines a sequence motif (such as GG4), enumerated by {@link DefinedMotif}.
 * @author S
 *
 */
public class Motif extends DTO {
	public String sequence;
	/**
	 * reference to the type this motif
	 */
	public DefinedMotif definition;
	/**
	 * the topology of this motif encoded as <code>double</code> value - this mapping is defined by the {@link Topology} enum
	 */
	public double topology;
	/*
	 * TODO maybe this should be similar to the pair annotation in the interaction object
	 */
	public int startResidueId;
	public int endResidueId;
}
