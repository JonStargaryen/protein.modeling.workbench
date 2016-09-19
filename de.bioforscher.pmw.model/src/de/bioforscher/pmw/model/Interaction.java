package de.bioforscher.pmw.model;

import org.osgi.dto.DTO;

/**
 * Describes interactions between residues, helices or sequence motifs.
 * @author S
 *
 */
public class Interaction extends DTO {
	/*
	 * TODO maybe wrap this in some pair object
	 */
	/**
	 * generic reference to both interaction partners - it is up to any class working with this information to retrieve the actually interacting objects by exploiting the given information on this interaction's type
	 */
	public int[] partners;
	public InteractionType type;
	/*
	 * TODO maybe we need some value(s) to describe the nature of this interaction in more detail
	 */
}
