package de.bioforscher.pmw.model;

import java.util.Arrays;
import java.util.List;

import org.osgi.dto.DTO;

/**
 * A wrapper used to expose constants (especially those defined by the model's enums) via REST.
 * @author S
 *
 */
public class APIConstants extends DTO {
	/** features (e.g. ASA, SSE, topology) which can be computed */
	public List<FeatureType> features = Arrays.asList(FeatureType.values());
	/** steps for the reconstruction process */
	public List<ReconstructionLevel> reconstructionLevels = Arrays.asList(ReconstructionLevel.values());
	/** known secondary structure elements - defined by BioJava's secondary structure annotator */
	public List<SecondaryStructure> secondaryStructures = Arrays.asList(SecondaryStructure.values());
	/**  */
	public List<Topology> topologies = Arrays.asList(Topology.values());
	/**  */
	public List<MotifType> motifTypes = Arrays.asList(MotifType.values());
	/**  */
	public List<InteractionType> interactionTypes = Arrays.asList(InteractionType.values());
}
