package de.bioforscher.pmw.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.osgi.dto.DTO;

public class Residue extends DTO {
	public String aminoAcid;
	public List<Atom> atoms;
	/**
	 * Entries consist of a value in index 0 and a normalized value accessible by index 1.
	 * This should be some EnumMap<FeatureType, Feature>, however MongoDB dies on complex keys and complex entries. This is still better than a List.
	 */
	public Map<String, double[]> features;
	/** never actually used, maybe drop it */
	public String insertionCode;
	/** the model/API uses this custom field of residues to identify them in the structure (it is unique across all chains, e.g. "C-123" could be utilized alternatively) */
	public int residueId;
	public int residueNumber;
	
	public Residue() {
		this.atoms = new ArrayList<>();
		// we could use a EnumMap here, however the MongoDB provider is not capable of processing their keys
//		this.features = new EnumMap<>(FeatureType.class);
		this.features = new HashMap<>();
	}
}