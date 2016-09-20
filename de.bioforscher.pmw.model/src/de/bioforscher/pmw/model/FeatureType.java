package de.bioforscher.pmw.model;

/**
 * States which features ({@code double} values describing each {@link Residue}) were previously computed for a {@link Protein} structure and can now be visualized, analyzed or be the basis for further computations.<br />
 * It can also be used to request the computation of certain features.
 *
 * @author S
 *
 */
public enum FeatureType {
	MOTIF_ANNOTATION(MotifType.values().length),
	SECONDARY_STRUCTURE(SecondaryStructure.values().length),
	ACCESSIBLE_SURFACE_AREA(),
	MEMBRANE_TOPOLOGY(Topology.values().length),
	INTERACTIONS(InteractionType.values().length);
	
	/**
	 * the value representing continous values
	 */
	public static final int CONTINUOUS_VALUE = 0;
	
	private transient int numberOfDiscreteValues;
	private transient boolean discrete;
	
	private FeatureType() {
		this(CONTINUOUS_VALUE);
	}
	
	private FeatureType(int numberOfDiscreteValues) {
		this.numberOfDiscreteValues = numberOfDiscreteValues;
		this.discrete = (numberOfDiscreteValues != CONTINUOUS_VALUE);
	}
	
	/**
	 * Gives access to the number of discrete values if this {@link FeatureType} represents an <code>enum</code>, otherwise <code>0</code> is returned.
	 * @return 0 if this feature will produce continuous values, otherwise the number of discrete values
	 */
	public int getNumberOfDiscreteValues() {
		return this.numberOfDiscreteValues;
	}
	
	/**
	 * States whether this feature's values are discrete (i.e. represented by an enum such as {@link InteractionType}) or continuous (e.g. in the case of the accessible surface area).
	 * @return
	 */
	public boolean isDiscrete() {
		return this.discrete;
	}
}
