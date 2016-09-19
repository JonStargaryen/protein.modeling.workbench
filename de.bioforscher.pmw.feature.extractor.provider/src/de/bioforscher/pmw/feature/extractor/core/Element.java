package de.bioforscher.pmw.feature.extractor.core;


import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * convenience copy of BioJava's Element enum to look up several features
 * @author S
 *
 */
public enum Element implements Serializable {

    // most frequently used elements first
    H(1, 1, 1.10f, 0.32f, 1, 1, 1, 1, 1, 1.008f, 0, 1, new int[] { 1 }, 2.20f),
    C(6, 2, 1.55f, 0.77f, 4, 4, 4, 4, 4, 12.011f, 2, -4, new int[] { -4, -3, -2, 0, -1, 1, 2, 3, 4 }, 2.55f),
    N(7, 2, 1.40f, 0.75f, 5, 2, 5, 3, 4, 14.007f, 2, -3, new int[] { -3, -2, -1, 0, 1, 2, 3, 4, 5 }, 3.04f),
    O(8, 2, 1.35f, 0.73f, 6, 1, 2, 2, 2, 16.000f, 2, -2, new int[] { -2, -1, 0, 1, 2 }, 3.44f),
    S(16, 3, 1.81f, 1.02f, 6, 2, 6, 2, 6, 32.060f, 10, -2, new int[] { -2, -1, 0, 1, 2, 3, 4, 5, 6 }, 2.58f),
    Se(34, 4, 0.90f, 1.16f, 6, 0, 12, 2, 6, 78.960f, 28, 4, new int[] { -2, 0, 1, 2, 4, 6 }, 2.55f);

    // should these be declared final?
    private int atomicNumber;
    private int period;
    private float VDWRadius; // in Angstroms
    private float covalentRadius; // in Angstroms
    private int valenceElectronCount;
    private int minimumValence;
    private int maximumValence;
    private int commonValence;
    private int maximumCovalentValence;
    private float atomicMass;
    private int coreElectronCount;
    private int oxidationState;
    private int[] allOxidationStates;
    private float paulingElectronegativity;

    private static final Map<String, Element> allElements;

    static {
        allElements = new HashMap<String, Element>();
        for (Element e : Element.values()) {
            allElements.put(e.toString().toLowerCase(), e);
        }
    }

    private Element(int atomicNumber, int period, float VDWRadius, float covalentRadius, int valenceElectronCount,
            int minimumValence, int maximumValence, int commonValence, int maximumCovalentValence, float atomicMass,
            int coreElectronCount, int oxidationState, int[] allOxidationStates, float paulingElectronegativity) {

        this.atomicNumber = atomicNumber;
        this.period = period;
        this.VDWRadius = VDWRadius;
        this.covalentRadius = covalentRadius;
        this.valenceElectronCount = valenceElectronCount;
        this.minimumValence = minimumValence;
        this.maximumValence = maximumValence;
        this.commonValence = commonValence;
        this.maximumCovalentValence = maximumCovalentValence;
        this.atomicMass = atomicMass;
        this.coreElectronCount = coreElectronCount;
        this.oxidationState = oxidationState;
        this.allOxidationStates = allOxidationStates;
        this.paulingElectronegativity = paulingElectronegativity;
    }

    /**
     * Returns a list of all oxidation states the element is found in. The set
     * is by Greenwood and Norman in "Chemistry of the Elements
     * (ISBN:0080379419).
     *
     * @return An array of oxidation states sorted from most negative to most
     *         positive.
     */
    public int[] getAllOxidationStates() {
        return this.allOxidationStates;
    }

    /**
     * Returns the atomic number of this Element.
     *
     * @return the atomic number of this Element.
     */
    public int getAtomicNumber() {
        return this.atomicNumber;
    }

    /**
     * Returns the period in the periodic table of this Element.
     *
     * @return the period in the periodic table of this Element.
     */
    public int getPeriod() {
        return this.period;
    }

    /**
     * Returns the van der Waals radius of this Element.
     *
     * @return the van der Waals radius of this Element, measured in Angstroms.
     */
    public float getVDWRadius() {
        return this.VDWRadius;
    }

    /**
     * Returns the covalent radius of this Element.
     *
     * @return covalent radius, measured in Angstroms.
     */
    public float getCovalentRadius() {
        return this.covalentRadius;
    }

    /**
     * Returns the number of valence electrons for this Element.
     *
     * @return the number of valence electrons for this Element.
     */
    public int getValenceElectronCount() {
        return this.valenceElectronCount;
    }

    /**
     * Returns the minimum valence for this Element.
     *
     * @return the minimum valence of this atom.
     */
    public int getMinimumValence() {
        return this.minimumValence;
    }

    /**
     * Returns the maximum valence for this Element.
     *
     * @return the maximum valence for this Element.
     */
    public int getMaximumValence() {
        return this.maximumValence;
    }

    /**
     * Returns the common valence for this Element.
     *
     * @return the common valence for this Element.
     */
    public int getCommonValence() {
        return this.commonValence;
    }

    /**
     * Returns the maximum valence for this Element.
     *
     * @return the maximum valence of this element.
     */
    public int getMaximumCovalentValence() {
        return this.maximumCovalentValence;
    }

    /**
     * Returns the atomic mass for this Element.
     *
     * @return the atomic mass for this Element, measured in g/mol.
     */
    public float getAtomicMass() {
        return this.atomicMass;
    }

    /**
     * Returns the number of core electrons for this Element.
     *
     * @return number of core electrons for this Element.
     */
    public int getCoreElectronCount() {
        return this.coreElectronCount;
    }

    /**
     * Returns a typical oxidation state for this Element. This information is
     * mostly useful for metals.
     *
     * @return a typical oxidation state for this Element.
     */
    public int getOxidationState() {
        return this.oxidationState;
    }

    /**
     * Returns the Pauling electronegativity for this Element.
     *
     * @return the Pauling electronegativity for this Element.
     */
    public float getPaulingElectronegativity() {
        return this.paulingElectronegativity;
    }

    /**
     * Returns the Element that corresponds to the specified element symbol. The
     * case of the element symbol is ignored. Example: FE, fe, Fe represent
     * iron.
     *
     * @param elementSymbol
     *            element symbol to specify Element.
     * @return the Element specified by the element symbol.
     */
    public static Element valueOfIgnoreCase(String elementSymbol) throws IllegalArgumentException {
        Element e = allElements.get(elementSymbol.toLowerCase());
        if (e != null) {
            return e;
        }
        throw new IllegalArgumentException("Invalid element symbol: " + elementSymbol);
    }

    /**
     * Returns <code>true</code> if this Element is Hydrogen.
     *
     * @return <CODE>true</CODE> if the Element is Hydrogen.
     */
    public boolean isHydrogen() {
        return this == H;
    }

    /**
     * Returns <CODE>true</CODE> is the Element is an not Hydrogen (or an
     * isotope of Hydrogen).
     * <p>
     * This method is the exact opposite of {@link #isHydrogen()}.
     * </p>
     *
     * @return <CODE>true</CODE> is Element is not Hydrogen.
     */
    public boolean isHeavyAtom() {
        return !isHydrogen();
    }
}

