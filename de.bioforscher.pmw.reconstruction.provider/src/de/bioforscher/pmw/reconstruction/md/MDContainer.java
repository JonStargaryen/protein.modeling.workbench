package de.bioforscher.pmw.reconstruction.md;

import java.util.HashMap;
import java.util.Map;

import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.model.Atom;
import de.bioforscher.pmw.model.Protein;

/**
 * A wrapper for {@link Protein} objects so MD-related values like forces and
 * velocities.
 * 
 * @author S
 *
 */
public class MDContainer {
	/**
	 * reference to the {@link ModelConverter} instance
	 */
	private ModelConverter modelConverter;

	/**
	 * defines the cooling strategy of the system
	 * 
	 * @author S
	 *
	 */
	public enum CoolingProtocol {
		LINEAR, SIGMOID, EXPONENT;
	}

	public enum Thermostat {
		BUSSI;
	}

	public enum Barostat {
		BERENDSEN;
	}

	/*
	 * constants
	 */
	private static final double DEFAULT_DT = 0.001;
	private static final int DEFAULT_NUMBER_OF_COOLING = 2000;
	private static final CoolingProtocol DEFAULT_COOLING_PROTOCOL = CoolingProtocol.LINEAR;

	/*
	 * variables
	 */
	/**
	 * reference to the {@link Protein}
	 */
	private Protein protein;
	/**
	 * keeps data - it is mapped from {@link Atom} objects to some yet to define
	 * {@link MDData} object (forces, velocities etc)
	 */
	private Map<Atom, MDData> data;
	/**
	 * the length of the dynamics time step in picoseconds (TODO what is the
	 * time scale?) - default 0.001
	 */
	private final double dt;
	/**
	 * start temperature in K - default 1000
	 */
	private final double startTemperatur;
	/**
	 * end temperature in K - default 0
	 */
	private final double endTemperatur;

	public double getDt() {
		return dt;
	}

	public double getStartTemperatur() {
		return startTemperatur;
	}

	public double getEndTemperatur() {
		return endTemperatur;
	}

	public int getNumberOfCoolingSteps() {
		return numberOfCoolingSteps;
	}

	public CoolingProtocol getCoolingProtocol() {
		return coolingProtocol;
	}

	/**
	 * number of dynamics steps for the cooling protocol - default 2000
	 */
	private final int numberOfCoolingSteps;
	/**
	 * decide which annealing cooling protocol to use
	 */
	private final CoolingProtocol coolingProtocol;
	/**
	 * the current temperature of the system
	 */
	private double loose;
	private double tight;
	private double currentTemperature;
	private double tauTemperature;
	@SuppressWarnings("unused")
	private boolean isobaric;

	/**
	 * Construct a new wrapping instance for a given protein. Initializes all
	 * values to default values.
	 * 
	 * @param protein
	 *            The coordinates to convert to a format processable by the MD
	 *            implementation.
	 */
	public MDContainer(ModelConverter modelConverter, Protein protein, double startTemperature, double endTemperature) {
		// set default values
		this.dt = DEFAULT_DT;
		this.numberOfCoolingSteps = DEFAULT_NUMBER_OF_COOLING;
		this.coolingProtocol = DEFAULT_COOLING_PROTOCOL;

		// assign handed over values
		this.data = new HashMap<>();
		this.modelConverter = modelConverter;
		this.protein = protein;
		this.startTemperatur = startTemperature;
		this.endTemperatur = endTemperature;

		// assign some default values - these are SA specific
		this.isobaric = true;
		this.loose = 100.0  * this.dt;
		this.tight = 10.0  * this.dt;
		this.tauTemperature = this.loose;
	}

	public MDData getData(Atom atom) {
		return this.data.get(atom);
	}

	public MDData setData(Atom atom, MDData data) {
		return this.data.put(atom, data);
	}

	public MDData getData(int atomIndex) {
		return getData(findAtomByAtomIndex(atomIndex));
	}

	private Atom findAtomByAtomIndex(int atomIndex) {
		// TODO we have to ensure that the given atomIndex actually corresponds
		// to the natural order of the atom list
		return this.modelConverter.getAtoms(this.protein).get(atomIndex);
	}

	public MDData setData(int atomIndex, MDData data) {
		return this.data.put(findAtomByAtomIndex(atomIndex), data);
	}

	public void setTemperature(double temperature) {
		this.currentTemperature = temperature;		
	}
	
	public double getTemperature() {
		return this.currentTemperature;
	}

	public double getTauTemperature() {
		return tauTemperature;
	}

	public void setTauTemperature(double tauTemperature) {
		this.tauTemperature = tauTemperature;
	}

	public double getTight() {
		return this.tight;
	}

	public double getLoose() {
		return this.loose;
	}
}
