package de.bioforscher.pmw.reconstruction.md.annealing;

import org.osgi.service.log.LogService;

import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.reconstruction.factory.ReconstructionAlgorithm;
import de.bioforscher.pmw.reconstruction.md.MDContainer;
import de.bioforscher.pmw.reconstruction.md.MDContainer.CoolingProtocol;

/**
 * Adapted implementation of a simulated annealing routine for molecular structures.
 * 
 * @see Tinker's anneal.f
 * 
 * @author S
 *
 */
public class SimulatedAnnealing implements ReconstructionAlgorithm {			
	private ModelConverter modelConverter;
	/**
	 * container/wrapper of the {@link Protein}
	 */
	private MDContainer container;

	private static final double DEFAULT_START_TEMPERATURE = 1000;
	private static final double DEFAULT_END_TEMPERATURE = 0;
	
	public SimulatedAnnealing(LogService logger, LinearAlgebra linearAlgebra, ModelConverter modelConverter) {
		this.modelConverter = modelConverter;
	}

	@Override
	public void reconstruct(Protein protein) {
		this.container = new MDContainer(this.modelConverter, protein, DEFAULT_START_TEMPERATURE, DEFAULT_END_TEMPERATURE);
		
		initializeSimulation(this.container);
		
		/**
		 * shakeup() - sets bond lengths and angles to expected values and fixates certain bonds
		 * mdinit() - ensure atoms have non-zero masses
		 */
		/**
		 * equilibration phase
		 */
		//TODO implement someday
		
		/**
		 * cooling phase
		 */
		for(int step = 1; step < this.container.getNumberOfCoolingSteps(); step++) {
			double ratio = computeTemperatureRatio(this.container.getCoolingProtocol(), step, this.container.getNumberOfCoolingSteps());
			this.container.setTemperature(this.container.getStartTemperatur() * (1.0 - ratio) + this.container.getEndTemperatur() * ratio);
			this.container.setTauTemperature(this.container.getLoose() * (1.0 - ratio) + this.container.getTight() * ratio);
			/**
			 * beeman(istep,dt)
			 */
		}
	}

	private void initializeSimulation(MDContainer container2) {
		// assign 1-3, 1-4, 1-5 connections of atoms		
	}

	private double computeTemperatureRatio(CoolingProtocol coolingProtocol, int step, int numberOfCoolingSteps) {
		double ratio = (double) step / (double) numberOfCoolingSteps;
		switch(coolingProtocol) {
		case LINEAR:
			return ratio;
		case SIGMOID:
			return sigmoid(3.5, ratio);
		case EXPONENT:
			return 1.0 - Math.exp(-5.0 * ratio);
		default:
			throw new UnsupportedOperationException("cooling protocol " + coolingProtocol.name() + " is unknown");
		}
	}

	/**
	 * compute the value of the normalized sigmoidal function
	 * "sigmoid" implements a normalized sigmoidal function on the
	 * interval [0,1]; the curves connect (0,0) to (1,1) and have
	 *  a cooperativity controlled by beta, they approach a straight
	 *  line as beta -> 0 and get more nonlinear as beta increases
	 * @param beta
	 * @param x
	 * @return
	 */
	private double sigmoid(double beta, double x) {
		if (beta == 0) {
			return x;
		}
		double expmax = 1.0 / (Math.exp(-beta) + 1.0);
        double expmin = 1.0 / (Math.exp(beta) + 1.0);
        double expterm = 1.0 / (Math.exp(beta*(2.0 * x - 1.0)) + 1.0);
        return (expmax - expterm) / (expmax - expmin);
	}
}
