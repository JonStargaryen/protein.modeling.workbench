package de.bioforscher.pmw.reconstruction.md.integrator;

import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.reconstruction.md.MDContainer;

/**
 * The default integrator of TINKER (and most other MD frameworks).
 * @author S
 *
 */
@SuppressWarnings(value="unused")
public class BeemanIntegrator implements Integrator {
	/**
	 * reference to the {@link LinearAlgebra} workhorse
	 */
	private LinearAlgebra linearAlgebra;
	/*
	 * some locally used variables
	 */
	private double factor;
	private double dt_x;
	private double part1;
	private double part2;

	public BeemanIntegrator(LinearAlgebra linearAlgebra) {
		this.linearAlgebra = linearAlgebra;
	}

	@Override
	public void integrate(MDContainer container) {
		this.factor = 8;
		this.dt_x = container.getDt() / this.factor;
		this.part1 = 0.5 * this.factor + 1.0;
		this.part2 = this.part1 - 2.0;
	}
}
