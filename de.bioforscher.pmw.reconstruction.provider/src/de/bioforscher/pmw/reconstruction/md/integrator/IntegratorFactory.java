package de.bioforscher.pmw.reconstruction.md.integrator;

import de.bioforscher.pmw.api.LinearAlgebra;

public class IntegratorFactory {
	private final LinearAlgebra linearAlgebra;
	
	public IntegratorFactory(LinearAlgebra linearAlgebra) {
		this.linearAlgebra = linearAlgebra;
	}

	public Integrator createBeemanIntegrator() {
		return new BeemanIntegrator(this.linearAlgebra);
	}
}
