package de.bioforscher.pmw.reconstruction.factory;

import org.osgi.service.log.LogService;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.reconstruction.backbone.bbq.BackboneBuildingFromQuadrilaterals;
import de.bioforscher.pmw.reconstruction.ca.mds.MultiDimensionalScalingReconstruction;
import de.bioforscher.pmw.reconstruction.md.annealing.SimulatedAnnealing;
import de.bioforscher.pmw.reconstruction.sidechain.pulchra.Pulchra;

public class CoordinateReconstructionAlgorithmFactory {
	private LogService logger;
	private LinearAlgebra linearAlgebra;
	private ModelConverter modelConverter;
	
	public CoordinateReconstructionAlgorithmFactory(LogService logger, LinearAlgebra linearAlgebra, ModelConverter modelConverter) {
		this.logger = logger;
		this.linearAlgebra = linearAlgebra;
		this.modelConverter = modelConverter;
	}
	
	public ReconstructionAlgorithm createMultiDimensionalScalingReconstruction() {
		return new MultiDimensionalScalingReconstruction(this.logger, this.linearAlgebra, this.modelConverter);
	}
	
	public ReconstructionAlgorithm createBackboneBuildingFromQuadrilaterals() {
		return new BackboneBuildingFromQuadrilaterals(this.logger, this.linearAlgebra, this.modelConverter);
	}
	
	public ReconstructionAlgorithm createPulchra() {
		return new Pulchra(this.logger, this.linearAlgebra, this.modelConverter);		
	}

	public ReconstructionAlgorithm createSimulatedAnnealing() {
		return new SimulatedAnnealing(this.logger, this.linearAlgebra, this.modelConverter);
	}
}
