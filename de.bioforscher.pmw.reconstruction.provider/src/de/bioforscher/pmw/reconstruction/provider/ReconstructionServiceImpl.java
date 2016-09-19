package de.bioforscher.pmw.reconstruction.provider;

import org.osgi.service.component.annotations.Activate;
import org.osgi.service.component.annotations.Component;
import org.osgi.service.component.annotations.Reference;
import org.osgi.service.log.LogService;

import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.api.ReconstructionService;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.ReconstructionLevel;
import de.bioforscher.pmw.reconstruction.factory.CoordinateReconstructionAlgorithmFactory;
import de.bioforscher.pmw.reconstruction.factory.ReconstructionAlgorithm;

/**
 * 
 */
@Component(name = "de.bioforscher.pmw.reconstruction")
public class ReconstructionServiceImpl implements ReconstructionService {
	@Reference
	private LinearAlgebra linearAlgebra;
	@Reference
	private ModelConverter modelConverter;
	private CoordinateReconstructionAlgorithmFactory factory;
	
	private LogService log;
	
	@Reference
	public void setLogService(LogService log) {
		this.log = log;
	}
	
	@Activate
	public void activate() {
		this.factory = new CoordinateReconstructionAlgorithmFactory(this.log, this.linearAlgebra, this.modelConverter);
	}
	
	@Override
	public void reconstruct(Protein protein, ReconstructionLevel reconstructionLevel) {
		// if nothing is to reconstruct: return
		if(reconstructionLevel.equals(protein.reconstructionLevel)) {
			throw new IllegalArgumentException("warning: no point in reconstructing a structure with reconstruction level '" +
					protein.reconstructionLevel + "' to " + 
					reconstructionLevel + " - request a level higher than that already present");
		}
		
		// more than 1 steps needed
		if(reconstructionLevel.ordinal() - protein.reconstructionLevel.ordinal() > 1) {
			// recursively delegate to lower reconstruction levels
			reconstruct(protein, ReconstructionLevel.values()[reconstructionLevel.ordinal() - 1]);
		}
		
		ReconstructionAlgorithm reconstructionAlgorithm;
		// set to true when atoms are added to the structure or atom indices are swapped
		boolean requiresReassignment = false;
		
		switch (reconstructionLevel) {
		case CALPHA:
			reconstructionAlgorithm = this.factory.createMultiDimensionalScalingReconstruction();
			requiresReassignment = true;
			break;
		case BACKBONE:
			reconstructionAlgorithm = this.factory.createBackboneBuildingFromQuadrilaterals();
			requiresReassignment = true;
			break;
		case SIDECHAIN:
			reconstructionAlgorithm = this.factory.createPulchra();
			requiresReassignment = true;
			break;
		case REFINED:
			reconstructionAlgorithm = this.factory.createSimulatedAnnealing();
			break;
		default:
			throw new UnsupportedOperationException("no reconstruction routine for " + reconstructionLevel);
		}
		
		int previousAtomCount = this.modelConverter.getAtoms(protein).size();
		this.log.log(LogService.LOG_DEBUG, "employing " + reconstructionAlgorithm.getClass().getSimpleName() + " to reconstruct " + reconstructionLevel.name());
		reconstructionAlgorithm.reconstruct(protein);
		this.log.log(LogService.LOG_DEBUG, "reconstructed protein from " + previousAtomCount + " to " + this.modelConverter.getAtoms(protein).size() + " atoms");

		if(requiresReassignment) {
			// reassign serials implies rearranging the atoms
			this.modelConverter.updatePdbSerials(protein);
		}
		
		// update pdb representation
		this.modelConverter.updatePdbRepresentation(protein);
		
		// assign reconstruction level
		protein.reconstructionLevel = reconstructionLevel;
	}
	
	@Override
	public void reconstruct(Protein protein) {
		this.reconstruct(protein, ReconstructionLevel.REFINED);
	}
}
