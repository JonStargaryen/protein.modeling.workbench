package de.bioforscher.pmw.reconstruction.ca.mds;

import java.util.List;

import org.osgi.service.log.LogService;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.Residue;
import de.bioforscher.pmw.reconstruction.factory.ReconstructionAlgorithm;

/**
 * A rather basic reconstruction algorithm for the placement of CA atoms by knowing their atomic distances.<br />
 * Finds a configuration of 3D coordinates explaining the given distances as well as possible.
 * 
 * @author S
 *
 */
public class MultiDimensionalScalingReconstruction implements ReconstructionAlgorithm {
	
	private LinearAlgebra linearAlgebra;
	private ModelConverter modelConverter;

	public MultiDimensionalScalingReconstruction(LogService logger, LinearAlgebra linearAlgebra, ModelConverter modelConverter) {
		this.linearAlgebra = linearAlgebra;
		this.modelConverter = modelConverter;
	}
	
	@Override
	public void reconstruct(Protein protein) {
		// compute distance map - TODO: remove this later on, as the distance map will then by provided by the predicted contacts
		@SuppressWarnings("deprecation")
		double[][] distanceMap = new de.bioforscher.pmw.reconstruction.ca.DistanceMapComposer(this.linearAlgebra, this.modelConverter).computeDistanceMap(protein);
		
		// forget previous coordinates
		this.modelConverter.removeAtoms(protein);
		// end remove
		
		MultiDimensionalScaling mds = new MultiDimensionalScaling();
		List<double[]> placedAtoms = mds.computeEmbedding(distanceMap);
		
		List<Residue> residues = this.modelConverter.getResidues(protein);
		if(placedAtoms.size() != residues.size()) {
			throw new IllegalArgumentException("reconstruction and protein size do not match (anymore?)");
		}
		
		for(int i = 0; i < placedAtoms.size(); i++) {
			this.modelConverter.createAtom(residues.get(i), ModelConverter.BACKBONE_CA_NAME, placedAtoms.get(i));
		}
	}
}
