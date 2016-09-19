package de.bioforscher.pmw.reconstruction.ca;

import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.model.Chain;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.Residue;

@Deprecated
public class DistanceMapComposer {
	private LinearAlgebra linearAlgebra;
	private ModelConverter modelConverter;

	public DistanceMapComposer(LinearAlgebra linearAlgebra, ModelConverter modelConverter) {
		this.linearAlgebra = linearAlgebra;
		this.modelConverter = modelConverter;
	}
	
	public double[][] computeDistanceMap(Protein protein) {		
		double[][] distanceMap = new double[protein.size][protein.size];
		for(Chain chain1 : protein.chains) {
			for(Residue residue1 : chain1.residues) {
				for(Chain chain2 : protein.chains) {
					for(Residue residue2 : chain2.residues) {
						double distance = this.linearAlgebra.distance(this.modelConverter.getCA(residue1).xyz,
								this.modelConverter.getCA(residue2).xyz);
						int i = residue1.residueId;
						int j = residue2.residueId;
						distanceMap[i][j] = distance;
						distanceMap[j][i] = distance;
					}
				}
			}
		}
		return distanceMap;
	}
}
