package de.bioforscher.pmw.api;

import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.ReconstructionLevel;

/**
 * Provides methods to place spatial coordinates of a protein from several previously extracted features such as contact maps or
 * coarsely placed coordinates which can then be further refined.
 * @author S
 *
 */
public interface ReconstructionService {
	/**
	 * Reconstructs a protein: for the currently available information the 3D structure of this protein is predicted to the requested level.
	 * Some levels may come with internal requirements and will throw exceptions, when these requirements 
	 * are not met and this service has no means of acquiring the needed information.
	 * The reconstruction levels have a defined order (provided by the value of the enum {@link ReconstructionLevel}).
	 * Requesting a lower or identical level to the current one of the structure will not trigger any computation.
	 * When the requested level is increasing by one, only that reconstruction procedure is performed.
	 * Otherwise - when multiple steps should be computed at once - each necessary routine is subsequently invoked.
	 * @param protein the container of all currently available information on this structure
	 * @param reconstructionLevel the target reconstruction level
	 */
	void reconstruct(Protein protein, ReconstructionLevel reconstructionLevel);
	
	/**
	 * standard reconstruction routine from predicted contacts to an {@link ReconstructionLevel#REFINED} structure, using default options along the way
	 * @param protein the container of all currently available information on this structure
	 */
	void reconstruct(Protein protein);
}
