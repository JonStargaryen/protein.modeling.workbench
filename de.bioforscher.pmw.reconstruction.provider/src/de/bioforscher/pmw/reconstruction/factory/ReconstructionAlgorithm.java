package de.bioforscher.pmw.reconstruction.factory;

import de.bioforscher.pmw.model.Protein;

/**
 * The contract which all reconstruction algorithms fulfill. They differ in the
 * level of reconstruction they are able to provide: some can convert distance
 * maps to CA coordianates, some place side chains after the backbone was
 * reconstructed. Obviously, no guarantees are made about the concrete
 * implementation. Many methods exist, but are not native Java code, thus, some
 * existing approaches were ported to Java and are available.
 * 
 * @author S
 *
 */
public interface ReconstructionAlgorithm {
	/**
	 * reconstruct a protein as good as possible according the rules of the impl
	 * 
	 * @param protein
	 *            the object to be reconstructed
	 */
	void reconstruct(Protein protein);
}
