package de.bioforscher.pmw.api;

import java.util.List;

import de.bioforscher.pmw.model.Alignment;
import de.bioforscher.pmw.model.Atom;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.Residue;

/**
 * High-level access to alignment functions.
 * @author S
 *
 */
public interface AlignmentService {
	/**
	 * Aligns two fragments. Both fragments must have the equal size (with respect to residues, not necessarily atoms). For each {@link Residue} the minimal set of atoms to align will be determined, this shared sets of atoms are what will actually be used to compute the {@link Alignment}.<br />
	 * <b>Important: the atoms/residues given in the arguments will not be manipulated.</b> To actually superimpose the fragments, use {@link LinearAlgebra#transform(double[], double[], double[][])} and the translation and rotation given by {@link Alignment#translationVector} respectively {@link Alignment#rotationMatrix} (for now, this is subject to change once the API structure settles).
	 * @param reference A fragment (collection of residues).
	 * @param fragmentToAlign Another fragment (collection of residues).
	 * @return The {@link Alignment} container describing the yielded alignment in detail.
	 */
	Alignment alignFragments(List<Residue> reference, List<Residue> fragmentToAlign);
	
	void transform(Protein protein, double[] translation, double[][] rotation);
	
	void transform(List<Atom> atoms, double[] translation, double[][] rotation);
}
