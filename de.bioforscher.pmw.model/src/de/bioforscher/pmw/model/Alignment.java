package de.bioforscher.pmw.model;

import java.util.List;

import org.osgi.dto.DTO;

/**
 * Describes the results of an alignment of two populations of {@link Residue} collections.
 * @author S
 *
 */
public class Alignment extends DTO {
//	/**
//	 * The instance of one residue list used to create this alignment. This is the reference of the alignment. The other residue list will be superimposed onto this.
//	 */
//	public List<Residue> residues1;
//	/**
//	 * The other instance of residues used as input. The translation vector and rotation matrix describe the rototranslation necessary to superimpose this to the reference residue list.
//	 */
//	public List<Residue> residues2;
	/**
	 * The aligned atoms of {@link Alignment#residues1}. Atoms are not part of these atom arrays, when they are not part of both fragments. They are identified by their atom name. The ordering in both atom arrays is identical.
	 */
	public List<double[]> atoms1;
	/**
	 * The aligned atoms of {@link Alignment#residues2}.
	 */
	public List<double[]> atoms2;
	/**
	 * The root-mean-square deviation (RMSD) of this alignment.
	 */
	public double rmsd;
	/**
	 * The translation vector necessary to move {@link Alignment#residues2} to the position of {@link Alignment#residues1}.
	 */
	public double[] translationVector;
	/**
	 * The rotation matrix necessary to rotate {@link Alignment#residues2} to the position of {@link Alignment#residues1}.
	 */
	public double[][] rotationMatrix;
	
	public static Alignment of(List<double[]> reference, List<double[]> fragmentToAlign, double[] translationVector, double[][] rotationMatrix) {
		Alignment alignment = new Alignment();
		alignment.atoms1 = reference;
		alignment.atoms2 = fragmentToAlign;
		alignment.translationVector = translationVector;
		alignment.rotationMatrix = rotationMatrix;
		return alignment;
	}
}