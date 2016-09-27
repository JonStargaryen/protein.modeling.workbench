package de.bioforscher.pmw.alignment.provider;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.jama.SingularValueDecomposition;
import org.osgi.service.component.annotations.Activate;
import org.osgi.service.component.annotations.Reference;
import de.bioforscher.pmw.api.AlignmentService;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.model.Alignment;
import de.bioforscher.pmw.model.Atom;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.Residue;

public class AlignmentServiceImplOld implements AlignmentService {
	@Reference
	private LinearAlgebra linearAlgebra;
	@Reference
	private ModelConverter modelConverter;

	@Activate
	public void activate() {
//		new StructureComposition(this.modelConverter, this);
	}

	@Override
	public Alignment alignFragments(List<Residue> reference, List<Residue> fragmentToAlign) {
		// ensure same size
		if (reference.size() != fragmentToAlign.size()) {
			throw new IllegalArgumentException("fragments must have the same number of residues - found "
					+ reference.size() + " and " + fragmentToAlign.size());
		}

		// create sets of shared atom arrays
		List<double[]> atoms1 = getCAsArray(reference);
		List<double[]> atoms2 = getCAsArray(fragmentToAlign);
		
		// calculate centroid
		double[] centroid1 = centroid(atoms1);
		double[] centroid2 = centroid(atoms2);

		// center atoms
		atoms1 = center(atoms1, centroid1);
		atoms2 = center(atoms2, centroid2);

		// create matrices and fill them
		Matrix atomMatrix1 = wrapInMatrix(atoms1);
		Matrix atomMatrix2 = wrapInMatrix(atoms2);

		// compute svd
		SingularValueDecomposition svd = atomMatrix2.transpose().times(atomMatrix1).svd();

		Matrix u = svd.getU();
		Matrix v = svd.getV();
		Matrix v_original = (Matrix) v.clone();
		Matrix u_transposed = u.transpose();

		Matrix rotation = v.times(u_transposed);
		Matrix rotation_transposed = rotation.transpose();

		// check if we have found a reflection
		if (rotation.det() < 0) {
			v = v_original.transpose();
			v.set(2, 0, (0 - v.get(2, 0)));
			v.set(2, 1, (0 - v.get(2, 1)));
			v.set(2, 2, (0 - v.get(2, 2)));

			Matrix nv_transposed = v.transpose();
			rotation = nv_transposed.times(u_transposed);
			rotation_transposed = rotation.transpose();
		}

		double[] cb_tmp = this.linearAlgebra.multiply(centroid2, rotation_transposed.getArray());
		double[] translationVector = this.linearAlgebra.add(centroid1, cb_tmp);
		double[][] rotationMatrix = rotation_transposed.getArray();

		Alignment alignment = Alignment.of(atoms1, atoms2, translationVector, rotationMatrix);
		alignment.rmsd = calculateRMSD(alignment);
		return alignment;
	}

	private Matrix wrapInMatrix(List<double[]> atoms) {
		double[][] matrix = new double[atoms.size()][3];
		for (int index = 0; index < atoms.size(); index++) {
			matrix[index] = atoms.get(index);
		}
		return new Matrix(matrix);
	}

	private List<double[]> center(List<double[]> atoms, final double[] centroid) {
		return atoms.stream().map(a -> this.linearAlgebra.subtract(a, centroid)).collect(Collectors.toList());
	}

	private double[] centroid(List<double[]> atoms) {
		double[] centroid = new double[3];
		centroid[0] = atoms.stream().mapToDouble(d -> d[0]).average().getAsDouble();
		centroid[1] = atoms.stream().mapToDouble(d -> d[1]).average().getAsDouble();
		centroid[2] = atoms.stream().mapToDouble(d -> d[2]).average().getAsDouble();
		return centroid;
	}

	/**
	 * get all Calphas and make them the foundation of the alignment
	 * 
	 * @param alignment
	 */
	private List<double[]> getCAsArray(List<Residue> residues) {
		return residues.stream().map(r -> Arrays.copyOf(this.modelConverter.getCA(r).xyz, 3)).collect(Collectors.toList());
	}

	private double calculateRMSD(Alignment alignment) {
		double rmsd = 0;
		for (int index = 0; index < alignment.atoms1.size(); index++) {
			//TODO these atoms are already translated, when moving this to the API, remind this!
			/**
			 * [-1.1305000000000014, -0.047999999999994714, -2.8354999999999997]
 			-> [1.8304715998598668, 0.13060744395015184, -2.439817599266646]
			 */
//			double[] atom2 = alignment.atoms2.get(index);
//			System.out.println(Arrays.toString(atom2));
//			atom2 = this.linearAlgebra.multiply(atom2, alignment.rotationMatrix);
//			System.out.println(" -> " + Arrays.toString(atom2));
//			rmsd += this.linearAlgebra.distanceFast(alignment.atoms1.get(index), atom2);
			double[] atom1 = alignment.atoms1.get(index);
			System.out.println(Arrays.toString(atom1));
			double[] atom2 = alignment.atoms2.get(index);
			System.out.println(" vs " + Arrays.toString(atom2));
			atom2 = this.linearAlgebra.transform(atom2, new double[] {0,0,0}, alignment.rotationMatrix);
			System.out.println(" vs -> " + Arrays.toString(atom2));
			
			rmsd += this.linearAlgebra.distanceFast(atom1, atom2);
		}
		return Math.sqrt(rmsd / alignment.atoms1.size());
	}

	/**
	 * 
	 * @param residue1
	 *            an atom container
	 * @param residue2
	 *            the other atom container
	 * @return set of atom names occuring in both residues
	 */
	Set<String> determineSharedAtomNames(Residue residue1, Residue residue2) {
		Set<String> residue1AtomNames = residue1.atoms.stream().map(a -> a.name).collect(Collectors.toSet());
		Set<String> residue2AtomNames = residue2.atoms.stream().map(a -> a.name).collect(Collectors.toSet());
		residue1AtomNames.retainAll(residue2AtomNames);
		return residue1AtomNames;
	}

	@Override
	public void transform(Protein protein, double[] translation, double[][] rotation) {
		transform(this.modelConverter.getAtoms(protein), translation, rotation);		
	}

	@Override
	public void transform(List<Atom> atoms, double[] translation, double[][] rotation) {
		atoms.parallelStream().forEach(a -> transform(a, translation, rotation));
	}
	
	private void transform(Atom atom, double[] translation, double[][] rotation) {
		atom.xyz = this.linearAlgebra.transform(atom.xyz, translation, rotation);
	}
}
