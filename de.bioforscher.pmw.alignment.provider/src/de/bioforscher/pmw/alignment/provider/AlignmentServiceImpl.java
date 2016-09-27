package de.bioforscher.pmw.alignment.provider;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.jama.SingularValueDecomposition;
import org.osgi.service.component.annotations.Activate;
import org.osgi.service.component.annotations.Component;
import org.osgi.service.component.annotations.Reference;

import de.bioforscher.pmw.api.AlignmentService;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.model.Alignment;
import de.bioforscher.pmw.model.Atom;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.Residue;

@Component(name = "de.bioforscher.pmw.alignment")
public class AlignmentServiceImpl implements AlignmentService {
	@Reference
	private LinearAlgebra linearAlgebra;
	@Reference
	private ModelConverter modelConverter;
	/**
	 * the atom names used to represent the individual amino acids
	 */
	private static final String[] ATOM_NAMES = { "CA" };
	
	@Activate
	public void activate() {
		new StructureComposition(modelConverter, this);
	}
	
	@Override
	public Alignment alignFragments(final List<Residue> reference, final List<Residue> fragmentToAlign) {
		List<double[]> referenceAtoms = extractAtoms(reference, ATOM_NAMES);
		List<double[]> fragmentAtoms = extractAtoms(fragmentToAlign, ATOM_NAMES);
		
		double[] referenceCentroid = calculateCentroid(referenceAtoms);
		double[] fragmentCentroid = calculateCentroid(fragmentAtoms);
//		System.out.println("old centroid1: " + Arrays.toString(referenceCentroid));
//		System.out.println("old centroid2: " + Arrays.toString(fragmentCentroid));
		// center atoms
		List<double[]> centeredReferenceAtoms = calculateCenteredAtoms(referenceAtoms, referenceCentroid);
		List<double[]> centeredFragmentAtoms = calculateCenteredAtoms(fragmentAtoms, fragmentCentroid);
//		System.out.println("new centroid1: " + Arrays.toString(calculateCentroid(centeredReferenceAtoms)));
//		System.out.println("new centroid2: " + Arrays.toString(calculateCentroid(centeredFragmentAtoms)));
		// wrap in matrices
		Matrix referenceMatrix = wrapInMatrix(centeredReferenceAtoms);
		Matrix fragmentMatrix = wrapInMatrix(centeredFragmentAtoms);
		
		// compose covariance matrix and calculate SVD
		SingularValueDecomposition svd = fragmentMatrix.transpose().times(referenceMatrix).svd();
		// R = (V * U')'
		Matrix ut = svd.getU().transpose();
		Matrix rotationMatrix = svd.getV().times(ut).transpose();
		// check if reflection
		if(rotationMatrix.det() < 0) {
			Matrix v = svd.getV().transpose();
			v.set(2, 0, (0 - v.get(2, 0)));
			v.set(2, 1, (0 - v.get(2, 1)));
			v.set(2, 2, (0 - v.get(2, 2)));
			rotationMatrix = v.transpose().times(ut).transpose();
		}
		double[][] rotation = rotationMatrix.getArray();
		
		// compute translation
		Matrix referenceCentroidMatrix = new Matrix(referenceCentroid, 1);
		Matrix fragmentCentroidMatrix = new Matrix(fragmentCentroid, 1);
		double[] translation = referenceCentroidMatrix.minus(fragmentCentroidMatrix.times(rotationMatrix)).getRowPackedCopy();
		
		// compute rmsd and prepare return type
		Alignment alignment = Alignment.of(referenceAtoms, fragmentAtoms, translation, rotation);
		alignment.rmsd = calculateRMSD(referenceAtoms, fragmentAtoms, translation, rotation);
		
		return alignment;
	}
	
	private Matrix wrapInMatrix(List<double[]> atoms) {
		double[][] matrix = new double[atoms.size()][3];
		for (int index = 0; index < atoms.size(); index++) {
			matrix[index] = atoms.get(index);
		}
		return new Matrix(matrix);
	}
	
	private double[] calculateCentroid(List<double[]> atoms) {
		return new double[] { atoms.stream().mapToDouble(d -> d[0]).average().getAsDouble(),
			atoms.stream().mapToDouble(d -> d[1]).average().getAsDouble(),
			atoms.stream().mapToDouble(d -> d[2]).average().getAsDouble() };
	}
	
	private List<double[]> calculateCenteredAtoms(final List<double[]> atoms, final double[] centroid) {
		return atoms.stream().map(a -> this.linearAlgebra.subtract(a, centroid)).collect(Collectors.toList());
	}
	
	/**
	 * Computes the root-mean-square deviation for 2 atom sets. <b>Important:</b> do not use already aligned atoms here, but the initial ones.
	 * @param atoms1
	 * @param atoms2
	 * @param translation
	 * @param rotation
	 * @return
	 */
	private double calculateRMSD(final List<double[]> atoms1, final List<double[]> atoms2, final double[] translation, final double[][] rotation) {
		double rmsd = 0;
		for (int index = 0; index < atoms1.size(); index++) {
			double[] atom1 = atoms1.get(index);
			double[] atom2 = this.linearAlgebra.transform(atoms2.get(index), translation, rotation);
			
			rmsd += this.linearAlgebra.distanceFast(atom1, atom2);
		}
		return Math.sqrt(rmsd / atoms1.size());
	}
	
	/**
	 * Extracts all atoms matching one of the given names.
	 * @param residues a collection of residues
	 * @param atomNames names to retain according to PDB naming - e.g. 'CA'
	 * @return all atoms of the given name
	 */
	private List<double[]> extractAtoms(final List<Residue> residues, final String... atomNames) {
		final List<String> atomNamesList = Arrays.asList(atomNames);
		return residues.stream().flatMap(r -> r.atoms.stream()).filter(a -> atomNamesList.contains(a.name)).map(a -> a.xyz).collect(Collectors.toList());
	}
	
	@Override
	public void transform(final Protein protein, final double[] translation, final double[][] rotation) {
		transform(this.modelConverter.getAtoms(protein), translation, rotation);
	}
	
	@Override
	public void transform(final List<Atom> atoms, final double[] translation, final double[][] rotation) {
		atoms.forEach(a -> {
//			System.out.println("moved " + Arrays.toString(a.xyz));
			a.xyz = this.linearAlgebra.transform(a.xyz, translation, rotation);
//			System.out.println(" to " + Arrays.toString(a.xyz));
		});
	}
}
