package de.bioforscher.pmw.feature.extractor.asa;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.slf4j.Logger;

import de.bioforscher.pmw.api.FeatureExtractor;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.feature.extractor.core.AbstractFeatureProvider;
import de.bioforscher.pmw.feature.extractor.core.Annotator;
import de.bioforscher.pmw.feature.extractor.core.Element;
import de.bioforscher.pmw.model.Atom;
import de.bioforscher.pmw.model.FeatureType;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.Residue;

/**
 * Computes the accessible surface area of each residue in a {@link Protein}.<br />
 * original doc:<br /><br />
 * 
 * Class to calculate Accessible Surface Areas based on
 * the rolling ball algorithm by Shrake and Rupley.
 *
 * The code is adapted from a python implementation at http://boscoh.com/protein/asapy
 * (now source is available at https://github.com/boscoh/asa).
 * Thanks to Bosco K. Ho for a great piece of code and for his fantastic blog.
 *
 * See
 * Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent of Protein Atoms.
 * Lysozyme and Insulin." JMB (1973) 79:351-371.
 * Lee, B., and Richards, F.M. "The interpretation of Protein Structures: Estimation of
 * Static Accessibility" JMB (1971) 55:379-400
 * @author duarte_j
 *
 */
public class AccessibleSurfaceAreaCalculator extends AbstractFeatureProvider implements Annotator {

	public AccessibleSurfaceAreaCalculator(FeatureExtractor featureExtractor, Logger logger, LinearAlgebra linearAlgebra, ModelConverter modelConverter) {
		super(featureExtractor, logger, linearAlgebra, modelConverter, new FeatureType[] { FeatureType.ACCESSIBLE_SURFACE_AREA });
	}
	
	// Bosco uses as default 960, Shrake and Rupley seem to use in their paper 92 (not sure if this is actually the same parameter)
	public static final int DEFAULT_N_SPHERE_POINTS = 960;
	public static final double DEFAULT_PROBE_SIZE = 1.4;

	// Chothia's amino acid atoms vdw radii
	public static final double TRIGONAL_CARBON_VDW = 1.76;
	public static final double TETRAHEDRAL_CARBON_VDW = 1.87;
	public static final double TRIGONAL_NITROGEN_VDW = 1.65;
	public static final double TETRAHEDRAL_NITROGEN_VDW = 1.50;
	public static final double SULFUR_VDW = 1.85;
	public static final double OXIGEN_VDW = 1.40;

	private double probe;
	
	private List<double[]> spherePoints;
	private double cons;
	private List<Atom> atoms;
	private List<Residue> residues;
	/**
	 * mapping between atoms (identified by their pdbSerial) and their atom radius (which depends on covalent bounds to the particular atom)
	 */
	private Map<Integer, Double> atomRadii;
	
	@Override
	protected void computeFeatureInternal(Protein protein) {
		this.probe = DEFAULT_PROBE_SIZE;

		this.residues = this.modelConverter.getResidues(protein);
		this.atoms = new ArrayList<>();
		this.atomRadii = new HashMap<>();
		// initialising the radii by looking them up through AtomRadii
		for(Residue residue : this.residues) {
			for(Atom atom : residue.atoms) {
				// skip hydrogen
				if(atom.element.equals("H") || atom.element.equals("D")) {
					continue;
				}
				
				this.atomRadii.put(atom.pdbSerial, determineRadius(residue, atom));
				// we add them explicitly here to ensure they are non-hydrogen atoms
				this.atoms.add(atom);
			}
		}
		
		// initialising the sphere points to sample
		this.spherePoints = generateSpherePoints(DEFAULT_N_SPHERE_POINTS);
		this.cons = 4.0 * Math.PI / DEFAULT_N_SPHERE_POINTS;
		
		this.residues.parallelStream().forEach(r -> {
			r.features.put(FeatureType.ACCESSIBLE_SURFACE_AREA.name(), wrapInArray(r.atoms.stream().mapToDouble(this::calcSingleAsa).sum()));
		});
	}

	/**
	 * Returns list of 3d coordinates of points on a sphere using the
	 * Golden Section Spiral algorithm.
	 * @param nSpherePoints the number of points to be used in generating the spherical dot-density
	 * @return
	 */
	private List<double[]> generateSpherePoints(int nSpherePoints) {
		List<double[]> points = new ArrayList<>();
		double inc = Math.PI * (3.0 - Math.sqrt(5.0));
		double offset = 2.0 / nSpherePoints;
		for (int k = 0 ; k < nSpherePoints; k++) {
			double y = k * offset - 1.0 + (offset / 2.0);
			double r = Math.sqrt(1.0 - y * y);
			double phi = k * inc;
			points.add(new double[] { Math.cos(phi) * r, y, Math.sin(phi) * r });
		}
		return points;
	}

	 /** Gets the van der Waals radius of the given atom following the values defined by
	 * Chothia (1976) J.Mol.Biol.105,1-14
	 * NOTE: the vdw values defined by the paper assume no Hydrogens and thus "inflates"
	 * slightly the heavy atoms to account for Hydrogens. Thus this method cannot be used
	 * in a structure that contains Hydrogens!
	 */
	private double determineRadius(Residue residue, Atom atom) {
		switch(atom.element) {
		case "H": case "D":
			return Element.H.getVDWRadius();
		case "O":
			return OXIGEN_VDW;
		case "S":
			return SULFUR_VDW;
		case "N":
			return atom.name.equals("NZ") ? TETRAHEDRAL_NITROGEN_VDW : TRIGONAL_NITROGEN_VDW;
		case "C":
			String atomName = atom.name;
			if(atomName.equals("C") || atomName.equals("CE1") || atomName.equals("CE2") || atomName.equals("CE3") ||
					atomName.equals("CH2") || atomName.equals("CZ") || atomName.equals("CZ2") || atomName.equals("CZ3")) {
				return TRIGONAL_CARBON_VDW;
			}
			if (atomName.equals("CA") || atomName.equals("CB") || atomName.equals("CE") || atomName.equals("CG1") || atomName.equals("CG2")) {
				return TETRAHEDRAL_CARBON_VDW;
			}
			switch(this.modelConverter.convertToOneLetterCode(residue.aminoAcid)) {
			case "F": case "W": case "Y": case "H": case "D": case "N":
				return TRIGONAL_CARBON_VDW;
			case "P": case "K": case "R": case "M": case "I": case "L":
				return TETRAHEDRAL_CARBON_VDW;
			case "Q": case "E":
				return atomName.equals("CD") ? TRIGONAL_CARBON_VDW : TETRAHEDRAL_CARBON_VDW;
			default:
				throw new IllegalArgumentException("unknown case for residue: " + residue);
			}
		default:
			throw new IllegalArgumentException("unknown case for atom: " + atom);
		}
	}
	

	/**
	 * Returns list of indices of atoms within probe distance to atom k.
	 * @param k index of atom for which we want neighbor indices
	 */
	private List<Atom> findNeighbors(Atom atom) {
		List<Atom> neighborAtoms = new ArrayList<>();
		double radius = this.atomRadii.get(atom.pdbSerial) + this.probe + this.probe;
//		for(Atom potentialNeighbor : this.atoms) {
		for(int i = 0; i < this.atoms.size(); i++) {
			Atom potentialNeighbor = this.atoms.get(i);
			if(potentialNeighbor.equals(atom)) {
				continue;
			}
			
			double distance = this.linearAlgebra.distance(potentialNeighbor.xyz, atom.xyz);
			if(distance < radius + this.atomRadii.get(potentialNeighbor.pdbSerial)) {
				neighborAtoms.add(potentialNeighbor);
			}
		}
		return neighborAtoms;
	}

	private double calcSingleAsa(Atom atom) {
		List<Atom> neighborAtoms = findNeighbors(atom);
		double radius = this.probe + this.atomRadii.get(atom.pdbSerial);
		int accessiblePoints = 0;

//		for (double[] point : this.spherePoints) {
		for(int i = 0; i < this.spherePoints.size(); i++) {
			double[] point = this.spherePoints.get(i);
			boolean isAccessible = true;
			double[] testPoint = this.linearAlgebra.add(this.linearAlgebra.multiply(point, radius), atom.xyz);
			for(Atom neighborAtom : neighborAtoms) {
				double r = this.atomRadii.get(neighborAtom.pdbSerial) + this.probe;
				double differenceSquared = this.linearAlgebra.distanceFast(testPoint, neighborAtom.xyz);
				if (differenceSquared < r * r) {
					isAccessible = false;
					break;
				}
			}
			if (isAccessible) {
				accessiblePoints++;
			}
		}

		return this.cons * accessiblePoints * radius * radius;
	}
}
