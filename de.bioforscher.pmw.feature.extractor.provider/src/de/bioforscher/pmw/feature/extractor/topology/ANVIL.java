package de.bioforscher.pmw.feature.extractor.topology;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.osgi.service.log.LogService;
import de.bioforscher.pmw.api.FeatureExtractor;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.feature.extractor.core.AbstractFeatureProvider;
import de.bioforscher.pmw.feature.extractor.core.Annotator;
import de.bioforscher.pmw.model.Atom;
import de.bioforscher.pmw.model.Chain;
import de.bioforscher.pmw.model.FeatureType;
import de.bioforscher.pmw.model.Membrane;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.Residue;
import de.bioforscher.pmw.model.Topology;

/**
 * An adaptation of the ANVIL algorithm which processes protein structures and places the membrane in the most suitable way.<br /><br />
 * 
 * Intellectual and implementation credit to:<br />
 * <b>[Postic, 2015] - Membrane positioning for high- and low-resolution protein structures through a binary classification approach - Guillaume Postic, Yassine Ghouzam, Vincent Guiraud, and Jean-Christophe Gelly, Protein Engineering, Design & Selection, 2015, 1–5, doi: 10.1093/protein/gzv063</b><br /><br />
 * 
 * ANVIL original documentation:<br />
 * <pre>© Univ. Paris Diderot & Inserm, 2015

guillaume.postic@univ-paris-diderot.fr

This software is a computer program whose purpose is to assign membrane
boundaries to a protein three-dimensional structure, by using the spatial
coordinates of the alpha carbons. It can also process coarse-grained
models of protein structures. The algorithm follows an approach that
treats the problem of membrane assignment as a binary classification.
The software is implemented in Python and requires the PyPy interpreter.
It also requires Naccess 2.1.1 (S.J. Hubbard, 1996) to calculate the
atomic accessible surface.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.</pre>
 */
public class ANVIL extends AbstractFeatureProvider implements Annotator {
	/**
	 * default values
	 */
	public static final List<String> MEMBRANE_AMINO_ACIDS = Arrays.asList("PHE", "MET", "ILE", "LEU", "VAL", "TRP", "ALA", "CYS", "GLY", "SER", "HIS");
	public static final List<String> SOLVENT_AMINO_ACIDS = Arrays.asList("ARG", "ASP", "LYS", "GLU", "ASN", "GLN", "PRO", "THR", "TYR");
	private static final int DEFAULT_NUMBER_OF_SPHERE_POINTS = 350;
	private static final double DEFAULT_STEP = 1.0;
	private static final double DEFAULT_MINTHICK = 20.0;
	private static final double DEFAULT_MAXTHICK = 40.0;
	private static final double DEFAULT_AFILTER = 40.0;
	private static final double DEFAULT_DENSITY_OF_MEMBRANE_POINTS = 2.0;
	
	/**
	 * parameters
	 */
	private double step;
	private double minthick;
	private double maxthick;	
	private double afilter;
	private int numberOfSpherePoints;
	private double density;

	/**
	 * stuff computed by the algorithm
	 */
	private Protein protein;
	private double[] centerOfMass;
	private double maximalExtent;
//	private FeatureType asaIndex;
	private int hphobtotal;
	private int hphiltotal;
	private PotentialMembrane membrane;
	
	public ANVIL(FeatureExtractor featureExtractor, LogService logger, LinearAlgebra linearAlgebra, ModelConverter modelConverter) {
		// provides membrane topology information, depends on ASA annotation
    	super(featureExtractor, logger, linearAlgebra, modelConverter, new FeatureType[] { FeatureType.MEMBRANE_TOPOLOGY }, FeatureType.ACCESSIBLE_SURFACE_AREA);
    }

	@Override
	protected void computeFeatureInternal(Protein protein) {
		//TODO: some day: expose options
		this.numberOfSpherePoints = DEFAULT_NUMBER_OF_SPHERE_POINTS;
		this.afilter = DEFAULT_AFILTER;
		this.minthick = DEFAULT_MINTHICK;
		this.maxthick = DEFAULT_MAXTHICK;
		this.step = DEFAULT_STEP;
		this.density = DEFAULT_DENSITY_OF_MEMBRANE_POINTS;
		
		this.protein = protein;
//		this.asaIndex = FeatureType.ACCESSIBLE_SURFACE_AREA;
//		this.asaIndex = this.modelConverter.determineFeatureIndex(protein, FeatureType.ACCESSIBLE_SURFACE_AREA);
		this.centerOfMass = centerOfMass(this.protein);
		// extent maximal extent some more
		this.maximalExtent = 1.2 * maximalExtent(protein, this.centerOfMass);
		int[] hphobHphil = hphobHphil();
		this.hphobtotal = hphobHphil[0];
		this.hphiltotal = hphobHphil[1];
		
		// determine the best possible membrane placement among all generated spherePoints
		this.logger.log(LogService.LOG_DEBUG, "initial membrane placement");
		PotentialMembrane initialMembrane = processSpherePoints(generateSpherePoints(this.numberOfSpherePoints));
//		System.out.println("initial quality: " + initialMembrane.qmax);
		
		// try to improve the results by finding a better axis relative to the currently winning membrane
		this.logger.log(LogService.LOG_DEBUG, "trying to refine");
		PotentialMembrane alternativeMembrane = processSpherePoints(findProximateAxes(initialMembrane));
//		System.out.println("refined quality: " + alternativeMembrane.qmax);
		
		this.logger.log(LogService.LOG_DEBUG, "membrane inclination did " + (initialMembrane.qmax > alternativeMembrane.qmax ? "not " : "") + "improve");
		this.membrane = initialMembrane.qmax > alternativeMembrane.qmax ? initialMembrane : alternativeMembrane;
		
		this.logger.log(LogService.LOG_DEBUG, "adjusting thickness");
		this.step = 0.3;
		double thickness = this.linearAlgebra.distance(this.membrane.c1, this.membrane.c2);
		this.logger.log(LogService.LOG_DEBUG, "membrane thickness is " + thickness + " A");
		//TODO: implement adaptation of membrane thickness
		
		assignTopology();
		placeMembraneMolecules(this.protein);
	}

	private double[] centerOfMass(Protein protein2) {
		int atomCount = 0;
		double[] coordinates = { 0, 0, 0 };
		for(Chain c : protein.chains) {
			for(Residue g : c.residues) {
				for(Atom a : g.atoms) {
					if(a.name.equals(ModelConverter.BACKBONE_CA_NAME)) {
						coordinates[0] += a.xyz[0];
						coordinates[1] += a.xyz[1];
						coordinates[2] += a.xyz[2];
						atomCount++;
					}
				}
			}
		}
		return new double[] {
			coordinates[0] / atomCount,
			coordinates[1] / atomCount,
			coordinates[2] / atomCount
		};
	}
	
	/**
	 * tests whether a point is in the corridor described by the given membrane
	 * TODO: streamline this cascade of call once we know if other classes will depend on this functionality
	 * @param pointToTest
	 * @param membrane
	 * @return true iff the point is embedded in the membrane
	 */
	private boolean isInSpace(Residue residue, Membrane membrane) {
		return isInSpace(this.modelConverter.getCA(residue).xyz, membrane);
	}
	
	private boolean isInSpace(double[] pointToTest, Membrane membrane) {
		return isInSpace(pointToTest, membrane.normalVector, membrane.planePoint1, membrane.planePoint2);
	}

	/**
	 * returns true iff the point is in between two planes defined by the vector diamVector and each passes through p1 respectively p2<br /><br />
	 * defined by ANVIL
	 * @return
	 */
	private boolean isInSpace(double[] pointToTest, double[] normalVector, double[] planePoint1, double[] planePoint2) {
		normalVector = this.linearAlgebra.normalize(normalVector);
		final double d1 = - this.linearAlgebra.dotProduct(normalVector, planePoint1);
		final double d2 = - this.linearAlgebra.dotProduct(normalVector, planePoint2);
		final double d = - this.linearAlgebra.dotProduct(normalVector, pointToTest);
		return d > Math.min(d1, d2) && d < Math.max(d1, d2);
	}

	/**
	 * computes the maximal extent of this protein in any given spatial direction to the center of mass of this structure
	 * @param protein
	 * @return the maximal distance occuring between the center of mass and any other atom
	 */
	private double maximalExtent(Protein protein, double[] centerOfMass) {
		return protein.chains.stream().flatMap(c ->
			c.residues.stream()).mapToDouble(r ->
			this.linearAlgebra.distance(this.modelConverter.getCA(r).xyz, centerOfMass)).max().getAsDouble();
	}
	
	/**
	 * computes the distance of the proposed membrane molecule to the protein - this can be used to ensure that no membrane molecules are placed within the protein
	 * @param atom
	 * @param protein
	 * @return the minimal squared distance of this atom to any CA atom of the protein
	 */
	private double minimalSquaredDistanceToProteinCAAtom(double[] atom, Protein protein) {
		//TODO this is wrong - either no membrane is added or no molecules too close to the protein are discarded
		return protein.chains.stream().flatMap(c -> c.residues.stream()).mapToDouble(r ->
			this.linearAlgebra.distanceFast(this.modelConverter.getCA(r).xyz, this.centerOfMass)).min().getAsDouble();
	}
	
	private void assignTopology() {
		Membrane membrane = new Membrane();
		membrane.centerOfMass = this.centerOfMass;
//		membrane.setMcc(this.membrane.qmax);
		membrane.normalVector = this.linearAlgebra.subtract(this.membrane.c1, this.membrane.c2);
		membrane.planePoint1 = this.membrane.c1;
		membrane.planePoint2 = this.membrane.c2;
		membrane.spherePoint = this.membrane.point;
		this.protein.membrane = membrane;
		for(Chain c : this.protein.chains) {
			for(Residue r : c.residues) {
//				System.out.println(r + " is positioned in the membrane: " + isInSpace);
				double value = isInSpace(r, membrane) ? Topology.TRANSMEMBRANE.ordinal() : Topology.NON_TRANSMEMBRANE.ordinal();
				r.features.put(FeatureType.MEMBRANE_TOPOLOGY.name(), wrapInArray(value));
			}
		}
	}

	/**
	 * places evenly distributed pseudo-atoms for sake of membrane visualization
	 * @param protein
	 */
	private void placeMembraneMolecules(Protein protein) {
		double radius = this.maximalExtent * this.maximalExtent;
		Membrane membrane = protein.membrane;
		double[] normalVector = membrane.normalVector;
		for(double[] layer : Arrays.asList(membrane.planePoint1, membrane.planePoint2)) {
			double d = - this.linearAlgebra.dotProduct(normalVector, layer);
			for(double i = -1000; i < 1000; i += this.density) {
				for(double j = -1000; j < 1000; j += this.density) {
					double[] atom = new double[] { i, j, 0 };
					try {
						atom[2] = -(d + i * normalVector[0] + j * normalVector[1]) / normalVector[2];
					} catch (ArithmeticException e) {
					}
					
					// distance cutoff is also squared
					if(this.linearAlgebra.distanceFast(atom, layer) <= radius && minimalSquaredDistanceToProteinCAAtom(atom, protein) > 12.0) {
						membrane.membraneMolecules.add(atom);
					}
				}
			}
		}
	}

	/**
	 * find the best possible solution among all proposed spherePoints
	 * @param spherePoints axes to check
	 * @return the best possible membrane, embedding as many residues as possible
	 */
	private PotentialMembrane processSpherePoints(List<double[]> spherePoints) {
		// best performing membrane
		PotentialMembrane membrane = null;
		// best performing membrane's score
		double qmax = 0;
		
		// construct slices of thickness 1.0 along the axis connecting the centerOfMass and the spherePoint
		for(int spIndex = 0; spIndex < spherePoints.size(); spIndex++) {
//			if(spIndex % 100 == 0) {
//				System.out.println(spIndex + " / " + spherePoints.size());
//			}
			double[] spherePoint = spherePoints.get(spIndex);
			double[] diam = this.linearAlgebra.multiply(this.linearAlgebra.subtract(this.centerOfMass, spherePoint), 2.0);
			double diamNorm = this.linearAlgebra.norm(diam);
			
			List<PotentialMembrane> qvartemp = new ArrayList<>();
			
			for(double i = 0; i < diamNorm - this.step; i += this.step) {
				double dPointC1 = i;
				double dPointC2 = i + this.step;
				
				double[] c1 = thales(diam, dPointC1, spherePoint);
				double[] c2 = thales(diam, dPointC2, spherePoint);
				
				// evaluate how well this membrane slice embeddeds the peculiar residues
				int[] hphobHphil = hphobHphil(true, diam, c1, c2);
				
				qvartemp.add(new PotentialMembrane(c1, c2, hphobHphil));
			}
			
			int jmax = (int) ((this.minthick / this.step) - 1);
	
			for(double width = 0; width < this.maxthick; width = (jmax + 1) * this.step) {
//				System.out.println(width + " / " + this.maxthick);
	            int imax = qvartemp.size() - 1 - jmax;
	            
	            for(int i = 0; i < imax; i++) {
	            	double[] c1 = qvartemp.get(i).c1;
	            	double[] c2 = qvartemp.get(i + jmax).c2;
//	            	System.out.println("distance between points is " + distance(c1, c2));
	           
	            	double hphob = 0;
	            	double hphil = 0;
	            	double total = 0;
	            	
	            	for(int j = 0; j < jmax; j++) {
	            		PotentialMembrane ij = qvartemp.get(i + j);
	            		if(j == 0 || j == jmax - 1) {
	            			hphob += 0.5 * ij.hphob;
	            			hphil += 0.5 * ij.hphil;
	            		} else {
	            			hphob += ij.hphob;
	            			hphil += ij.hphil;
	            		}
	            		total += ij.total;
	            	}
	            	
	            	if(hphob > 0) {
	            		double qvaltest = qValue(hphil, hphob, this.hphiltotal, this.hphobtotal);
//	            		System.out.println("membrane " + c1 + " -> " + c2 + ": " + qvaltest);
	            		if(qvaltest > qmax) {
	            			qmax = qvaltest;
	            			membrane = new PotentialMembrane(spherePoint, c1, c2, hphob, hphil, total, qmax);
	            		}
	            	}
	            }
	            jmax++;
			}
		}
		
		return membrane;
	}

	/**
	 * creates a number of axes close to that of the given membrane
	 * @param membrane a well-positioned, but potentially not optimal, membrane
	 * @return 350 axes close to that of the membrane
	 */
	private List<double[]> findProximateAxes(PotentialMembrane membrane) {
		List<double[]> points = generateSpherePoints(30000);
		Collections.sort(points, new Comparator<double[]>() {			
			@Override
			public int compare(double[] d1, double[] d2) {
				return Double.compare(distance(d1), distance(d2));
			}
			
			private double distance(double[] d) {
				return ANVIL.this.linearAlgebra.distance(d, membrane.point);
			}
		});
		return points.subList(0, this.numberOfSpherePoints);
	}

	/**
	 * counts the total number of exposed hydrophobic respectively hydrophilic residues
	 * @return
	 */
	private int[] hphobHphil() {
		// delegate the more fine-grained impl when not interested in the placement of a residue relative
		// to the potential membrane plane
		return hphobHphil(false, null, null, null);
	}
	
	/**
	 * counts how well hydrophobic and hydrophilic residues exposed to the solvent are embedded by the membrane
	 * @param checkMembranePlane when false residues only need to be exposed in order to count (this is used for the initial 'global' counting of hydrophobic/hydrophilic residues within the structure)
	 * @param diam parameters describing the membrane placement
	 * @param c1 parameters describing the membrane placement
	 * @param c2 parameters describing the membrane placement
	 * @return [countOfHydrophobicResidues, countOfHydrophilicResidues]
	 */
	private int[] hphobHphil(boolean checkMembranePlane, double[] diam, double[] c1, double[] c2) {
		int[] hphobHphil = { 0, 0 };

	    for(Chain chain : this.protein.chains) {
	    	for(Residue residue : chain.residues) {
	    		// skip residues with too low ASA values - in the original code this is 
				// checked after determining whether the residue is within the membrane 
				// but this check should be way faster and, thus, reduce computational load
	    		if(residue.features.get(FeatureType.ACCESSIBLE_SURFACE_AREA.name())[0] < this.afilter) {
	    			continue;
	    		}
	    		
	    		// give the option to ignore the membrane placement
	    		if(checkMembranePlane && !isInSpace(this.modelConverter.getCA(residue).xyz, diam, c1, c2)) {
	    			continue;
	    		}
	    		
	    		if(MEMBRANE_AMINO_ACIDS.contains(residue.aminoAcid)) {
	    			hphobHphil[0]++;
	    		} else {
	    			hphobHphil[1]++;
	    		}
	    	}
	    }
	    
	    return hphobHphil;
	}

	/**
	 * evaluates the quality of any given membrane slice
	 * @param hphil
	 * @param hphob
	 * @param hphiltotal
	 * @param hphobtotal
	 * @return the qValue - the higher, the better is the embedding described by this membrane
	 */
	private double qValue(double hphil, double hphob, double hphiltotal, double hphobtotal) {
		if(hphobtotal < 1) {			
			hphobtotal=0.1;
		}
	    if(hphiltotal < 1) {
	        hphiltotal += 1;
	    }
	    double partTotal = hphob + hphil;

	    return (hphob * (hphiltotal - hphil) - hphil * (hphobtotal - hphob)) /
	        (partTotal * hphobtotal * hphiltotal * (hphobtotal + hphiltotal - partTotal));
	}

	/**
	 * return a point S so that S = V * d + P 
	 * @param vector V
	 * @param distance d
	 * @param point P
	 * @return a point
	 */
	private double[] thales(double[] vector, double distance, double[] point) {
		return this.linearAlgebra.add(point, this.linearAlgebra.multiply(vector, distance / this.linearAlgebra.norm(vector)));
	}
	
	/**
	 * generates a defined number of points on a sphere with radiues <code>maximalExtent</code> around <code>centerOfMass</code>
	 * @return
	 */
	private List<double[]> generateSpherePoints(int numberOfSpherePoints) {
		List<double[]> points = new ArrayList<>();
		
		double oldPhi = 0;
		for(int k = 1; k < numberOfSpherePoints + 1; k++) {
			double h = -1 + 2*(k-1)/((double) (numberOfSpherePoints-1));
			double theta = Math.acos(h);
			double phi = (k == 1 || k == numberOfSpherePoints) ? 0 : (oldPhi + 3.6/Math.sqrt(numberOfSpherePoints*(1-h*h))) % (2*Math.PI);
			
			double[] point = new double[] {
				this.maximalExtent * Math.sin(phi) * Math.sin(theta) + this.centerOfMass[0],
				this.maximalExtent * Math.cos(theta) + this.centerOfMass[1],
				this.maximalExtent * Math.cos(phi) * Math.sin(theta) + this.centerOfMass[2]
			};
			points.add(point);
			oldPhi = phi;
		}
		
		return points;
	}
}

class PotentialMembrane {
	double[] point;
	double[] c1;
	double[] c2;
	double hphob;
	double hphil;
	double total;
	double qmax;
	
	public PotentialMembrane(double[] c1, double[] c2, int[] hphobHphil) {
		this.c1 = c1;
		this.c2 = c2;
		this.hphob = hphobHphil[0];
		this.hphil = hphobHphil[1];
		this.total = this.hphob + this.hphil;
//		System.out.println("constructing potential membrane with\nc1: " + Arrays.toString(c1) + "\nc2: " + Arrays.toString(c2) + "\nphob: " + this.hphob + "\nphil: " + this.hphil);
	}

//	public double[] computeNormalVector() {
//		return this.linearAlgebra.substract(this.c1, this.c2);
//	}

	public PotentialMembrane(double[] spherePoint, double[] c1, double[] c2, double hphob, double hphil,
			double total, double qmax) {
		this.point = spherePoint;
		this.c1 = c1;
		this.c2 = c2;
		this.hphob = hphob;
		this.hphil = hphil;
		this.total = total;
		this.qmax = qmax;
	}
}