package de.bioforscher.pmw.reconstruction.backbone.bbq;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.osgi.service.log.LogService;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.model.Atom;
import de.bioforscher.pmw.model.Chain;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.Residue;
import de.bioforscher.pmw.reconstruction.factory.ReconstructionAlgorithm;

/**
 * Reduced implementation of the BBQ algorithm.<br />
 * Intellectual and implementation credit to:<br />
 * <b>[Gront, 2007] - Backbone Building from Quadrilaterals: A Fast and Accurate Algorithm for Protein Backbone Reconstruction from Alpha Carbon Coordinates</b><br /><br />
 *
 * A algorithm which can reconstruct backbone atoms from a CA trace.<br />
 * This is done by finding suitable known backbone conformations (stored in a library), whereby <i>suitable</i> is determined by a hash function based on the distances of 4 consecutive residues.
 */
public class BackboneBuildingFromQuadrilaterals implements ReconstructionAlgorithm {
	private LinearAlgebra linearAlgebra;
	private ModelConverter modelConverter;
	private LogService logger;
	
	private static final int QUADRILATERAL_SETOFF = 3;
	private static final int[] BINS = { 1, 32, 1024, 33, 31, 1025, 1023, 1056, 992, 1057, 1055, 993, 991, 2, 128, 2048, 130, 126, 2050, 2046, 2176, 1920, 2178, 2174, 1922, 1918 };
	private static final String LIBRARY = "de/bioforscher/pmw/reconstruction/backbone/bbq/bbq-library.dat";
	private static final double BIN_SIZE = 0.2;
	private static final int MAX_INDEX = 107221;
	private static Map<Integer, List<double[]>> quadrilateralLookup;
	
	public BackboneBuildingFromQuadrilaterals(LogService logger, LinearAlgebra linearAlgebra, ModelConverter modelConverter) {
		if(quadrilateralLookup == null) {
			initializeLibrary();
		}
		
		this.linearAlgebra = linearAlgebra;
		this.modelConverter = modelConverter;
		this.logger = logger;
	}
	
	private Residue createProxyResidue(List<Residue> residues) {
		Residue proxyResidue = new Residue();
		Atom proxyAtom = this.modelConverter.createAtom(ModelConverter.BACKBONE_CA_NAME, this.linearAlgebra.subtract(residues.get(0).atoms.get(0).xyz, residues.get(1).atoms.get(0).xyz));
		proxyResidue.atoms.add(proxyAtom);
		return proxyResidue;
	}
	
	private int determineFallbackIndex(int index) {	
		int[] arrayOfInt = getBinNeighborhood(index);
	    for (int j = 0; j < arrayOfInt.length; j++) {
	      if ((arrayOfInt[j] >= 0) && 
	        (arrayOfInt[j] < MAX_INDEX)) {
	        if (quadrilateralLookup.get(arrayOfInt[j]) != null) {
	          return arrayOfInt[j];
	        }
	      }
	    }
		
	    int bmu = -1;
	    double bmuDistance = Double.MAX_VALUE;

	    for(int i : quadrilateralLookup.keySet()) {
	    	double distance = distanceToBin(index, i);
//	    	System.out.println(distance + " for " + i + " vs " + index);
	    	if(distance < bmuDistance) {
	    		bmuDistance = distance;
	    		bmu = i;
	    	}
	    }
//	    System.out.println("computing fallback: bmu index is " + bmu);
	    return bmu;
	}
	
	private int[] getBinNeighborhood(int index) {
		int[] arrayOfInt = new int[BINS.length << 1];
	    for (int i = 0; i < BINS.length; i++) {
	      arrayOfInt[(2 * i)] = (index + BINS[i]);
	      arrayOfInt[(2 * i + 1)] = (index - BINS[i]);
	    }
//	    System.out.println("bin neighborhood: " + Arrays.toString(arrayOfInt));
	    return arrayOfInt;
	}
	
	private synchronized void initializeLibrary() {
		quadrilateralLookup =  new HashMap<>();
		InputStream is = Thread.currentThread().getContextClassLoader().getResourceAsStream(LIBRARY);
		new BufferedReader(new InputStreamReader(is)).lines().forEach(line -> {
			String[] tmp = line.trim().split("\\s+");
			quadrilateralLookup.put(Integer.parseInt(tmp[0]), parseAtoms(tmp));
		});
	}
	
	private int lookupSuitableQuadrilateral(double[] ca1, double[] ca2, double[] ca3, double[] ca4) {
		double d13 = this.linearAlgebra.distance(ca1, ca3);
	    double d14 = this.linearAlgebra.distance14(ca1, ca2, ca3, ca4);
	    double d24 = this.linearAlgebra.distance(ca2, ca4);
	    
	    return lookupSuitableQuadrilateral(d13, d24, d14);
	}
	
	private double distanceToBin(int index1, int index2) {
		int i = index1 - index2 & 0x1F;
		int j = (index1 >> 5) - (index2 >> 5) & 0x1F;
		index1 = (index1 >> 10) - (index2 >> 10);

		return Math.sqrt(i * i + j * j + index1 * index1);
	}
	
	private int lookupSuitableQuadrilateral(double d13, double d24, double d14) {
	    int index = 0;
		int i = (int) ((d13 - 3) / BIN_SIZE);
	    int j = (int) ((d24 - 3) / BIN_SIZE);
	    int k;
	    if (d14 > 0) {
	      k = (int) ((d14 - 3) / BIN_SIZE);
	    } else {
	      k = (int) ((-d14 - 3) / BIN_SIZE);
	      index = 65536;
	    }
	    index += i;
	    if (k > 0) {
	      index += (j << 5);
	      index += (k << 10);
	    }
	    return index;
	}
	
	/**
	 * 
	 * @param r1 the quadrilateral
	 * @param r2
	 * @param r3
	 * @param r4
	 * @return the most suitable quadrilateral (by convention the atoms are in order CA, O, N, CB)
	 */
	private List<double[]> lookupSuitableQuadrilateral(Residue r1, Residue r2, Residue r3, Residue r4) {
		double[] ca1 = r1.atoms.get(0).xyz;
		double[] ca2 = r2.atoms.get(0).xyz;
		double[] ca3 = r3.atoms.get(0).xyz;
		double[] ca4 = r4.atoms.get(0).xyz;
		int index = lookupSuitableQuadrilateral(ca1, ca2, ca3, ca4);

		List<double[]> sc = quadrilateralLookup.getOrDefault(index, quadrilateralLookup.get(determineFallbackIndex(index)));
//		if(sc == null) {
//			index = determineFallbackIndex(index);
//			sc = quadrilateralLookup.get(index);
//		}
		
		return sc;
	}

	private List<double[]> parseAtoms(String[] strArray) {
		List<double[]> atoms = new ArrayList<>();
		// C
		atoms.add(new double[] { Double.parseDouble(strArray[2]), Double.parseDouble(strArray[3]), Double.parseDouble(strArray[4]) });
		// O
		atoms.add(new double[] { Double.parseDouble(strArray[6]), Double.parseDouble(strArray[7]), Double.parseDouble(strArray[8]) });
		// N
		atoms.add(new double[] { Double.parseDouble(strArray[10]), Double.parseDouble(strArray[11]), Double.parseDouble(strArray[12]) });
		// CB
		atoms.add(new double[] { Double.parseDouble(strArray[14]), Double.parseDouble(strArray[15]), Double.parseDouble(strArray[16]) });
		return atoms;
	}
	
	@Override
	public void reconstruct(Protein protein) {
		this.logger.log(LogService.LOG_DEBUG, "bbq sanity check: " + quadrilateralLookup.size());
		for(Chain chain : protein.chains) {
			List<Residue> residues = chain.residues;
			if(residues.size() < 4) {
				System.err.println("cannot handle chains with less than 4 residues");
				continue;
			}
			// reconstruct first quadrilaterals
			reconstructFragment(createProxyResidue(residues), residues.get(0), residues.get(1), residues.get(2));
			//TODO: reconstruct first residue and all its backbone atoms 'truly'
			
			// traverse all possible quadrilaterals
			for(int residueIndex = 1; residueIndex < chain.residues.size() - QUADRILATERAL_SETOFF + 1; residueIndex++) {
				reconstructFragment(residues.get(residueIndex - 1), residues.get(residueIndex), residues.get(residueIndex + 1), residues.get(residueIndex + 2));
			}
			
			//TODO: same for the last residues
		}
	}
	
	private void reconstructFragment(Residue residue1, Residue residue2, Residue residue3, Residue residue4) {
		List<double[]> fragmentScaffold = lookupSuitableQuadrilateral(residue1, residue2, residue3, residue4);
		double[] c = fragmentScaffold.get(0);
		double[] o = fragmentScaffold.get(1);
		double[] n = fragmentScaffold.get(2);
		
		// compute rototranslation to the local coordinate system
		double[] translation = residue2.atoms.get(0).xyz;
		double[][] rotation = this.linearAlgebra.rotation(this.modelConverter.getCA(residue1).xyz,
				this.modelConverter.getCA(residue2).xyz,
				this.modelConverter.getCA(residue3).xyz);
		
		// create and add atom
		this.modelConverter.createAtom(residue1, "C", this.linearAlgebra.rototranslate(c, translation, rotation));
		this.modelConverter.createAtom(residue1, "O", this.linearAlgebra.rototranslate(o, translation, rotation));
		this.modelConverter.createAtom(residue2, "N", this.linearAlgebra.rototranslate(n, translation, rotation));

		// cannot perform these options for glycines without CB
		//TODO: trust BBQ or PULCHRA more?
//		if(!residue2.getAminoAcid().equals('G')) {
//			double[] cb = fragmentScaffold.get(3);
//			CoordinateManipulations.placeAtom(residue2, "CB", rototranslate(cb, translation, rotation));
//		}
	}
}
