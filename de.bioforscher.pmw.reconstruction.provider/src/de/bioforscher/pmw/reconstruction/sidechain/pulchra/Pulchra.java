package de.bioforscher.pmw.reconstruction.sidechain.pulchra;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.slf4j.Logger;

import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.Residue;
import de.bioforscher.pmw.reconstruction.factory.ReconstructionAlgorithm;

/**
 * Reduced implementation of the Pulchra algorithm.<br />
 * Intellectual and implementation credit to:<br />
 * <b>[Rotkiewicz, 2008] - Fast procedure for reconstruction of full-atom protein models from reduced representations, 10.1002/jcc.20906</b><br /><br />
 * 
 * Places side chains by finding suitable rotamers from a library.<br /><br />
 * 
 * A rather rudimentary adaptation of the PULCHRA algorithm for side chain
 * placement.<br />
 * Actually, several other features are available in the native impl such as
 * backbone placement or several refinement steps. Long term, they may be added
 * to this implementation. For now only plain side chain reconstruction is
 * supported.<br />
 * <br />
 * 
 * original comment:
 * 
 * <pre>
 * PULCHRA Protein Chain Restoration Algorithm
 * 
 * Version 3.04 July 2007 Contact: Piotr Rotkiewicz, piotr -at- pirx -dot- com
 * </pre>
 * 
 * @see J Comput Chem. 2008 July 15; 29(9): 1460–1465. doi:10.1002/jcc.20906.
 */
public class Pulchra implements ReconstructionAlgorithm {
	private LinearAlgebra linearAlgebra;
	private ModelConverter modelConverter;
	
	private static final double BIN_SIZE = 0.3;
	private static final String BASE_PATH = "de/bioforscher/pmw/reconstruction/sidechain/pulchra/";
	private static final String ROT_STAT_IDX_LIBRARY = BASE_PATH + "rot_data_idx.h";
	private static final String ROT_STAT_COORDS_LIBRARY = BASE_PATH + "rot_data_coords.h";
	@SuppressWarnings("unused")
	private static final String ROT_STAT_LIBRARY = BASE_PATH + "nco_data.h";

	private static Map<String, String[]> sideChainAtomNames;
	private static List<int[]> rotStatIdx;
	private static List<double[]> rotStatCoords;
//	private static Map<int[], double[]> ncoStat;
//	private static Map<int[], double[]> ncoStatPro;
	
	private List<Residue> residues;

	public Pulchra(Logger logger, LinearAlgebra linearAlgebra, ModelConverter modelConverter) {
		if(rotStatIdx == null) {
			initializeLibrary();
		}
		
		this.linearAlgebra = linearAlgebra;
		this.modelConverter = modelConverter;
	}

	@Override
	public void reconstruct(Protein protein) {
//		System.out.println("pulchra sanity check: " + rotStatIdx.size() + " - " + rotStatCoords.size());
		this.residues = this.modelConverter.getResidues(protein);
		buildSidechains();
	}

	private void buildSidechains() {
		// TODO again reconstruct missing leading/tailing residues
		for (int residueIndex = 2; residueIndex < this.residues.size() - 1; residueIndex++) {			
			Residue r_p2 = this.residues.get(residueIndex - 2);
			Residue r_p = this.residues.get(residueIndex - 1);
			Residue r_c = this.residues.get(residueIndex);
			if(r_c.aminoAcid.equals("G") /*|| r_c.getAminoAcid().equals('A')*/ || r_c.aminoAcid.equals("X")) {
				// nothing to reconstruct for glycine
				//TODO: for alanine, the cb is already placed by BBQ - longterm, this may not the case for other backbone reconstruction algorithms
				//TODO: some fallback for X?
				continue;
			}
			Residue r_n = this.residues.get(residueIndex + 1);
			
			String aminoAcidName = r_c.aminoAcid;

			double[] ca_p2 = this.modelConverter.getCA(r_p2).xyz;
			double[] ca_p = this.modelConverter.getCA(r_p).xyz;
			double[] ca_c = this.modelConverter.getCA(r_c).xyz;
			double[] ca_n = this.modelConverter.getCA(r_n).xyz;
			
			double d13 = this.linearAlgebra.distance(ca_p2, ca_c);
			double d24 = this.linearAlgebra.distance(ca_p, ca_n);
			double d14 = this.linearAlgebra.distance14(ca_p2, ca_p, ca_c, ca_n);
			int[] residueBinning = binResidues(d13, d24, d14);
			int bin13_1 = residueBinning[0];
			int bin13_2 = residueBinning[1];
			int bin14 = residueBinning[2];
			
			// find closest rotamer conformation
			int[] bestMatchingRotamer = null;
			double bestMatchingRotamerDistance = Double.MAX_VALUE;
			int aminoAcidIndex = getAminoAcidIndex(r_c);
			for(int[] rotamer : rotStatIdx) {
				// check whether rotamer describes the correct amino acid
				if(rotamer[0] != aminoAcidIndex) {
					continue;
				}
				
				double rotamerDistance = Math.abs(rotamer[1] - bin13_1) +
						Math.abs(rotamer[2] - bin13_2) + 
						0.2 * Math.abs(rotamer[3] - bin14);
				if(rotamerDistance < bestMatchingRotamerDistance) {
					bestMatchingRotamerDistance = rotamerDistance;
					bestMatchingRotamer = rotamer;
				}
			}
			
			// new rebuild
			double[] translation = ca_c;
			double[][] rotation = this.linearAlgebra.rotation(ca_p, ca_c, ca_n);
			
			int pos = bestMatchingRotamer[5];
			int nsc = getHeavySideChainAtomCount(r_c);
			// all atoms within the coordinate file describing this residue
			for(int j = 0; j < nsc; j++) {
//				System.out.println(j + " : " + aminoAcidName + " : " + determineAtomName(aminoAcidName, j) + " : " + Arrays.toString(rotStatCoords.get(pos + j + 1)));
				this.modelConverter.createAtom(r_c, determineAtomName(aminoAcidName, j), this.linearAlgebra.rototranslate(rotStatCoords.get(pos + j + 1), translation, rotation));
			}
		}
	}
	
	public static int getAminoAcidIndex(Residue residue) {
		//TODO: move
		final String aminoAcids = "GASCVTIPMDNLKEQRHFYWX";
		return aminoAcids.indexOf(residue.aminoAcid);
	}
	
	public static int getHeavySideChainAtomCount(Residue residue) {
		//TODO: move
		final int[] heavySideChainAtomCount = { 0, 1, 2, 2, 3, 3, 4, 3, 4, 4, 4, 4, 5, 5, 5, 7, 6, 7, 8, 10, 0 };
		return heavySideChainAtomCount[getAminoAcidIndex(residue)];
	}
	
	private int[] binResidues(double d13_1, double d13_2, double d14) {
		int bin13_1 = (int) ((d13_1 - 4.6) / BIN_SIZE);
		int bin13_2 = (int) ((d13_2 - 4.6) / BIN_SIZE);
		int bin14 = (int) ((d14 + 11.0) / BIN_SIZE);
		
		bin13_1 = capToInterval(bin13_1, 0, 9);
		bin13_2 = capToInterval(bin13_2, 0, 9);
		bin14 = capToInterval(bin14, 0, 73);
		
		return new int[] { bin13_1, bin13_2, bin14 };
	}
	
	/**
	 * cap a variable x by a lower bound x_i and x_a
	 * @param value
	 * @param lowerBound
	 * @param upperBound
	 * @return the original value if it lies in the interval [x_i,x_a], else x_i or x_a
	 */
	private int capToInterval(int value, int lowerBound, int upperBound) {
		List<Integer> values = Arrays.asList(value, lowerBound, upperBound);
		Collections.sort(values);
		return values.get(1);
	}
	
	// TODO: move
	private String determineAtomName(String aminoAcidName, int j) {
		return sideChainAtomNames.get(aminoAcidName)[j];
	}
	
	private synchronized void initializeLibrary() {
		// parse sidechain indices
		rotStatIdx = new ArrayList<>();
		InputStream idxIs = Thread.currentThread().getContextClassLoader().getResourceAsStream(ROT_STAT_IDX_LIBRARY);
		new BufferedReader(new InputStreamReader(idxIs)).lines()
			// filter lines which do not contain information
			.filter(line -> line.contains("},")).forEach(line -> {
			// remove padding
			String[] tmp = line.replace("{", "").replace("}", "").split(",");
			rotStatIdx.add(parseRotStatIdxLine(tmp));
		});
		
		// parse side chain coordinates
		rotStatCoords = new ArrayList<>();
		InputStream coordsIs = Thread.currentThread().getContextClassLoader().getResourceAsStream(ROT_STAT_COORDS_LIBRARY);
		new BufferedReader(new InputStreamReader(coordsIs)).lines()
			// filter lines which do not contain information
			.filter(line -> line.endsWith(",")).forEach(line -> {
			String[] tmp = line.replace("{", "").replace("}", "").split(",");
			rotStatCoords.add(parseRotStatCoordsLine(tmp));
		});
		
		// parse backbone library for proline
//		ncoStat = new HashMap<>();
		//TODO: impl
		
		// parse backbone lirbary for non-prolines
//		ncoStatPro = new HashMap<>();
		//TODO: impl
		
		// TODO: move this some where else
		sideChainAtomNames = new HashMap<>();
		sideChainAtomNames.put("G", new String[] {});
		sideChainAtomNames.put("A", new String[] { "CB" });
		sideChainAtomNames.put("S", new String[] { "CB", "OG" });
		sideChainAtomNames.put("C", new String[] { "CB", "SG" });
		sideChainAtomNames.put("V", new String[] { "CB", "CG1", "CG2" });
		sideChainAtomNames.put("T", new String[] { "CB", "OG1", "CG2" });
		sideChainAtomNames.put("I", new String[] { "CB", "CG1", "CG2", "CD1" });
		sideChainAtomNames.put("P", new String[] { "CB", "CG", "CD" });
		sideChainAtomNames.put("M", new String[] { "CB", "CG", "SD", "CE" });
		sideChainAtomNames.put("D", new String[] { "CB", "CG", "OD1", "OD2" });
		sideChainAtomNames.put("N", new String[] { "CB", "CG", "OD1", "ND2" });
		sideChainAtomNames.put("L", new String[] { "CB", "CG", "CD1", "CD2" });
		sideChainAtomNames.put("K", new String[] { "CB", "CG", "CD", "CE", "NZ" });
		sideChainAtomNames.put("E", new String[] { "CB", "CG", "CD", "OE1", "OE2" });
		sideChainAtomNames.put("Q", new String[] { "CB", "CG", "CD", "OE1", "NE2" });
		sideChainAtomNames.put("R", new String[] { "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2" });
		sideChainAtomNames.put("H", new String[] { "CB", "CG", "ND1", "CD2", "CE1", "NE2" });
		sideChainAtomNames.put("F", new String[] { "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ" });
		sideChainAtomNames.put("Y", new String[] { "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH" });
		sideChainAtomNames.put("W", new String[] { "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2" });
	}
	
	private double[] parseRotStatCoordsLine(String[] tmp) {
		return new double[] {
			Double.valueOf(tmp[0].trim()),
			Double.valueOf(tmp[1].trim()),
			Double.valueOf(tmp[2].trim())
		};
	}
	
	private int[] parseRotStatIdxLine(String[] tmp) {
		return new int[] {
			Integer.valueOf(tmp[0].trim()),
			Integer.valueOf(tmp[1].trim()),
			Integer.valueOf(tmp[2].trim()),
			Integer.valueOf(tmp[3].trim()),
			Integer.valueOf(tmp[4].trim()),
			Integer.valueOf(tmp[5].trim())
		};
	}
}
