package de.bioforscher.pmw.model.converter.provider;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.UUID;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.osgi.service.component.annotations.Activate;
import org.osgi.service.component.annotations.Component;
import org.osgi.service.component.annotations.Reference;
import org.osgi.service.log.LogService;

import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.model.Atom;
import de.bioforscher.pmw.model.Chain;
import de.bioforscher.pmw.model.Project;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.ReconstructionLevel;
import de.bioforscher.pmw.model.Residue;
import de.bioforscher.pmw.model.converter.parser.PDBConverter;
import de.bioforscher.pmw.model.converter.parser.SimplePDBConverter;

/**
 * A version of the converter not depended on BioJava.
 * 
 * @author S
 *
 */
@Component(name = "de.bioforscher.pmw.model.converter")
public class ModelConverterImpl implements ModelConverter {
	/**
	 * can be used to sort atoms in a residue container to match the ordering in a PDB file
	 */
	private static final Comparator<Atom> ATOM_NAME_COMPARATOR = new Comparator<Atom>() {
		final String[] ORDER = new String[] {
				BACKBONE_N_NAME,
				BACKBONE_CA_NAME,
				BACKBONE_C_NAME,
				BACKBONE_O_NAME,
				SIDECHAIN_CB_NAME
			};
		
		@Override
		public int compare(Atom a1, Atom a2) {
			return orderOf(a1) - orderOf(a2);
		}
		
		private int orderOf(Atom a) {
			for(int i = 0; i < ORDER.length; i++) {
				if(ORDER[i].equals(a.name)) {
					return i;
				}
			}
			return ORDER.length;
		}
	};

	/**
	 * connection between 1-letter and 3-letter amino acid name codes
	 */
	private Map<String, String> aminoAcidNameMapping;
	@Reference
	private LogService logger;
	
	@Activate
	public void activate() {
		this.aminoAcidNameMapping = new HashMap<>();
		this.aminoAcidNameMapping.put("ALA", "A");
		this.aminoAcidNameMapping.put("ARG", "R");
		this.aminoAcidNameMapping.put("ASN", "N");
		this.aminoAcidNameMapping.put("ASP", "D");
		this.aminoAcidNameMapping.put("CYS", "C");
		this.aminoAcidNameMapping.put("GLN", "Q");
		this.aminoAcidNameMapping.put("GLU", "E");
		this.aminoAcidNameMapping.put("GLY", "G");
		this.aminoAcidNameMapping.put("HIS", "H");
		this.aminoAcidNameMapping.put("ILE", "I");
		this.aminoAcidNameMapping.put("LEU", "L");
		this.aminoAcidNameMapping.put("LYS", "K");
		this.aminoAcidNameMapping.put("MET", "M");
		this.aminoAcidNameMapping.put("PHE", "F");
		this.aminoAcidNameMapping.put("PRO", "P");
		this.aminoAcidNameMapping.put("SER", "S");
		this.aminoAcidNameMapping.put("THR", "T");
		this.aminoAcidNameMapping.put("TRP", "W");
		this.aminoAcidNameMapping.put("TYR", "Y");
		this.aminoAcidNameMapping.put("VAL", "V");
		this.aminoAcidNameMapping.put(UNKNOWN_AMINO_ACID_THREE_LETTER_CODE, UNKNOWN_AMINO_ACID_ONE_LETTER_CODE);
	}
	
	private Stream<Atom> asAtomStream(Residue residue) {
		return residue.atoms.stream();
	}

	private Stream<Residue> asResidueStream(Chain chain) {
		return chain.residues.stream();
	}

	private void assignResidueIds(Protein protein) {
		int residueId = 0;
		for (Chain chain : protein.chains) {
			for (Residue residue : chain.residues) {
				residue.residueId = residueId;
				residueId++;
			}
		}
	}
	
	@Override
	public String convertToOneLetterCode(String threeLetterCode) {
		return this.aminoAcidNameMapping.getOrDefault(threeLetterCode, UNKNOWN_AMINO_ACID_ONE_LETTER_CODE);
	}
	
	@Override
	public String convertToThreeLetterCode(String oneLetterCode) {
		//TODO it would be nice to get rid of this lambda and replace it with a method reference
		return this.aminoAcidNameMapping.entrySet().stream().filter(e -> e.getValue().equals(oneLetterCode)).findFirst().get().getKey();
	}

	@Override
	public void createAtom(Residue residue, String name, double[] coordinates) {
		Atom atom = createAtom(name, coordinates);
		// ensure no old atom describing is still present
		removeAtomByName(residue, name);
		residue.atoms.add(atom);
	}

	@Override
	public Atom createAtom(String name, double[] xyz) {
		Atom a = new Atom();
		a.element = name.substring(0, 1);
		a.name = name;
		a.xyz = xyz;
		a.occupancy = ModelConverter.DEFAULT_OCCUPANCY;
		a.tempFactor = ModelConverter.DEFAULT_TEMP_FACTOR;
		return a;
	}
	
	@Override
	public Project createModelingProject(byte[] pdbFileContent) throws IOException {
		PDBConverter parser = new SimplePDBConverter();
		Protein protein = parser.parsePDBFile(new ByteArrayInputStream(pdbFileContent));
		updatePdbRepresentation(protein);
		// // FIXME: for test purposes
		// // protein.setReconstructionLevel(ReconstructionLevel.NONE);
		protein.size = getResidues(protein).size();
		protein.reconstructionLevel = ReconstructionLevel.VALIDATED;
		assignResidueIds(protein);
		return createModelingProject(protein,
				getResidues(protein).stream().map(r -> String.valueOf(r.aminoAcid)).collect(Collectors.joining()));
	}

	private Project createModelingProject(Protein protein, String sequence) {
		Project project = new Project();
		project.date = new Date().getTime();
		project._id = UUID.randomUUID().toString();
		project.proteins.add(protein);
//		project.sequence = sequence;
		return project;
	}
	
	public Protein createProteinByPDBId(String pdbId) throws IOException {
		this.logger.log(LogService.LOG_INFO, "downloading PDB structure with id '" + pdbId + "'");
		PDBConverter parser = new SimplePDBConverter();
		Protein protein = parser.parsePDBFile(new URL(String.format(PDB_FETCH_URL, pdbId)).openStream());
		updatePdbRepresentation(protein);
		protein.size = getResidues(protein).size();
		protein.reconstructionLevel = ReconstructionLevel.VALIDATED;
		assignResidueIds(protein);
		return protein;
	}
	
	public Protein createProtein(File file) throws IOException {
		Protein protein = new SimplePDBConverter().parsePDBFile(file);
		updatePdbRepresentation(protein);
		protein.size = getResidues(protein).size();
		protein.reconstructionLevel = ReconstructionLevel.VALIDATED;
		assignResidueIds(protein);
		if(protein.name == null) {
			protein.name = file.getName().split("\\.")[0];
		}
		return protein;
	}

//	@Override
//	public Project createModelingProject(String proteinName, String proteinTitle, String sequence) {
//		Protein protein = new Protein();
//		protein.name = proteinName;
//		protein.title = proteinTitle != null ? proteinTitle : DEFAULT_PROTEIN_TITLE;
//		Chain chain = new Chain();
//		chain.chainId = DEFAULT_CHAIN_ID;
//		for (int resNum = 0; resNum < sequence.length(); resNum++) {
//			Residue residue = new Residue();
////			residue.aminoAcid = String.valueOf(sequence.charAt(resNum));
//			residue.aminoAcid = convertToThreeLetterCode(String.valueOf(sequence.charAt(resNum)));
//			residue.residueNumber = resNum + 1;
//			chain.residues.add(residue);
//		}
//		protein.chains.add(chain);
//		protein.size = sequence.length();
//		protein.reconstructionLevel = ReconstructionLevel.NONE;
//		assignResidueIds(protein);
//		return createModelingProject(protein, sequence);
//	}
	
	@Override
	public Project createModelingProject(String projectName, String sequence) {
		Protein protein = new Protein();
		protein.name = "general information"; // others should be named 'model 1,2,...,n'
		Chain chain = new Chain();
		chain.chainId = DEFAULT_CHAIN_ID;
		for (int resNum = 0; resNum < sequence.length(); resNum++) {
			Residue residue = new Residue();
			residue.aminoAcid = convertToThreeLetterCode(String.valueOf(sequence.charAt(resNum)));
			residue.residueNumber = resNum + 1;
			residue.residueId = resNum;
			chain.residues.add(residue);
		}
		protein.chains.add(chain);
		protein.size = sequence.length();
		protein.reconstructionLevel = ReconstructionLevel.NONE;
		Project project = new Project();
		project.date = new Date().getTime();
		project._id = UUID.randomUUID().toString();
		project.name = projectName;
		project.proteins.add(protein);
		return project;
	}

	@Override
	public Atom getAtomByName(Residue residue, String name) {
		return residue.atoms.stream().filter(r -> r.name.equals(name)).findFirst().get();
	}

	@Override
	public List<Atom> getAtoms(Protein protein) {
		return getResidues(protein).stream().flatMap(this::asAtomStream).collect(Collectors.toList());
	}

	@Override
	public Atom getC(Residue residue) {
		// this is not really elegant - quite boiler-plate
		return residue.atoms.stream().filter(this::isC).findFirst().get();
	}

	@Override
	public Atom getCA(Residue residue) {
//		try {
			return residue.atoms.stream().filter(this::isCA).findFirst().get();
//		} catch (NoSuchElementException e) {
//			System.out.println("atoms:");
//			residue.atoms.forEach(System.out::println);
//			e.printStackTrace();
//			throw new NoSuchElementException();
//		}
	}

	@Override
	public Atom getCAByResidueIndex(Protein protein, int index) throws NoSuchElementException {
		return this.getCA(this.getResidues(protein).get(index));
	}
	
	@Override
	public Atom getH(Residue residue) {
		return residue.atoms.stream().filter(this::isH).findFirst().get();
	}
	
	@Override
	public Atom getN(Residue residue) {
		return residue.atoms.stream().filter(this::isN).findFirst().get();
	}
	
	@Override
	public Atom getO(Residue residue) {
		return residue.atoms.stream().filter(this::isO).findFirst().get();
	}

	@Override
	public List<Residue> getResidues(Protein protein) {
		return protein.chains.stream().flatMap(this::asResidueStream).collect(Collectors.toList());
	}

	private boolean isC(Atom atom) {
		return atom.name.equals(BACKBONE_C_NAME);
	}

	private boolean isCA(Atom atom) {
		return atom.name.equals(BACKBONE_CA_NAME);
	}
	
	private boolean isH(Atom atom) {
		return atom.name.equals(BACKBONE_H_NAME);
	}

	private boolean isN(Atom atom) {
		return atom.name.equals(BACKBONE_N_NAME);
	}
	
	private boolean isO(Atom atom) {
		return atom.name.equals(BACKBONE_O_NAME);
	}

	private void rearrangeAtoms(Residue residue) {
		Collections.sort(residue.atoms, ATOM_NAME_COMPARATOR);		
	}
	
	@Override
	public boolean removeAtomByName(Residue residue, String name) {
		return residue.atoms.removeIf(a -> a.name.equals(name));
	}

	@Override
	public void removeAtoms(Protein protein) {
		this.getResidues(protein).forEach(this::removeAtoms);
	}
	
	private void removeAtoms(Residue residue) {
		residue.atoms.clear();
	}
	
	@Override
	public void updatePdbRepresentation(Protein protein) {
		new SimplePDBConverter().updatePdbRepresentation(protein);
	}
	
	@Override
	public void updatePdbSerials(Protein protein) {
		int atomCount = 0;
		for(Chain chain : protein.chains) {
			for(Residue residue : chain.residues) {
				rearrangeAtoms(residue);
				for(Atom atom : residue.atoms) {
					atom.pdbSerial = atomCount;
					atomCount++;
				}
			}
		}		
	}
}