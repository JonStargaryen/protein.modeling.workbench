package de.bioforscher.pmw.api;

import java.io.IOException;
import java.io.File;
import java.util.List;
import java.util.NoSuchElementException;

import de.bioforscher.pmw.model.Atom;
import de.bioforscher.pmw.model.Project;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.Residue;

/**
 * The service for model manipulation. Basically, any operations which will
 * create or modify the model are provided by this service. Also, lots of
 * convenience methods are specified which allow retrieval of individual atoms
 * or residues by defined rules. More or less, every low-level modification to the model
 * is provided by this class.<br />
 * Low-level computations are part of {@link LinearAlgebra}.<br />
 * High-level access to certain capabilities can be found in {@link ContactPredictor},
 * {@link FeatureExtractor} and {@link ReconstructionService}.<br />
 * <br />
 * Parts of the API are:
 * <ul>
 * <li><b>create</b> new model project instances: either by sequence or PDB file</li>
 * <li><b>update</b> PDB representation of {@link Protein} and {@link Residue} objects - each stores its relevant information in a PDB-formatted <code>ATOM</code> record</li>
 * <li><b>convert</b> of 3-letter and 1-letter codes (the model stores 3-letter codes)</li>
 * <li><b>get</b> convenient access to {@link Residue} and {@link Atom} objects of a {@link Protein} - backbone atoms have dedicated functions</li>
 * <li><b>remove</b> defined atoms or residues from the structure
 * </ul>
 *
 * @author S
 *
 */
public interface ModelConverter {
	/*
	 * default values of features not present in reconstructed structure
	 */
	/**
	 * the default occupancy of a reconstructed {@link Atom}
	 */
	float DEFAULT_OCCUPANCY = 1.0f;
	/**
	 * the default temp factor of a reconstructed {@link Atom}
	 */
	float DEFAULT_TEMP_FACTOR = 0.0f;
	/**
	 * how name a {@link Protein} created solely be its sequence
	 */
	String DEFAULT_PROTEIN_NAME = "Protein";
	/**
	 * how describe a {@link Protein} created solely be its sequence
	 */
	String DEFAULT_PROTEIN_TITLE = "NO INFORMATION AVAILABLE";
	/**
	 * how name a {@link Chain} of a {@link Protein} created solely be its sequence
	 */
	String DEFAULT_CHAIN_ID = "A";
	/**
	 * 1-letter code of non-standard amino acids
	 */
	String UNKNOWN_AMINO_ACID_ONE_LETTER_CODE = "X";
	/**
	 * 3-letter code of non-standard amino acids
	 */
	String UNKNOWN_AMINO_ACID_THREE_LETTER_CODE = "XXX";

	/**
	 * handle to the backbone CA name constant
	 */
	String BACKBONE_CA_NAME = "CA";
	/**
	 * handle to the backbone O name constant
	 */
	String BACKBONE_O_NAME = "O";
	/**
	 * handle to the backbone N name constant
	 */
	String BACKBONE_N_NAME = "N";
	/**
	 * handle to the backbone C name constant
	 */
	String BACKBONE_C_NAME = "C";
	/**
	 * handle to the backbone H name constant
	 */
	String BACKBONE_H_NAME = "H";
	/**
	 * handle to the CB name constant
	 */
	String SIDECHAIN_CB_NAME = "CB";

	/**
	 * the PDB URL which can be used to fetch structures by ID (format this using the id and you are good to go)
	 */
	String PDB_FETCH_URL = "https://files.rcsb.org/download/%s.pdb";
	
	/**
	 * convenience function to convert 3letter-code to 1letter-code
	 *
	 * @param aminoAcidName
	 * @return
	 */
	String convertToOneLetterCode(String threeLetterCode);
	
	/**
	 * convenience function to convert 1letter-code to 3letter-code
	 *
	 * @param aminoAcidName
	 * @return
	 */
	String convertToThreeLetterCode(String oneLetterCode);
	
	/**
	 * add an atom to the specified {@link Residue}
	 * @param residue the parent to which this atom will be added
	 * @param name the atom name (C, CA, CB, O, N)
	 * @param coordinates the spatial coordinates of the atom
	 */
	void createAtom(Residue residue, String name, double[] coordinates);

	/**
	 * creates an instance of an atom with the given name and coordinates
	 * @param name the atom name (such as C, CA, CB, N etc)
	 * @param xyz the spatial coordinates as <code>double[]</code>
	 * @return
	 */
	Atom createAtom(String name, double[] xyz);

	/**
	 * creates a project instance by a given PDB structure
	 * 
	 * @param pdbFileContent
	 *            binary representation of a PDB file
	 * @return
	 * @throws IOException
	 *             when parsing the file fails
	 */
	Project createModelingProject(byte[] pdbFileContent) throws IOException;
	
	/**
	 * creates a protein instance by fetching the given pdbId from the PDB and parsing it
	 * @param pdbId 4 digit pdb code
	 * @return the file content
	 * @throws Exception
	 */
	Protein createProteinByPDBId(String pdbId) throws Exception;
	
	Protein createProtein(File file) throws Exception;
	
//	/**
//	 * uses the specified sequence to initialize a {@link Protein}-scaffold
//	 * which can be subsequently refined<br />
//	 *
//	 * @param proteinName
//	 *            the name which will be set to the {@link Protein}s name
//	 *            attribute (e.g. a PDB code or job name)
//	 * @param sequence
//	 *            the sequence to be processed - currently, a single chain in
//	 *            plain String-format can be handled
//	 * @return a {@link Project} containing a {@link Protein} will represent the sequence by a single
//	 *         {@link Chain} consisting of {@link Residue}s with correct
//	 *         residue-numbering and amino acid-attributes
//	 */
//	Project createModelingProject(String proteinName, String proteinTitle, String sequence);
	
	/**
	 * Starts a new modeling project consisting of one protein containing all input information and data shared by all structure models about to be created.
	 * @param projectName the project's name
	 * @param sequence the sequence to be processed
	 * @return a {@link Project} containing a {@link Protein} will represent the sequence by a single 
	 * {@link Chain} consisting of {@link Residue}s with correct residue-numbering and amino acid-attributes
	 */
	Project createModelingProject(String projectName, String sequence);
	
	/**
	 * convenience method to get a residue's atom by its unique name (defined as constant in {@link ModelConverter})
	 * <b>important:</b> in its current implementation use the dedicated functions for specific atom name if possible
	 * @param residue the container to process
	 * @param name the unique atom name (C, N, O, CA)
	 * @return the atom if present - else a {@link NoSuchElementException} is thrown
	 */
	Atom getAtomByName(Residue residue, String name);
	
	/**
	 * returns all atoms of this protein
	 * 
	 * @param protein
	 *            the container
	 * @return a collection of all atoms present in the protein (however, they have to be
	 *         part of {@link Residue}s)
	 */
	List<Atom> getAtoms(Protein protein);
	
	Atom getC(Residue residue);
	
	Atom getCA(Residue residue);
	
	/**
	 * 
	 * @param protein the context of the query
	 * @param residueId the id of the residue whose CA atom should be returned
	 * @return
	 */
	Atom getCAByResidueIndex(Protein protein, int residueId);
	
	Atom getH(Residue residue);
	
	Atom getN(Residue residue);

	Atom getO(Residue residue);

	/**
	 * returns all residues of this protein
	 * 
	 * @param protein
	 *            the container
	 * @return a collection of all residues present in the protein
	 */
	List<Residue> getResidues(Protein protein);

	/**
	 * delete atoms from a residue by a given name
	 * @param residue the parent from which atoms will potentially be removed
	 * @param name what atoms will be hunted
	 * @return true if something was removed
	 */
	boolean removeAtomByName(Residue residue, String name);
	
	/**
	 * removes all atoms entries (and thus their coordinates) from the structure
	 */
	void removeAtoms(Protein protein);
	
//	/**
//	 * creates the pdbRepresentation of this {@link Protein} and updates the
//	 * corresponding field
//	 * 
//	 * @param protein
//	 *            the object whose pdbRepresentation attribute will be updated
//	 * @return the updated pdbRepresentation - it is still set to the field of
//	 *         the {@link Protein} object
//	 */
//	String updatePdbRepresentation(Protein protein);
	/**
	 * creates the pdbRepresentation of this {@link Protein} and updates the
	 * corresponding fields
	 * 
	 * @param protein
	 *            the object whose atom's pdbRepresentation attributes will be updated
	 */
	void updatePdbRepresentation(Protein protein);
	
	/**
	 * Ensures correct numbering of all present atoms within the structure and correct ordering in the final PDB file: N -> CA -> C -> O -> CB+.
	 * @param protein the protein to be processed
	 */
	void updatePdbSerials(Protein protein);
}
