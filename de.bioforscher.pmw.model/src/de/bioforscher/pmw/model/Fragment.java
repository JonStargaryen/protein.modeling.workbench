package de.bioforscher.pmw.model;

import java.util.List;

import org.osgi.dto.DTO;

/**
 * A fragment (i.e. a collection of sequentially neighbored residues).
 * @author S
 *
 */
public class Fragment extends DTO {
	public DefinedMotif sequenceMotif;
	public String sequence;
	/**
	 * describes the origin of this fragment
	 */
	public String _id;
	public List<Residue> residues;
	
	//TODO probably we have to consider membrane topology
	
	/**
	 * Use this function to create fragments. It will ensure consistent naming.
	 * @param protein
	 * @param chain
	 * @param residues the residues part of this motif
	 * @param sequenceMotif
	 * @return
	 */
	public static Fragment of(Protein protein, Chain chain, List<Residue> residues, String sequence, DefinedMotif sequenceMotif) {
		Fragment fragment = new Fragment();
		fragment._id = protein.name + "_" + chain.chainId + "_" + residues.get(0).residueNumber + "-" +
				residues.get(residues.size() - 1).residueNumber + "_" + sequence + "_" + sequenceMotif.toString();
		fragment.sequence = sequence;
		fragment.residues = residues;
		fragment.sequenceMotif = sequenceMotif;
		return fragment;
	}
	
	@Override
	public String toString() {
		return this.getClass().getSimpleName() + " sequenceMotif='" + this.sequenceMotif + "' sequence='" + this.sequence + "'";
	}
}
