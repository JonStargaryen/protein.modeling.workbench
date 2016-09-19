package de.bioforscher.pmw.model;

import java.util.ArrayList;
import java.util.List;

import org.osgi.dto.DTO;

public class Chain extends DTO {
	public String chainId;
	public List<Residue> residues;
	
	public Chain() {
		this.residues = new ArrayList<>();
	}
}
