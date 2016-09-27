package de.bioforscher.pmw.model;

import java.util.ArrayList;
import java.util.List;

import org.osgi.dto.DTO;

public class Protein extends DTO {
	public List<FeatureType> availableFeatures;
	public List<Chain> chains;
	public Membrane membrane;
	public String name;
	public ReconstructionLevel reconstructionLevel;
	public int size;
	public String title;
	public List<Motif> motifs;
	public List<Contact> contacts;
	public List<Interaction> interactions;
	//TODO redundant to motifs - keep one
	public List<Fragment> fragments;
	
	public Protein() {
		this.availableFeatures = new ArrayList<>();
		this.chains = new ArrayList<>();
		this.motifs = new ArrayList<>();
		this.contacts = new ArrayList<>();
		this.interactions = new ArrayList<>();
		this.fragments = new ArrayList<>();
	}
	
	@Override
	public String toString() {
		return this.getClass().getSimpleName() + " name='" + this.name + "' size='" + this.size + "'";
	}
}
