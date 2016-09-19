package de.bioforscher.pmw.model;

import org.osgi.dto.DTO;

public class Atom extends DTO {
	public String element;
	public String name;
	public float occupancy;
	public int pdbSerial;
	public float tempFactor;
	public double[] xyz;
	public String pdbRepresentation;
}
