package de.bioforscher.pmw.model;

import java.util.Arrays;

import org.osgi.dto.DTO;

public class Atom extends DTO {
	public String element;
	public String name;
	public float occupancy;
	public int pdbSerial;
	public float tempFactor;
	public double[] xyz;
	public String pdbRepresentation;
	
	@Override
	public String toString() {
		return this.getClass().getSimpleName() + " name='" + this.name + "' coords='" + Arrays.toString(this.xyz) + "' element='" + this.element + "'";
	}
}
