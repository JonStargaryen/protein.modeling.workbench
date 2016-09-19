package de.bioforscher.pmw.model;

import java.util.ArrayList;
import java.util.List;

import org.osgi.dto.DTO;

public class Membrane extends DTO {
	public double[] centerOfMass;
	public List<double[]> membraneMolecules;
	public double[] normalVector;
	public double[] planePoint1;
	public double[] planePoint2;
	public double[] spherePoint;
	
	public Membrane() {
		this.membraneMolecules = new ArrayList<>();
	}
}
