package de.bioforscher.pmw.alignment.provider;

import java.io.File;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.DoubleSummaryStatistics;
import java.util.List;
import java.util.Locale;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import de.bioforscher.pmw.api.AlignmentService;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.model.Alignment;
import de.bioforscher.pmw.model.Atom;
import de.bioforscher.pmw.model.Protein;

public class StructureComposition {
	private ModelConverter modelConverter;
	private AlignmentService alignmentService;
	private Protein reference;
	private List<Double> rmsds;
	
	public StructureComposition(ModelConverter modelConverter, AlignmentService alignmentService) {
		this.modelConverter = modelConverter;
		this.alignmentService = alignmentService;
//		System.out.println("motif,length,averageRmsd");
		Stream.of(new File("D:/fragments/").listFiles()).forEach(this::handleDirectory);
	}
	
	private void handleDirectory(File file) {
//		System.out.println(file.getName());
		this.reference = null;
		this.rmsds = new ArrayList<>();
		Stream.of(file.listFiles()).forEach(this::handleFile);
		DoubleSummaryStatistics stats = this.rmsds.stream().mapToDouble(Double::new).summaryStatistics();
		System.out.println(file.getName() + "," + file.getName().charAt(2) + "," + stats.getAverage());
	}
	
	private void handleFile(File file) {
//		System.out.println(file.getName());
		try {
			Protein protein = this.modelConverter.createProtein(file);
			List<Atom> atoms = this.modelConverter.getResidues(protein).stream().map(r -> this.modelConverter.getCA(r)).collect(Collectors.toList());
			if(this.reference == null) {
//				System.out.println("using " + file.getName() + " as reference for alignment");
				this.reference = protein;
			} else {
				Alignment alignment = this.alignmentService.alignFragments(this.modelConverter.getResidues(this.reference), this.modelConverter.getResidues(protein));
				this.rmsds.add(alignment.rmsd);
				this.alignmentService.transform(atoms, alignment.translationVector, alignment.rotationMatrix);
//				System.out.println(reference.name + " <> " + protein.name + " : " + alignment.rmsd);
			}
			new File("D:/fragments-superimposed/" + file.getParentFile().getName() + "/").mkdir();
			updatePdbRepresentation(atoms);
			Files.write(new File("D:/fragments-superimposed/" + file.getParentFile().getName() + "/" + file.getName()).toPath(), atoms.stream().map(a -> a.pdbRepresentation).collect(Collectors.joining("\n")).getBytes());
		} catch (IllegalArgumentException e) {
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void updatePdbRepresentation(List<Atom> atoms) {
		atoms.forEach(a -> a.pdbRepresentation = composePdbRepresentation(a));
	}
	
	// Locale should be english, e.g. in DE separator is "," -> PDB files have
	// "." !
	public static DecimalFormat d3 = (DecimalFormat) NumberFormat.getInstance(Locale.US);
	static {
		d3.setMaximumIntegerDigits(4);
		d3.setMinimumFractionDigits(3);
		d3.setMaximumFractionDigits(3);
	}
	public static DecimalFormat d2 = (DecimalFormat) NumberFormat.getInstance(Locale.US);
	static {
		d2.setMaximumIntegerDigits(3);
		d2.setMinimumFractionDigits(2);
		d2.setMaximumFractionDigits(2);
	}
	
	private String formatAtomName(Atom a) {
		String fullName = null;
		String name = a.name;
		String element = a.element;

		// RULES FOR ATOM NAME PADDING: 4 columns in total: 13, 14, 15, 16

		// if length 4: nothing to do
		if (name.length() == 4) {
			fullName = name;
		}

		// if length 3: they stay at 14
		else if (name.length() == 3) {
			fullName = " " + name;
		}

		// for length 2 it depends:
		// carbon, oxygens, nitrogens, phosphorous stay at column 14
		// elements with 2 letters (e.g. NA, FE) will go to column 13
		else if (name.length() == 2) {
			if (element.equals("C") || element.equals("N") || element.equals("O") || element.equals("P")
					|| element.equals("S")) {
				fullName = " " + name + " ";
			} else {
				fullName = name + "  ";
			}
		}

		// for length 1 (e.g. K but also C, O) they stay in column 14
		else if (name.length() == 1) {
			fullName = " " + name + "  ";
		}
		// if (fullName.length()!=4)
		// logger.warn("Atom name "+fullName+"to be written in PDB format does
		// not have length 4. Formatting will be incorrect");

		return fullName;
	}
	
	private String composePdbRepresentation(Atom atom) {
		String record = "ATOM  ";
		// format output ...
		String resName = "???";
		String pdbcode = String.valueOf("0");

		int seri = atom.pdbSerial;
		String serial = String.format("%5d", seri);
		String fullName = formatAtomName(atom);

		Character altLoc = ' ';
		String resseq = String.format("%4s", pdbcode) + " ";

		String x = String.format("%8s", d3.format(atom.xyz[0]));
		String y = String.format("%8s", d3.format(atom.xyz[1]));
		String z = String.format("%8s", d3.format(atom.xyz[2]));
		String occupancy = String.format("%6s", d2.format(atom.occupancy));
		String tempfactor = String.format("%6s", d2.format(atom.tempFactor));
		String leftResName = String.format("%3s", resName);

		StringBuffer s = new StringBuffer();
		s.append(record);
		s.append(serial);
		s.append(" ");
		s.append(fullName);
		s.append(altLoc);
		s.append(leftResName);
		s.append(" ");
		s.append("?");
		s.append(resseq);
		s.append("   ");
		s.append(x);
		s.append(y);
		s.append(z);
		s.append(occupancy);
		s.append(tempfactor);

		String e = atom.element;

		return String.format("%-76s%2s", s.toString(), e) + " ";
	}
}
