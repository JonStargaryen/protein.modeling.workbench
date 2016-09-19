package de.bioforscher.pmw.model.converter.parser;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Locale;
import de.bioforscher.pmw.model.Atom;
import de.bioforscher.pmw.model.Chain;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.Residue;

/**
 * A minimal PDB parser which discards most "useless" information and rather only parses and returns information strictly captured by the model. Allows also to write ATOM records to individual atoms. Most code taken from BioJava's parsing/writing capabilities - however, significantly reduced as why do not depend on many features.
 * @author S
 *
 */
public class SimplePDBConverter implements PDBConverter {
	private Protein protein;
	private Chain currentChain;
	private Residue currentResidue;
	// the model/API uses this custom field of residues to identify them in the structure (it is unique across all chains, e.g. "C-123" could be utilized alternatively)
	private int residueId;
	private static final String HEADER_PREFIX = "HEADER";
	private static final String TITLE_PREFIX = "TITLE";
	private static final String ATOM_PREFIX = "ATOM";
	
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
	
	@Override
	public Protein parsePDBFile(InputStream inputStream) throws IOException {
		this.protein = new Protein();
		this.residueId = 0;
		// 'initialize' title field as it tends to be split over multiple lines - thus, we have to append previous results when we find further entries
		this.protein.title = "";
		
		// parse file
		BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
		String line;
		while((line = br.readLine()) != null) {
			parseLine(line);
		}
		
		return this.protein;
	}
	
	private void parseLine(String line) {
		// indices taken from: ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
		// their column definition has certain offset to the definition of String#substring(int, int)

		// found header record
		// 63 - 66       IDcode        idCode            This identifier is unique within the PDB.
		if(line.startsWith(HEADER_PREFIX)) {
			protein.name = line.substring(62, 66);
		}
		
		// handling title records
		// 11 - 80       String        title         Title of the experiment.
		if(line.startsWith(TITLE_PREFIX)) {
			// trim to omit tailing white-spaces
			// extra whitespace to ensure that words are separated
			// maybe some StringJoiner is the way to go
			protein.title += " " + line.substring(10, 80).trim();
			protein.title = protein.title.trim();
		}
		
		// parsing atom record - information we need is marked with an '*' - indirectly needed information (chain/residue) marked with an '#'
		// some information will inform us about changing chain/residue
		/*	COLUMNS        DATA TYPE     FIELD        DEFINITION
			-------------------------------------------------------------------------------------
			1 - 6          Record name   "ATOM  "
		*	7 - 11   	   Integer       serial       Atom serial number.
		*	13 - 16        Atom          name         Atom name.
			17             Character     altLoc       Alternate location indicator.
		#	18 - 20        Residue name  resName      Residue name.
		#	22             Character     chainID      Chain identifier.
		#	23 - 26        Integer       resSeq       Residue sequence number.
			27             AChar         iCode        Code for insertion of residues.
		*	31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		*	39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
		*	47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		*	55 - 60        Real(6.2)    occupancy     Occupancy.
		*	61 - 66        Real(6.2)    tempFactor    Temperature factor.
		*	77 - 78        LString(2)   element       Element symbol, right justified.
			79 - 80        LString(2)   charge        Charge on the atom */
		if(line.startsWith(ATOM_PREFIX)) {
			String chainId = line.substring(21, 22);
			int resNum = Integer.parseInt(line.substring(22, 26).trim());
			if(this.currentChain == null || !this.currentChain.chainId.equals(chainId)) {
				// chain changed - create new chain object and set reference
				this.currentChain = new Chain();
				
				// parse chain
				this.currentChain.chainId = chainId;
				
				this.protein.chains.add(this.currentChain);
			}
			
			if(this.currentResidue == null || this.currentResidue.residueNumber != resNum) {
				// residue changed - create new residue object and set reference
				this.currentResidue = new Residue();
				
				// parse residue
				this.currentResidue.aminoAcid = line.substring(17, 20).trim();
				this.currentResidue.insertionCode = null;
				this.currentResidue.residueId = this.residueId;
				this.currentResidue.residueNumber = resNum;
				
				this.currentChain.residues.add(this.currentResidue);
			}
			
			// we append the current residue container with additional atoms
			Atom atom = new Atom();
			
			// parse atom
			atom.element = line.substring(76, 78).trim();
			atom.name = line.substring(12, 16).trim();
			atom.occupancy = Float.valueOf(line.substring(54, 60).trim());
			atom.pdbRepresentation = line;
			atom.pdbSerial = Integer.valueOf(line.substring(6, 11).trim());
			atom.tempFactor = Float.valueOf(line.substring(60, 66).trim());
			atom.xyz = new double[] { Double.valueOf(line.substring(30, 38).trim()),
					Double.valueOf(line.substring(38, 46).trim()),
					Double.valueOf(line.substring(46, 54).trim())
				};
			
			this.currentResidue.atoms.add(atom);
		}
	}
	
	@Override
	public void updatePdbRepresentation(Protein protein) {
		for (Chain chain : protein.chains) {
			for (Residue residue : chain.residues) {
				for (Atom atom : residue.atoms) {
					atom.pdbRepresentation = composePdbRepresentation(atom, residue, chain);
				}
			}
		}
	}
	
	private String composePdbRepresentation(Atom atom, Residue residue, Chain chain) {
		String record = "ATOM  ";
		// format output ...
		String resName = residue.aminoAcid;
		String pdbcode = String.valueOf(residue.residueNumber);

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
		s.append(chain.chainId);
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
}
