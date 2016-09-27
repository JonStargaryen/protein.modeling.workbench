package de.bioforscher.pmw.alignment.provider;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.stream.Stream;

/**
 * Composes all motif sequences for given fragments into one file.
 * @author S
 *
 */
public class RunSequenceComposition {
	private StringBuilder output;
	
	public static void main(String[] args) {
		new RunSequenceComposition().run();
	}
	
	private void run() {
		Stream.of(new File("D:/fragments/").listFiles()).forEach(this::handleDirectory);
	}
	
	private void handleDirectory(File file) {
		System.out.println(file.getName());
		this.output = new StringBuilder();
		Stream.of(file.listFiles()).forEach(this::handleFile);
		try {
			Files.write(new File("D:/fragments-stat/fragments-seq-reduced/" + file.getName() + ".seq").toPath(), this.output.toString().getBytes());
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void handleFile(File file) {
		System.out.println(file.getName());
		
		// normal / no mapping
//		this.output.append(file.getName().split("_")[3] + "\n");
		
		// map to reduced alphabet - c = charged, p = polar, h = hydrophobic
		this.output.append(mapToReducedAlphabet(file.getName().split("_")[3]) + "\n");
	}
	
	String mapToReducedAlphabet(String rawSequence) {
		return rawSequence.charAt(0) + rawSequence.substring(1, rawSequence.length() - 1).replace('R', 'c').replace('K', 'c').replace('D', 'c').replace('E', 'c')
				.replace('Q', 'p').replace('N', 'p').replace('H', 'p').replace('S', 'p').replace('T', 'p').replace('Y', 'p').replace('C', 'p').replace('M', 'p').replace('W', 'p')
				.replace('A', 'h').replace('I', 'h').replace('L', 'h').replace('F', 'h').replace('V', 'h').replace('P', 'h').replace('G', 'h') + rawSequence.charAt(rawSequence.length() - 1);
//		http://www.proteinstructures.com/Structure/Structure/amino-acids.html
//			• Arginine - Arg - R
//			• Lysine - Lys - K
//			• Aspartic acid - Asp - D
//			• Glutamic acid - Glu - E
//
//			Polar (may participate in hydrogen bonds):
//			• Glutamine - Gln - Q
//			• Asparagine - Asn - N
//			• Histidine - His - H
//			• Serine - Ser - S
//			• Threonine - Thr - T
//			• Tyrosine - Tyr - Y
//			• Cysteine - Cys - C
//			• Methionine - Met - M
//			• Tryptophan - Trp - W
//
//			Hydrophobic (normally buried inside the protein core):
//			• Alanine - Ala - A
//			• Isoleucine - Ile - I
//			• Leucine - Leu - L
//			• Phenylalanine - Phe - F
//			• Valine - Val - V
//			• Proline - Pro - P
//			• Glycine - Gly - G 
	}
}
