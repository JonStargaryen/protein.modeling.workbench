package de.bioforscher.pmw.model.converter.parser;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;

import de.bioforscher.pmw.model.Protein;

public interface PDBConverter {
	/**
	 * Parses a PDB file and stores all handled information in a {@link Protein} container.
	 * @param inputStream
	 * @return
	 * @throws IOException
	 */
	Protein parsePDBFile(InputStream inputStream) throws IOException;
	
	/**
	 * @see {@link #parsePDBFile(byte[])}
	 */
	default Protein parsePDBFile(File file) throws IOException {
		return parsePDBFile(Files.newInputStream(file.toPath()));
	}
	
	/**
	 * @see {@link #parsePDBFile(byte[])}
	 */
	default Protein parsePDBFile(String filepath) throws IOException {
		return parsePDBFile(new File(filepath));
	}
	
	void updatePdbRepresentation(Protein protein);
}
