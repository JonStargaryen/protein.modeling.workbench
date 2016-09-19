package parser.test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.nio.charset.StandardCharsets;
import java.util.Objects;

import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.converter.parser.PDBConverter;
import de.bioforscher.pmw.model.converter.parser.SimplePDBConverter;
import junit.framework.TestCase;

public class SimplePDBParserFunctionalTest extends TestCase {
	// local file to parse
	private static final String PDB_PATH = "parser/test/4cha.pdb";
	// reference to the list of all PDB structures
	private static final String PDB_ID_LIST_URL = "ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt";
	
	public void testSimplePDBParser() {
		try {
			PDBConverter parser = new SimplePDBConverter();
			Protein protein = parser.parsePDBFile(getResourceAsFilepath(PDB_PATH));
			assertValidPDBFile(protein);
		} catch (IOException e) {
			fail("failed with IOException: " + e.getLocalizedMessage());
		}
	}
	
	public void testOnWholePDB() throws IOException {
		URLConnection conn = new URL(PDB_ID_LIST_URL).openConnection();
		try (BufferedReader reader = new BufferedReader(new InputStreamReader(conn.getInputStream(), StandardCharsets.UTF_8))) {
		    reader.lines().forEach(line -> {
		    	String id = line.split("\t")[0];
		    	System.out.println(id);
		    	try {
					assertValidPDBFile(new SimplePDBConverter().parsePDBFile(new URL("https://files.rcsb.org/view/" + id + ".pdb").openStream()));
				} catch (IOException e) {
					fail("failed with IOException: " + e.getLocalizedMessage());
				}
		    });
		}
	}
	
	private void assertValidPDBFile(Protein protein) {
		// ensure name is 4 letter pdb code
		assert(protein.name.length() == 4);
		System.out.println("name: " + protein.name);
		
		// ensure title is parsed
		assert(protein.title.length() > 0);
		System.out.println("title: " + protein.title);
		
		// ensure atoms were parsed - at this point there is no possibility to decide whether the all ATOM records were parsed
		System.out.println(protein.chains.stream().mapToInt(c -> c.residues.size()).sum() + " residues in " + protein.chains.size() + " chains");
	}
	
	private String getResourceAsFilepath(String filename) {
		ClassLoader ccl = Thread.currentThread().getContextClassLoader();
		Objects.requireNonNull(ccl);
		URL resource = ccl.getResource(filename);
		Objects.requireNonNull(resource);
		return resource.getPath();
	}
}
