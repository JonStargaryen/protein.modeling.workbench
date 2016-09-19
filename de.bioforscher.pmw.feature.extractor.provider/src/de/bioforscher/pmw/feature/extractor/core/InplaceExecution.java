package de.bioforscher.pmw.feature.extractor.core;

import java.io.File;
import java.io.IOException;

import de.bioforscher.pmw.model.Protein;

@Deprecated
public interface InplaceExecution extends ProteinWriter {

    /**
     * provides a way to 'silently' execute certain operations such as calling
     * DSSP, eGOR, Topcons etc<br />
     * to realize this, the {@link Protein} is converted to BioJava and written
     * to a temp file. furthermore, the service is called and its results are
     * written to a 2nd temp file which is subsequently parsed. when finished,
     * both temporary resources are discarded and the information was assigned
     * to the protein structure
     *
     * @param protein
     * @return
     */
    default boolean executeInplace(Protein protein) {
        boolean success = false;
        File tmpPdbFile = null;
        File tmpResultFile = null;
        try {
            tmpPdbFile = writeProtein(protein);
            System.out.println("wrote protein file to " + tmpPdbFile.getAbsolutePath());
            tmpResultFile = createTempFile();
            execute(protein, tmpPdbFile, tmpResultFile);
            parse(protein, tmpResultFile);
            success = true;
        } catch (Exception e) {
            e.printStackTrace();
            success = false;
        } finally {
            tmpPdbFile.delete();
            tmpResultFile.delete();
        }
        return success;
    }

    void execute(Protein protein, File inputFile, File outputFile) throws IOException;

    void parse(Protein protein, File outputFile) throws IOException;
}
