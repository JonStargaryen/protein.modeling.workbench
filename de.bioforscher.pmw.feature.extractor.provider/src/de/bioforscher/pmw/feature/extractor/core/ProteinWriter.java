package de.bioforscher.pmw.feature.extractor.core;

import java.io.File;
import java.io.IOException;
import de.bioforscher.pmw.model.Protein;

@Deprecated
public interface ProteinWriter {

    String TEMP_FILE_PREFIX = "java-pmw-";
    String TEMP_FILE_SUFFIX = ".tmp";

    default File writeProtein(Protein protein) throws IOException {
        File tmpFile = createTempFile();
        //TODO change to some existing value if we would ever use this again
//        Files.write(tmpFile.toPath(), protein.pdbRepresentation.getBytes());
        return tmpFile;
    }

    default File createTempFile() throws IOException {
        return File.createTempFile(TEMP_FILE_PREFIX, TEMP_FILE_SUFFIX);
    }
}
