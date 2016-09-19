package de.bioforscher.pmw.feature.extractor.algorithm.impl;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.NoSuchElementException;
import java.util.function.Predicate;

import org.slf4j.Logger;

import de.bioforscher.pmw.api.FeatureExtractor;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.feature.extractor.core.AbstractFeatureProvider;
import de.bioforscher.pmw.feature.extractor.core.Annotator;
import de.bioforscher.pmw.feature.extractor.core.InplaceExecution;
import de.bioforscher.pmw.model.Chain;
import de.bioforscher.pmw.model.FeatureType;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.Residue;

/**
 * a wrapper for the DSSP executable - adapting this code seems quite ambitious - for now everything is just delegated to original DSSP - thus, the runtime must be available on the target platform
 * @author S
 *
 */
@Deprecated
public class DSSP extends AbstractFeatureProvider implements InplaceExecution, Annotator {
    
	public DSSP(FeatureExtractor featureExtractor, Logger logger, LinearAlgebra linearAlgebra, ModelConverter modelConverter) {
		super(featureExtractor, logger, linearAlgebra, modelConverter, new FeatureType[] { FeatureType.ACCESSIBLE_SURFACE_AREA, FeatureType.SECONDARY_STRUCTURE });
	}

	@Override
	protected void computeFeatureInternal(Protein protein) {	
		executeInplace(protein);
	}

    private Predicate<? super Residue> findResidue(String line) {
        return r -> r.residueNumber == Integer.parseInt(line.substring(5, 10).trim());
    }

    private Predicate<? super Chain> findChain(String line) {
        return c -> c.chainId.equals(line.substring(11, 12));
    }

	@Override
	public void execute(Protein protein, File inputFile, File outputFile) throws IOException {
		Process process = null;
		//TODO: this path should be part of some config file
		final String[] windowsCommand = new String[] { "d:/eclipse/eclipse/plugins/dssp-2.0.4-win32.exe", inputFile.getAbsolutePath(),
                outputFile.getAbsolutePath() };
		final String[] unixCommand = new String[] { "/usr/local/bin/dssp", "-i", inputFile.getAbsolutePath(), "-o",
                outputFile.getAbsolutePath() };
        try {
            process = Runtime.getRuntime().exec(System.getProperty("os.name").toLowerCase().indexOf("win") >= 0 ?
            		windowsCommand : unixCommand);
            process.waitFor();
        } catch (IOException | InterruptedException e) {
            throw new IOException("executing dssp failed " + e.getMessage() + "\rprobably no DSSP executable is installed");
        } finally {
            try {
                process.destroy();
            } catch (NullPointerException e) {
                throw new RuntimeException("no dssp process could be created");
            }
        }
	}

	@Override
	public void parse(Protein protein, File outputFile) throws IOException {
		for (String line : Files.readAllLines(outputFile.toPath())) {
            if (line.length() != 136 || line.charAt(13) == '!') {
                continue;
            }
            try {
                Chain chain = protein.chains.stream().filter(findChain(line)).findFirst().get();
                Residue residue = chain.residues.stream().filter(findResidue(line)).findFirst().get();
                residue.features.put(FeatureType.SECONDARY_STRUCTURE.name(), wrapInArray(Character.getNumericValue(line.substring(16, 17).charAt(0))));
                residue.features.put(FeatureType.ACCESSIBLE_SURFACE_AREA.name(), wrapInArray(Double.parseDouble(line.substring(34, 38).trim())));
            } catch (NoSuchElementException e) {
//                e.printStackTrace();
                System.err.println("parsing dssp results failed for: " + line);
            }
        }
	}
}
