package de.bioforscher.pmw.fragment.provider;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.stream.Collectors;

import org.osgi.service.component.annotations.Activate;
import org.osgi.service.component.annotations.Component;
import org.osgi.service.component.annotations.Reference;
import org.osgi.service.log.LogService;

import de.bioforscher.pmw.api.FeatureExtractor;
import de.bioforscher.pmw.api.FragmentLibrary;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.api.ModelPersistence;
import de.bioforscher.pmw.model.FeatureType;
import de.bioforscher.pmw.model.Fragment;
import de.bioforscher.pmw.model.Protein;

@Component(name = "de.bioforscher.pmw.fragmentlibrary")
public class FragmentLibraryImpl implements FragmentLibrary {
	/**
	 * custom list of non-redundant alpha-helical TM proteins - produced by PDB's drilldown search
	 */
	public static final String PDB_IDS = "de/bioforscher/pmw/fragment/provider/pdb_ids.dat";
	/**
	 * standard PDBTM id list - non-redundant, alpha-helical
	 */
	public static final String PDBTM_IDS = "de/bioforscher/pmw/fragment/provider/pdbtm_ids.dat";

	@Reference
	private ModelPersistence modelPersistence;
	@Reference
	private LogService logger;
	@Reference
	private ModelConverter modelConverter;
	@Reference
	private FeatureExtractor featureExtractor;
	
	@Activate
	public void activate() throws Exception {
		try {
			System.out.println("loading fragment library - " + this.modelPersistence.retrieveFragments().size() + " known fragments");
			if(this.modelPersistence.retrieveFragments().size() == 0) {
				System.out.println("creating new fragment library");
//				createFragmentLibrary();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	@SuppressWarnings("unused")
	private void createFragmentLibrary() {
		InputStream is = Thread.currentThread().getContextClassLoader().getResourceAsStream(PDBTM_IDS);
		new BufferedReader(new InputStreamReader(is)).lines().parallel().forEach(this::handleIdLine);
	}

	private void handleIdLine(String line) {
		if(line.startsWith("#")) {
			System.out.println(line);
		} else {
			try {
				String pdbId = line.split("_")[0];
				String chainId = line.split("_")[1];
				Protein protein = this.modelConverter.createProteinByPDBId(pdbId);
				System.out.println(protein.name + " - " + protein.size + " residues");
				fragmentize(protein, chainId);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	private void fragmentize(Protein protein, String chainId) {
		this.featureExtractor.computeFeatures(protein, FeatureType.MOTIF_ANNOTATION);

		final String basePath = "D:/fragments/";
		
		for(Fragment fragment : protein.fragments) {
			if(!fragment._id.split("_")[1].equals(chainId)) {
				continue;
			}
			
			// create/ensure dir
			File motifDir = new File(basePath + fragment.sequenceMotif);
			if(!motifDir.exists()) {
				System.out.println("creating dir: " + motifDir.getAbsolutePath());
				motifDir.mkdir();
			}
			
			File fragmentFile = new File(motifDir + "/" + fragment._id + ".pdb");
			String fragmentFileContent = fragment.residues.stream().flatMap(r -> r.atoms.stream()).map(a -> a.pdbRepresentation).collect(Collectors.joining(System.lineSeparator()));
			try {
				System.out.println("writing file: " + fragmentFile.getAbsolutePath());
				Files.write(fragmentFile.toPath(), fragmentFileContent.getBytes());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}		
	}
}
