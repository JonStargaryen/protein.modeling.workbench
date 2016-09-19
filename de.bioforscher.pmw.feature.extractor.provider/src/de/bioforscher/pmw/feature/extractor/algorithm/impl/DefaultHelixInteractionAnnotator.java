package de.bioforscher.pmw.feature.extractor.algorithm.impl;

import org.slf4j.Logger;

import de.bioforscher.pmw.api.FeatureExtractor;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.feature.extractor.core.AbstractFeatureProvider;
import de.bioforscher.pmw.feature.extractor.core.Annotator;
import de.bioforscher.pmw.model.FeatureType;
import de.bioforscher.pmw.model.Protein;

@Deprecated
public class DefaultHelixInteractionAnnotator extends AbstractFeatureProvider implements Annotator {

	public DefaultHelixInteractionAnnotator(FeatureExtractor featureExtractor, Logger logger, LinearAlgebra linearAlgebra, ModelConverter modelConverter) {
		super(featureExtractor, logger, linearAlgebra, modelConverter, new FeatureType[] { null, null });
//		super(new FeatureType[] { FeatureType.HELIX_HELIX_INTERACTION_ANNOTATION }, FeatureType.HELIX_ANNOTATION, FeatureType.RESIDUE_CONTACTS);
	}

	@Override
	protected void computeFeatureInternal(Protein protein) {
//        determineHelixHelixInteractions(protein);
	}
	
	/**
     * scans the previously annotated helices for interactions<br />
     * <br />
     * according to <i>TMPad: an integrated structural database for
     * helix-packing folds in transmembrane proteins - Nucleic Acids Research,
     * 2011, Vol. 39, Database issue D347–D355 doi:10.1093/nar/gkq1255</i> a
     * helix-helix interaction can be observed, when either: <br />
     * <ul>
     * <li>distances between Cb of each helix < 6A</li>
     * <li>the distance between any two heavy atoms (one from each helix) is
     * less than the sum of their VDW radii plus 0.6 A</li>
     * </ul>
     *
     * @param protein
     * @return true when at least 1 interacting helix pair was found
     */
//    private boolean determineHelixHelixInteractions(Protein protein) {
//        // when less than 2 helices could be observed, no helix-helix
//        // interactions can occur
//        int numberOfHelices = protein.helices.size();
//        if (numberOfHelices < 2) {
//            return false;
//        }
//
//        for (int i = 0; i < numberOfHelices - 1; i++) {
//            TransmembraneHelix helix1 = protein.getHelices().get(i);
//            for (int j = i + 1; j < numberOfHelices; j++) {
//                TransmembraneHelix helix2 = protein.getHelices().get(j);
//                checkWhetherHelicesInteract(protein, helix1, helix2);
//            }
//        }
//        
//        return protein.getHelixHelixInteractions().size() > 0;
//    }
//
//    private void checkWhetherHelicesInteract(Protein protein, TransmembraneHelix helix1, TransmembraneHelix helix2) {
//        List<ResidueResidueInteraction> residueInteractions = new ArrayList<>();
//        for (Residue residue1 : helix1.getResidues()) {
//            Optional<Atom> cb1 = retrieveAtomByName(residue1, BETA_CARBON_NAME);
//            outer: for (Residue residue2 : helix2.getResidues()) {
//                Optional<Atom> cb2 = retrieveAtomByName(residue2, BETA_CARBON_NAME);
//
//                // interaction according to criterion 1 - cb-cb <6 A
//                if (computeDistance(cb1, cb2) < HELIX_HELIX_INTERACTION_CB_CUTOFF) {
//                    residueInteractions.add(createResidueResidueInteraction(residue1, residue2,
//                            ResidueResidueInteractionType.CBETA_DISTANCE));
//                    continue outer;
//                }
//
//                // interaction according to criterion 2 - vdw-vdw <0.6 A
//                for (int i = 0; i < residue1.getAtoms().size() - 1; i++) {
//                    for (int j = i + 1; j < residue2.getAtoms().size(); j++) {
//                        Atom a1 = residue1.getAtoms().get(i);
//                        Atom a2 = residue2.getAtoms().get(j);
//                        if (computeDistance(Optional.of(a1), Optional.of(a2)) - lookUpVanDerWaalsRadius(a1)
//                                - lookUpVanDerWaalsRadius(a2) < HELIX_HELIX_INTERACTION_VDW_CUTOFF) {
//                            residueInteractions.add(createResidueResidueInteraction(residue1, residue2,
//                                    ResidueResidueInteractionType.VDW_DISTANCE));
//                            continue outer;
//                        }
//                    }
//                }
//            }
//        }
//
//        HelixHelixInteraction helixHelixInteractions = ProteinFactory.eINSTANCE.createHelixHelixInteraction();
//        helixHelixInteractions.getResidueResidueInteractions().addAll(residueInteractions);
//        helixHelixInteractions.getInteractingHelices()
//                .addAll(Arrays.asList(new TransmembraneHelix[] { helix1, helix2 }));
//        protein.getHelixHelixInteractions().add(helixHelixInteractions);
//    }
}
