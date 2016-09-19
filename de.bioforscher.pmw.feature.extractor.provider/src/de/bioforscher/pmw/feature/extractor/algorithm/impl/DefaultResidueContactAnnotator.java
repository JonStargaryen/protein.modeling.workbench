package de.bioforscher.pmw.feature.extractor.algorithm.impl;

import java.util.ArrayList;
import java.util.List;

import org.osgi.service.log.LogService;
import org.slf4j.Logger;

import de.bioforscher.pmw.api.FeatureExtractor;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.feature.extractor.core.AbstractFeatureProvider;
import de.bioforscher.pmw.feature.extractor.core.Annotator;
import de.bioforscher.pmw.model.Atom;
import de.bioforscher.pmw.model.Chain;
import de.bioforscher.pmw.model.FeatureType;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.Residue;

@SuppressWarnings("unused")
@Deprecated
public class DefaultResidueContactAnnotator extends AbstractFeatureProvider implements Annotator {
	/**
     * marks this residue as part of a sequence motif
     */
    public static final double RESIDUE_CONTACT_VALUE = 1.0;
    
    /**
     * some arbitrary value which should separate residues so no sequential neighbors interact
     */
    public static final int SEQUENTIAL_SEPARATION = 7;
    
    /**
     * for faster distance calculations
     */
    private static final double HELIX_HELIX_INTERACTION_CB_CUTOFF_SQUARED = HELIX_HELIX_INTERACTION_CB_CUTOFF * HELIX_HELIX_INTERACTION_CB_CUTOFF;
	
	public DefaultResidueContactAnnotator(FeatureExtractor featureExtractor, LogService logger, LinearAlgebra linearAlgebra, ModelConverter modelConverter) {
		super(featureExtractor, logger, linearAlgebra, modelConverter, new FeatureType[] { /*FeatureType.RESIDUE_CONTACTS*/ });
	}

	@Override
	protected void computeFeatureInternal(Protein protein) {
//        List<Interaction> residueInteractions = new ArrayList<>();
//        for(Chain chain : protein.chains) {
//	        for (int residueIndex1 = 0; residueIndex1 < chain.residues.size() - /*1*/SEQUENTIAL_SEPARATION; residueIndex1++) {
//	        	Residue residue1 = chain.residues.get(residueIndex1);
//	        	Atom cb1 = this.modelConverter.findAtomByName(residue1, BETA_CARBON_NAME);
//	            outer: for (int residueIndex2 = residueIndex1 + /*1*/SEQUENTIAL_SEPARATION; residueIndex2 < chain.residues.size(); residueIndex2++) {
//	            	Residue residue2 = chain.residues.get(residueIndex2);
//	            	// move on if sequential separation criterion is not met
////	            	if(Math.abs(residue1.getResidueNumber() - residue2.getResidueNumber()) < SEQUENTIAL_SEPARATION) {
////	            		continue;
////	            	}
//	            	
//	            	Atom cb2 = this.modelConverter.findAtomByName(residue2, BETA_CARBON_NAME);
//	                // interaction according to criterion 1 - cb-cb <6 A
//	                if (this.linearAlgebra.distanceFast(cb1.xyz, cb2.xyz) < HELIX_HELIX_INTERACTION_CB_CUTOFF_SQUARED) {
//	                    residueInteractions.add(createResidueResidueInteraction(residue1, residue2, InteractionType.CBETA_DISTANCE));
//	                    continue outer;
//	                }
//	
//	                // interaction according to criterion 2 - vdw-vdw <0.6 A
//	                for (int i = 0; i < residue1.atoms.size() - 1; i++) {
//	                    for (int j = i + 1; j < residue2.atoms.size(); j++) {
//	                        Atom a1 = residue1.atoms.get(i);
//	                        Atom a2 = residue2.atoms.get(j);
//	                        if (this.linearAlgebra.distance(a1.xyz, a2.xyz)) - lookUpVanDerWaalsRadius(a1)
//	                                - lookUpVanDerWaalsRadius(a2) < HELIX_HELIX_INTERACTION_VDW_CUTOFF) {
//	                            residueInteractions.add(createResidueResidueInteraction(residue1, residue2,
//	                                    InteractionType.VDW_DISTANCE));
//	                            continue outer;
//	                        }
//	                    }
//	                }
//	            }
//	        }
//        }
//        protein.interactions.addAll(residueInteractions);
//        
//        System.out.println("found " + residueInteractions.size() + " residue residue contacts");
	}
}
