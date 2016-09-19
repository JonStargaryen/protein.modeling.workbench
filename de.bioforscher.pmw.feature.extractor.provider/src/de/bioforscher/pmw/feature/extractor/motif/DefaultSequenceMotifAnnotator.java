package de.bioforscher.pmw.feature.extractor.motif;

import java.util.List;
import java.util.stream.Collectors;

import org.osgi.service.log.LogService;
import de.bioforscher.pmw.api.FeatureExtractor;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.feature.extractor.core.AbstractFeatureProvider;
import de.bioforscher.pmw.feature.extractor.core.Annotator;
import de.bioforscher.pmw.model.Chain;
import de.bioforscher.pmw.model.DefinedMotif;
import de.bioforscher.pmw.model.FeatureType;
import de.bioforscher.pmw.model.Motif;
import de.bioforscher.pmw.model.Protein;
import de.bioforscher.pmw.model.Residue;
import de.bioforscher.pmw.model.Topology;

public class DefaultSequenceMotifAnnotator extends AbstractFeatureProvider implements Annotator {
    
	public DefaultSequenceMotifAnnotator(FeatureExtractor featureExtractor, LogService logger, LinearAlgebra linearAlgebra, ModelConverter modelConverter) {
		super(featureExtractor, logger, linearAlgebra, modelConverter, new FeatureType[] {FeatureType.MOTIF_ANNOTATION});
	}

	/**
     * marks this residue as part of a sequence motif
     */
    public static final double SEQUENCE_MOTIF_VALUE = 1.0;
    
	@Override
	protected void computeFeatureInternal(Protein protein) {
		for (Chain chain : protein.chains) {
            int chainLength = chain.residues.size();
            for (int resNum = 0; resNum < chainLength; resNum++) {
                Residue startResidue = chain.residues.get(resNum);
                // even though we express 1-letter-codes as Strings, excessive matching is probably much faster when done with chars
                char startAminoAcid = this.modelConverter.convertToOneLetterCode(startResidue.aminoAcid).charAt(0);
                for (DefinedMotif candidate : DefinedMotif.values()) {
                    // get motif length
                    int motifLength = Integer.parseInt(candidate.name().substring(2));
                    char motifStart = candidate.name().charAt(0);
                    char motifEnd = candidate.name().charAt(1);

                    // chain not long enough to cover the proposed motif
                    if (resNum + motifLength >= chainLength) {
                        continue;
                    }

                    // start amino acid names does not match the motif
                    if (startAminoAcid != motifStart) {
                        continue;
                    }

                    Residue endResidue = chain.residues.get(resNum + motifLength);

                    // end amino acid does not match
                    if (this.modelConverter.convertToOneLetterCode(endResidue.aminoAcid).charAt(0) != motifEnd) {
                        continue;
                    }

                    Motif motif = new Motif();
                    List<Residue> sublist = chain.residues.subList(resNum, resNum + motifLength + 1);
                    motif.startResidueId = startResidue.residueId;
                    motif.endResidueId = endResidue.residueId;
                    markResiduesAsPartOfSequenceMotif(sublist);
                    // when membrane topology information is available, annotate the motifs topology
                    if(protein.availableFeatures.contains(FeatureType.MEMBRANE_TOPOLOGY)) {
                    	double startResidueTopology = startResidue.features.get(FeatureType.MEMBRANE_TOPOLOGY)[0];
                    	motif.topology = startResidueTopology == endResidue.features.get(FeatureType.MEMBRANE_TOPOLOGY)[0]
                            ? startResidueTopology : Topology.TRANSITION.ordinal();
                    }
                    motif.sequence = extractSequence(sublist);
                    motif.definition = candidate;
                    protein.motifs.add(motif);
                }
            }
        }
	}
	
	private String extractSequence(List<Residue> sublist) {
		return sublist.stream().map(r -> this.modelConverter.convertToOneLetterCode(r.aminoAcid)).collect(Collectors.joining());
	}
	
	private void markResiduesAsPartOfSequenceMotif(List<Residue> sublist) {
		sublist.forEach(r -> {
	        r.features.put(FeatureType.MOTIF_ANNOTATION.name(), wrapInArray(SEQUENCE_MOTIF_VALUE));
		});		
	}
}
