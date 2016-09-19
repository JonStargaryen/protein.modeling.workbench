package de.bioforscher.pmw.feature.extractor.algorithm.impl;

import org.osgi.service.log.LogService;
import de.bioforscher.pmw.api.FeatureExtractor;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.feature.extractor.core.AbstractFeatureProvider;
import de.bioforscher.pmw.feature.extractor.core.Annotator;
import de.bioforscher.pmw.model.FeatureType;
import de.bioforscher.pmw.model.Protein;

@Deprecated
public class DefaultHelixAnnotator extends AbstractFeatureProvider implements Annotator {

	public DefaultHelixAnnotator(FeatureExtractor featureExtractor, LogService logger, LinearAlgebra linearAlgebra, ModelConverter modelConverter) {
		super(featureExtractor, logger, linearAlgebra, modelConverter, new FeatureType[] { null });
//		super(new FeatureType[] { FeatureType.HELIX_ANNOTATION });
	}

	@Override
	protected void computeFeatureInternal(Protein protein) {
//		annotateHelices(protein);
	}
	
//	/**
//    *
//    * @param protein
//    * @return true when at least 1 helix was found
//    */
//	private boolean annotateHelices(Protein protein) {
//       int helixCount = 0;
//       for (Chain c : protein.chains) {
//           // can't find tm-spanning helices in chains too short to give them a
//           // home
//           if (c.residues.size() < MINIMAL_LENGTH_OF_A_MEMBRANE_SPANNING_HELIX) {
//               continue;
//           }
//
//           int chainSpecificHelixCount = 0;
//           int helixStart = 0;
//           int helixEnd = 0;
//           boolean inHelix = false;
//           for (int i = 0; i < c.residues.size(); i++) {
//               Residue r = c.residues.get(i);
//               // helix starting
//               if (inHelix == false && r.getSecondaryStructure().equals(SecondaryStructure.HELIX)) {
//                   inHelix = true;
//                   helixStart = i;
//               }
//
//               // helix ending
//               if (inHelix == true
//                       && (!r.getSecondaryStructure().equals(SecondaryStructure.HELIX) || i == c.getSize() - 1)) {
//                   inHelix = false;
//                   helixEnd = i;
//
//                   // helix too short to be membrane spanning
//                   if (helixEnd - helixStart < MINIMAL_LENGTH_OF_A_MEMBRANE_SPANNING_HELIX) {
//                       continue;
//                   }
//
//                   // composed proposed helix object
//                   TransmembraneHelix helix = ProteinFactory.eINSTANCE.createTransmembraneHelix();
//                   helix.setName(c.getChainId() + "-TM-" + (chainSpecificHelixCount + 1));
//                   List<Residue> sublist = c.getResidues().subList(helixStart, helixEnd);
//                   helix.setSequence(extractSequence(sublist));
//                   helix.getResidues().addAll(sublist);
//
//                   // if the found helix is not in a TM region, it should not
//                   // considered
//                   if (noTransmembraneHelix(protein, helix)) {
//                       continue;
//                   }
//
//                   // found legit helix
//                   helixCount++;
//                   chainSpecificHelixCount++;
//                   // helix descriptors can only be computed when coordinates are present
//                   if(!protein.getReconstructionLevel().equals(ReconstructionLevel.NONE)) {
//                	   computeHelixDescriptors(protein, helix);
//                   }
//                   System.out.println("helix '" + helix.getName() + "' from " + sublist.get(0).getResidueNumber()
//                           + " to " + sublist.get(sublist.size() - 1) + " in chain " + c.getChainId()
//                           + "\n\tsequence: " + helix.getSequence() + "\n\tlength: " + helix.getLength()
//                           + " A" /*
//                                   * + "\n\tangle: " + helix.getTiltAngle() +
//                                   * "°"
//                                   */);
//                   protein.getHelices().add(helix);
//               }
//           }
//       }
//       System.out.println("found " + helixCount + " membrane spanning helices");
//       return helixCount > 0;
//   }
//
//   /**
//    * computes the length of a helix (the distance between to start residue's
//    * Ca and the end one's is measured) as well as the tilt angle relative to
//    * the membrane's plane
//    *
//    * @param protein
//    *            the containing protein
//    * @param helix
//    *            the helix to measure
//    */
//   private void computeHelixDescriptors(Protein protein, TransmembraneHelix helix) {
//       helix.setLength(computeDistance(retrieveAtomByName(helix.getResidues().get(0), ALPHA_CARBON_NAME),
//               retrieveAtomByName(helix.getResidues().get(helix.getResidues().size() - 1), ALPHA_CARBON_NAME)));
//       // TODO: implement tilt angle computation
//       helix.setTiltAngle(0.0);
//       // TODO: center of mass
//       // TODO: helical axis
//   }
}
