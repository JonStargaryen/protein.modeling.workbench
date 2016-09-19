package de.bioforscher.pmw.feature.extractor.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.DoubleSummaryStatistics;
import java.util.List;
import java.util.stream.Stream;

import org.slf4j.Logger;
import de.bioforscher.pmw.api.FeatureExtractor;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.model.Atom;
import de.bioforscher.pmw.model.FeatureType;
import de.bioforscher.pmw.model.Protein;

public abstract class AbstractFeatureProvider implements FeatureProvider {

	protected FeatureExtractor featureExtractor;
	protected Logger logger;
	protected ModelConverter modelConverter;
	protected LinearAlgebra linearAlgebra;
	private final FeatureType[] PROVIDED_FEATURES;
	private final FeatureType[] REQUIRED_FEATURES;

	public AbstractFeatureProvider(FeatureExtractor featureExtractor, Logger logger, LinearAlgebra linearAlgebra,
			ModelConverter modelConverter, FeatureType[] providedFeatures, FeatureType... requiredFeatures) {
		this.featureExtractor = featureExtractor;
		this.logger = logger;
		this.linearAlgebra = linearAlgebra;
		this.modelConverter = modelConverter;
		this.PROVIDED_FEATURES = providedFeatures;
		this.REQUIRED_FEATURES = requiredFeatures;
	}

	protected void assignBaseline(Protein protein, FeatureType[] featureType) {
		this.modelConverter.getResidues(protein).forEach(r -> {
			Stream.of(featureType).forEach(ft -> {
				r.features.put(ft.name(), wrapInArray(0.0));
			});
		});
	}

	protected double[] wrapInArray(double value) {
		return new double[] { value, 0.0 };
	}

	// protected Feature wrapInFeature(double value) {
	// Feature feature = new Feature();
	// feature.value = value;
	// return feature;
	// }

	@Override
	public FeatureType[] getRequiredFeatures() {
		return REQUIRED_FEATURES;
	}

	@Override
	public FeatureType[] getProvidedFeatures() {
		return PROVIDED_FEATURES;
	}

	// /**
	// * checks whether this helix is TM - this is realized by determining
	// whether
	// * topology information is present in the first place - later on,
	// *
	// * @return
	// */
	// protected boolean noTransmembraneHelix(Protein protein,
	// TransmembraneHelix helix) {
	// // when no topology information is present, every helix is considered a
	// // valid helix
	// if (!protein.availableFeatures.contains(FeatureType.MEMBRANE_TOPOLOGY)) {
	// return false;
	// }
	//
	// // TODO: implement
	// return false;
	// }

	/**
	 * uses BioJava's information to look up the Van-der-Waals-radii of certain
	 * elements
	 *
	 * @param atom
	 *            the atom to measure
	 * @return the vdw-radius in A
	 */
	protected double lookUpVanDerWaalsRadius(Atom atom) {
		return Element.valueOfIgnoreCase(atom.element).getVDWRadius();
	}
	
//	protected Optional<Atom> retrieveAtomByName(Residue residue, String atomName) {
//		return residue.atoms.stream().filter(r -> r.name.equals(atomName)).findFirst();
//	}

	// protected ResidueResidueInteraction
	// createResidueResidueInteraction(Residue residue1, Residue residue2,
	// ResidueResidueInteractionType type) {
	// // create global residue-residue interaction object
	// ResidueResidueInteraction residueResidueInteraction =
	// ProteinFactory.eINSTANCE
	// .createResidueResidueInteraction();
	// residueResidueInteraction.getInteractingResidues().addAll(Arrays.asList(new
	// Residue[] { residue1, residue2 }));
	// residueResidueInteraction.setInteractionType(type);
	// // create local annotation
	// replaceFeature(residue1, FeatureType.RESIDUE_CONTACTS,
	// DefaultResidueContactAnnotator.RESIDUE_CONTACT_VALUE);
	// replaceFeature(residue2, FeatureType.RESIDUE_CONTACTS,
	// DefaultResidueContactAnnotator.RESIDUE_CONTACT_VALUE);
	//
	//// System.out.println(residueResidueInteraction.getInteractingResidues().get(0)
	// + " <=> "
	//// + residueResidueInteraction.getInteractingResidues().get(1) + " : "
	//// + residueResidueInteraction.getInteractionType());
	// return residueResidueInteraction;
	// }

	@Override
	public FeatureType[] checkForMissingRequirements(Protein protein) {
		List<FeatureType> missingFeatures = new ArrayList<>();
		for (FeatureType featureType : REQUIRED_FEATURES) {
			if (!protein.availableFeatures.contains(featureType))
				missingFeatures.add(featureType);
		}
		return missingFeatures.toArray(new FeatureType[0]);
	}

	@Override
	public void computeFeature(Protein protein) {
		// fail when required features are missing
		if (checkForMissingRequirements(protein).length != 0) {
			throw new RuntimeException(getClass().getSimpleName() + " depends on the features: "
					+ Arrays.deepToString(checkForMissingRequirements(protein)) + " to be present");
		}

		// set baseline so visualization is not messed up
		assignBaseline(protein, PROVIDED_FEATURES);

		// delegate to concrete implementation
		computeFeatureInternal(protein);

		normalizeValues(protein);

		setFeatureFlags(protein, PROVIDED_FEATURES);
	}

	protected void normalizeValues(Protein protein) {
		for (FeatureType featureType : this.getProvidedFeatures()) {
			DoubleSummaryStatistics summary = this.modelConverter.getResidues(protein).stream()
					.map(r -> r.features.get(featureType.name())).mapToDouble(r -> r[0]).summaryStatistics();
			final double min = summary.getMin();
			final double max = summary.getMax();
//			System.out.println(featureType.name() + " - min : " + min + " - max : " + max + " - sum : " + summary.getSum());

			this.modelConverter.getResidues(protein).stream().map(r -> r.features.get(featureType.name()))
					.forEach(f -> {
						f[1] = (f[0] - min) / (max - min);
						// some safety net, so no NaNs are propagated to the
						// front-end
						if (Double.isNaN(f[1])) {
							f[1] = 0.0;
						}
					});
		}
	}

	protected abstract void computeFeatureInternal(Protein protein);

	protected void setFeatureFlags(Protein protein, FeatureType... featureTypes) {
		protein.availableFeatures.addAll(Arrays.asList(featureTypes));
	}
}
