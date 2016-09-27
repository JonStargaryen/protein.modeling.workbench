package de.bioforscher.pmw.feature.extractor.provider;

import org.osgi.service.component.annotations.Activate;
import org.osgi.service.component.annotations.Component;
import org.osgi.service.component.annotations.Reference;
import org.osgi.service.log.LogService;
import de.bioforscher.pmw.api.FeatureExtractor;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.feature.extractor.core.AbstractFeatureProvider;
import de.bioforscher.pmw.feature.extractor.core.FeatureProvider;
import de.bioforscher.pmw.feature.extractor.factory.FeatureProviderAlgorithmFactory;
import de.bioforscher.pmw.model.FeatureType;
import de.bioforscher.pmw.model.Protein;

@Component(name = "de.bioforscher.pmw.feature.extractor")
public class FeatureExtractorImpl implements FeatureExtractor {
	/**
	 * service references
	 */
	@Reference
	private LinearAlgebra linearAlgebra;
	@Reference
	private ModelConverter modelConverter;
	@Reference
	private LogService logger;
	private FeatureProviderAlgorithmFactory factory;
	
	@Activate
	public void activate() {
		this.factory = new FeatureProviderAlgorithmFactory(this, this.logger, this.linearAlgebra, this.modelConverter);
	}
	
	@Override
	public void computeFeatures(Protein protein, FeatureType... featuresToCompute) {
		//TODO: added later-on, disabled for testing - however, this should probably be pivoted to a more fine-grained dependency on reconstruction levels 
//		boolean hasCoordinates = !protein.getReconstructionLevel().equals(ReconstructionLevel.NONE);
		for(FeatureType featureOption : featuresToCompute) {
			// retrieve suitable FeatureProvider
			FeatureProvider featureProvider = null;
			switch (featureOption) {
			case MOTIF_ANNOTATION:
				featureProvider = this.factory.createDefaultSequenceMotifAnnotator();
				break;
			case ACCESSIBLE_SURFACE_AREA:
				featureProvider = this.factory.createAccessibleSurfaceAreaCalculator();
				break;
			case SECONDARY_STRUCTURE:
				featureProvider = this.factory.createSecondaryStructureAnnotator();
				break;
			case MEMBRANE_TOPOLOGY:
				featureProvider = this.factory.createAnvilInstance();
				break;
			case INTERACTIONS:
				featureProvider = new AbstractFeatureProvider(this, this.logger, this.linearAlgebra, this.modelConverter, new FeatureType[] { FeatureType.INTERACTIONS }) {
					@Override
					protected void computeFeatureInternal(Protein protein) {
						// TODO Auto-generated method stub	
					}
				};
				break;
			default:
				//TODO this exception will not be propagate to the front-end at all
				throw new UnsupportedOperationException(featureOption.name() + " is not yet implemented");
			}
			
			// recursively generate required features
			FeatureType[] missingFeatures = featureProvider.checkForMissingRequirements(protein);
			computeFeatures(protein, missingFeatures);
			this.logger.log(LogService.LOG_INFO, "using " + featureProvider.getClass().getSimpleName() + " to generate " + featureOption);
			featureProvider.computeFeature(protein);
		}
	}
}
