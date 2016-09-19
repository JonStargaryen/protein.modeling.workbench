package de.bioforscher.pmw.feature.extractor.factory;

import org.slf4j.Logger;

import de.bioforscher.pmw.api.FeatureExtractor;
import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.feature.extractor.algorithm.impl.DSSP;
import de.bioforscher.pmw.feature.extractor.algorithm.impl.DefaultHelixAnnotator;
import de.bioforscher.pmw.feature.extractor.algorithm.impl.DefaultHelixInteractionAnnotator;
import de.bioforscher.pmw.feature.extractor.algorithm.impl.DefaultResidueContactAnnotator;
import de.bioforscher.pmw.feature.extractor.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.pmw.feature.extractor.core.FeatureProvider;
import de.bioforscher.pmw.feature.extractor.motif.DefaultSequenceMotifAnnotator;
import de.bioforscher.pmw.feature.extractor.sse.SecondaryStructureElementAnnotator;
import de.bioforscher.pmw.feature.extractor.topology.ANVIL;

@SuppressWarnings("deprecation")
public class FeatureProviderAlgorithmFactory {
	private FeatureExtractor featureExtractor;
	private Logger logger;
	private LinearAlgebra linearAlgebra;
	private ModelConverter modelConverter;
	
	public FeatureProviderAlgorithmFactory(FeatureExtractor featureExtractor, Logger logger, LinearAlgebra linearAlgebra, ModelConverter modelConverter) {
		this.featureExtractor = featureExtractor;
		this.logger = logger;
		this.linearAlgebra = linearAlgebra;
		this.modelConverter = modelConverter;
	}
	
	public FeatureProvider createSecondaryStructureAnnotator() {
		return new SecondaryStructureElementAnnotator(this.featureExtractor, this.logger, this.linearAlgebra, this.modelConverter);
	}
	
	public FeatureProvider createAccessibleSurfaceAreaCalculator() {
		return new AccessibleSurfaceAreaCalculator(this.featureExtractor, this.logger, this.linearAlgebra, this.modelConverter);
	}
	
	public FeatureProvider createDefaultHelixAnnotator() {
		return new DefaultHelixAnnotator(this.featureExtractor, this.logger, this.linearAlgebra, this.modelConverter);
	}
	
	public FeatureProvider createDefaultHelixInteractionAnnotator() {
		return new DefaultHelixInteractionAnnotator(this.featureExtractor, this.logger, this.linearAlgebra, this.modelConverter);
	}
	
	public FeatureProvider createDefaultResidueContactAnnotator() {
		return new DefaultResidueContactAnnotator(this.featureExtractor, this.logger, this.linearAlgebra, this.modelConverter);
	}
	
	public FeatureProvider createDefaultSequenceMotifAnnotator() {
		return new DefaultSequenceMotifAnnotator(this.featureExtractor, this.logger, this.linearAlgebra, this.modelConverter);
	}
	
	public FeatureProvider createDSSPInstance() {
		return new DSSP(this.featureExtractor, this.logger, this.linearAlgebra, this.modelConverter);
	}

	public FeatureProvider createAnvilInstance() {
		return new ANVIL(this.featureExtractor, this.logger, this.linearAlgebra, this.modelConverter);
	}
}
