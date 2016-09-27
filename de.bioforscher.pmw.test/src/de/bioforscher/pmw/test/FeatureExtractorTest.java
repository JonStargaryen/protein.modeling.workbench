package de.bioforscher.pmw.test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.osgi.framework.BundleContext;
import org.osgi.framework.FrameworkUtil;
import org.osgi.util.tracker.ServiceTracker;

import de.bioforscher.pmw.api.FeatureExtractor;

public class FeatureExtractorTest {
	private final BundleContext context = FrameworkUtil.getBundle(this.getClass()).getBundleContext();
	private FeatureExtractor featureExtractor;

	@Test
	public void shouldAlignmentFragments() {
		System.out.println("ja! " + this.featureExtractor);
	}
	
	@Before
	public void setup() throws Exception {
		Assert.assertNotNull(this.context);
		this.featureExtractor = getService(FeatureExtractor.class);
		Assert.assertNotNull(this.featureExtractor);
	}

	private <T> T getService(Class<T> clazz) throws InterruptedException {
		ServiceTracker<T, T> st = new ServiceTracker<>(this.context, clazz, null);
		System.out.println(st);
		st.open();
		return st.waitForService(1000);
	}
}
