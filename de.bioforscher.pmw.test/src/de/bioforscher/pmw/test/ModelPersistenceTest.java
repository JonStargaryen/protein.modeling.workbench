package de.bioforscher.pmw.test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.osgi.framework.BundleContext;
import org.osgi.framework.FrameworkUtil;
import org.osgi.util.tracker.ServiceTracker;

import de.bioforscher.pmw.api.ModelPersistence;

public class ModelPersistenceTest {
	private final BundleContext context = FrameworkUtil.getBundle(this.getClass()).getBundleContext();
	private ModelPersistence modelPersistence;

	@Test
	public void shouldAlignmentFragments() {
		System.out.println("ja! " + this.modelPersistence);
	}
	
	@Before
	public void setup() throws Exception {
		Assert.assertNotNull(this.context);
		this.modelPersistence = getService(ModelPersistence.class);
		Assert.assertNotNull(this.modelPersistence);
	}

	private <T> T getService(Class<T> clazz) throws InterruptedException {
		ServiceTracker<T, T> st = new ServiceTracker<>(this.context, clazz, null);
		System.out.println(st);
		st.open();
		return st.waitForService(10000);
	}
}
