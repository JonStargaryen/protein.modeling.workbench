package de.bioforscher.pmw.test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.osgi.framework.BundleContext;
import org.osgi.framework.FrameworkUtil;
import org.osgi.util.tracker.ServiceTracker;

import de.bioforscher.pmw.api.AlignmentService;

public class AlignmentServiceTest {
	private final BundleContext context = FrameworkUtil.getBundle(this.getClass()).getBundleContext();
	private AlignmentService alignmentService;

	@Test
	public void shouldAlignmentFragments() {
		System.out.println("ja! " + this.alignmentService);
	}
	
	@Before
	public void setup() throws Exception {
		Assert.assertNotNull(this.context);
		this.alignmentService = getService(AlignmentService.class);
		Assert.assertNotNull(this.alignmentService);
	}

	private <T> T getService(Class<T> clazz) throws InterruptedException {
		ServiceTracker<T, T> st = new ServiceTracker<>(this.context, clazz, null);
		st.open();
		return st.waitForService(10000);
	}
}
