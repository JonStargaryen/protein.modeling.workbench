package de.bioforscher.pmw.test;

import java.util.Arrays;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.osgi.framework.BundleContext;
import org.osgi.framework.FrameworkUtil;
import org.osgi.util.tracker.ServiceTracker;

import de.bioforscher.pmw.api.LinearAlgebra;

public class LinearAlgebraTest {
	private final BundleContext context = FrameworkUtil.getBundle(this.getClass()).getBundleContext();
	private LinearAlgebra linearAlgebra;
	private static final double[] VECTOR = new double[] { 101, 102, 103 };
	private static final double[] TRANSLATION_VECTOR = new double[] { 100, 100, 100 };
	private static final double[][] ROTATION_MATRIX = new double[][] {{ 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }};

	@Test
	public void shouldAddNumbers() {
		double result = this.linearAlgebra.add(VECTOR, VECTOR)[0];
		Assert.assertEquals(202.0, result, 0.0);
	}
	
	@Test
	public void shouldTransformVector() {
		double[] result = this.linearAlgebra.transform(VECTOR, TRANSLATION_VECTOR, ROTATION_MATRIX);
		System.out.println(Arrays.toString(result));
		
	}

	@Before
	public void setup() throws Exception {
		Assert.assertNotNull(this.context);
		this.linearAlgebra = getService(LinearAlgebra.class);
		Assert.assertNotNull(this.linearAlgebra);
	}

	private <T> T getService(Class<T> clazz) throws InterruptedException {
		ServiceTracker<T, T> st = new ServiceTracker<>(this.context, clazz, null);
		st.open();
		return st.waitForService(1000);
	}
}
