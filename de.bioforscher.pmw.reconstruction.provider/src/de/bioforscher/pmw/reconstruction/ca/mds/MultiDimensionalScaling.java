package de.bioforscher.pmw.reconstruction.ca.mds;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.biojava.nbio.structure.jama.EigenvalueDecomposition;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * The more or less (more less than more though) generic implementation of MDS.
 * <br />
 * This can be used to map data of some dimension (really high for gene
 * expression levels or quite low for distance maps) to some other dimension
 * which can obviously different and will usually be 2 or 3. This can be
 * exploited to placed atom coordinates for proteins in the 3D space merely by
 * knowing their interatomic distances. Neat!<br />
 * To actually use this for reconstruction purposes use the wrapper
 * {@link MultiDimensionalScalingReconstruction}.
 * 
 * @author S
 *
 */
public class MultiDimensionalScaling {
	/**
	 * the default dimension of the space the MDS will embed into<br />
	 * for sake of protein structure reconstruction, we aim at m=3
	 */
	public static final int DEFAULT_TARGET_DIMENSION = 3;

	/**
	 * the number of data points (determines the number of projected data points as well as the dimension of eigenvectors)
	 */
	private int numberOfDataPoints;

	/**
	 * the target dimension to which data points will be mapped (normally, 2 or
	 * 3) - also known as m
	 */
	private int targetDimension;

	/**
	 * the distance map provided as input, describes an arbitrary (yet
	 * symmetric) distance between an element i and j
	 */
	private Matrix distanceMap;

	/**
	 * the matrix A ('proximityMap' or D²) whereby each entry derived from the
	 * input distances a_ij = -0.5*d_ij^2 thus, dim(distanceMap) =
	 * dim(proximityMap)
	 */
	private Matrix proximityMap;

	/**
	 * the matrix B ('centeringMap') b_ij = a_ij - a_i* - a_j* + a_** a_i* - row
	 * average a_j* - column average a_** - overall average
	 */
	private Matrix centeringMap;

	/**
	 * the m (target dimension) normalized eigenvectors which can subsequently be used to compose the embedded data points
	 */
	private List<double[]> normalizedEigenvectors;

	/**
	 * embedding in the space of dimension m
	 */
	private List<double[]> embedding;

	/**
	 * transforms this distance map into an embedded in the 3D space
	 * 
	 * @param distanceMap
	 *            a symmetric distance map of dim <tt>nxn</tt> - so <tt>n</tt>
	 *            data points (e.g. atoms) are described by these distances
	 * @return a list of size <tt>n</tt> and each entry has dimensionality 3
	 *         (i.e. {@link SpatialCoordinate})
	 */
	public List<double[]> computeEmbedding(double[][] distanceMap) {
		// wrap the call, so the JAMA matrix object is not only way to use this
		// class
		return computeEmbedding(distanceMap, DEFAULT_TARGET_DIMENSION);
	}

	public List<double[]> computeEmbedding(double[][] distanceMap, int targetDimension) {
		return computeEmbedding(new Matrix(distanceMap), targetDimension);
	}

	public List<double[]> computeEmbedding(Matrix distanceMap, int targetDimension) {
		this.numberOfDataPoints = distanceMap.getRowDimension();
		this.targetDimension = targetDimension;
		if (this.targetDimension > this.numberOfDataPoints) {
			throw new IllegalArgumentException("target dimension must not exceed number of data points");
		}

		this.distanceMap = distanceMap;

		this.proximityMap = computeSquaredProximityMap(this.distanceMap);

		this.centeringMap = this.computeConfiguration(this.proximityMap);

		this.normalizedEigenvectors = new ArrayList<>();
		this.embedding = new ArrayList<>();
		EigenvalueDecomposition evd = this.centeringMap.eig();
		// we are looking for the m biggest eigenvalues - they are at the last
		// elements of the matrix
		Matrix eigenvectors = evd.getV();
		// Matrix eigenvalues = evd.getD();
		double[] eigenvalues = evd.getRealEigenvalues();
//		System.out.println(Arrays.toString(eigenvalues));
		Map<Integer, Double> eigenvalueMap = new HashMap<>();
		for (int eigenvalueIndex = 0; eigenvalueIndex < eigenvalues.length; eigenvalueIndex++) {
			eigenvalueMap.put(eigenvalueIndex, eigenvalues[eigenvalueIndex]);
		}
		List<Entry<Integer, Double>> sortedEigenvalues = entriesSortedByValues(eigenvalueMap).subList(0,
				targetDimension);

		// normalize eigenvectors
		for (Entry<Integer, Double> sortedEigenvalue : sortedEigenvalues) {
			if (sortedEigenvalue.getValue() <= 0) {
				throw new IllegalArgumentException("eigenvalue is negative: " + sortedEigenvalue.getValue());
			}
//			System.out.println("one of the largest eigenvalues has index '" + sortedEigenvalue.getKey() + "' and value '" + sortedEigenvalue.getValue() + "'");
			int index = this.numberOfDataPoints * sortedEigenvalue.getKey();
			double[] eigenvector = Arrays.copyOfRange(eigenvectors.getColumnPackedCopy(), index,
					index + this.numberOfDataPoints);
			this.normalizedEigenvectors.add(normalize(eigenvector, Math.sqrt(sortedEigenvalue.getValue())));
		}

		// compose embedded data points from normalized eigenvectors
		for (int dataPointIndex = 0; dataPointIndex < this.numberOfDataPoints; dataPointIndex++) {
			double[] dataPoint = new double[this.targetDimension];
			for (int dataPointDimension = 0; dataPointDimension < this.targetDimension; dataPointDimension++) {
				dataPoint[dataPointDimension] = this.normalizedEigenvectors.get(dataPointDimension)[dataPointIndex];
			}
			this.embedding.add(dataPoint);
		}

		return this.embedding;
	}

	// TODO: maybe move to commons
	private static <K, V extends Comparable<? super V>> List<Entry<K, V>> entriesSortedByValues(Map<K, V> map) {
		List<Entry<K, V>> sortedEntries = new ArrayList<Entry<K, V>>(map.entrySet());
//		Collections.sort(sortedEntries, new Comparator<Entry<K, V>>() {
//			@Override
//			public int compare(Entry<K, V> e1, Entry<K, V> e2) {
//				return e2.getValue().compareTo(e1.getValue());
//			}
//		});
		sortedEntries.sort((Entry<K, V> e1, Entry<K, V> e2)->e1.getValue().compareTo(e2.getValue()));

		return sortedEntries;
	}

	private double[] normalize(double[] vector, double normalizationFactor) {
		for (int dimension = 0; dimension < vector.length; dimension++) {
			vector[dimension] *= normalizationFactor;
		}
		return vector;
	}

	private Matrix computeConfiguration(Matrix proximityMap) {
		// TODO: what is the fastest way to map these values? naive,
		// native/System.arraycopy, stream?
		Matrix centeringMap = new Matrix(this.numberOfDataPoints, this.numberOfDataPoints);

		double[] rowAverage = new double[this.numberOfDataPoints];
		double[] columnAverage = new double[this.numberOfDataPoints];
		double overallAverage = 0;
		// assess rows and overall average
		for (int row = 0; row < this.numberOfDataPoints; row++) {
			double tempRowAverage = 0;
			for (int column = 0; column < this.numberOfDataPoints; column++) {
				double entry = proximityMap.get(row, column);
				tempRowAverage += entry;
				overallAverage += entry;
			}
			rowAverage[row] = tempRowAverage / this.numberOfDataPoints;
		}
		overallAverage /= this.numberOfDataPoints * this.numberOfDataPoints;

		// assess columns
		for (int column = 0; column < this.numberOfDataPoints; column++) {
			double tempColumnAverage = 0;
			for (int row = 0; row < this.numberOfDataPoints; row++) {
				tempColumnAverage += proximityMap.get(row, column);
			}
			columnAverage[column] = tempColumnAverage / this.numberOfDataPoints;
		}

		for (int row = 0; row < this.numberOfDataPoints; row++) {
			for (int column = 0; column < this.numberOfDataPoints; column++) {
				// b_ij = a_ij - a_i* - a_j* + a_**
				centeringMap.set(row, column,
						proximityMap.get(row, column) - rowAverage[row] - columnAverage[column] + overallAverage);
			}
		}

		return centeringMap;
	}

	private Matrix computeSquaredProximityMap(Matrix distanceMap) {
		Matrix proximityMap = new Matrix(this.numberOfDataPoints, this.numberOfDataPoints);

		for (int row = 0; row < this.numberOfDataPoints; row++) {
			for (int column = 0; column < this.numberOfDataPoints; column++) {
				// a_ij = -0.5*d_ij^2
				double entry = distanceMap.get(row, column);
				proximityMap.set(row, column, -0.5 * entry * entry);
			}
		}

		return proximityMap;
	}
}
