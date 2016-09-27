package de.bioforscher.pmw.api;

/**
 * A service which implements a potpourri of linear algebraic calculations. All operations are low-level on the level of <code>double[]</code> and, by contract,
 * do never happen in place (which means they do not modify any of their arguments, but rather return a new instance to wrap the result).<br />
 * Access is provided by strictly low-level <code>double[]</code> or to process model entries such as {@link Atom} use their <code>xyz</code> field.<br />
 * There is no real point in wrapping this in an OSGi service (in contrast to a class with static methods), maybe this even degenerates performance significantly as all these functions tend to be invoked numerous time during computation/reconstruction steps.
 * @author S
 *
 */
public interface LinearAlgebra {
	/**
	 * add 2 vectors
	 * @param v1
	 * @param v2
	 * @return the vector sum
	 */
	double[] add(double[] v1, double[] v2);
	
	/**
	 * Returns the angle in radians between this vector and the vector
	 * parameter; the return value is constrained to the range [0,PI].
	 * @return the angle in radians in the range [0,PI]
	 */
	double angle(double[] v1, double[] v2);
	
	/**
	 * Calculates the distance between 2 points. Calls {@link LinearAlgebra#distanceFast(double[], double[])} and takes the square root of the result.
	 * You are encouraged to use {@link LinearAlgebra#distanceFast(double[], double[])} where appropriate (e.g. to sort distances or to check some distance relative to some threshold: in that case square the threshold value once).
	 * @param v1
	 * @param v2
	 * @return
	 */
	double distance(double[] v1, double[] v2);
	
	/**
	 * specific for BBQ and PULCHRA implementation - based on this term, the distance between residue 1 and 4 of four consecutive residues can be computed<br />
	 * the sign of the number will indicate the handedness of this particular conformation
	 * @see LinearAlgebra#isLeftHandedTwist(double[], double[], double[], double[])
	 * @param v1 four
	 * @param v2 consecutive
	 * @param v3 backbone
	 * @param v4 CA atoms
	 * @return
	 */
	double distance14(double[] v1, double[] v2, double[] v3, double[] v4);
	
	/**
	 * computes the rotation matrix needed to place the quadrilateral scaffold correctly
	 * @param v1
	 * @param v2
	 * @param v3
	 * @return
	 */
	double[][] rotation(double[] v1, double[] v2, double[] v3);
	
	/**
	 * the squared distance between 2 points
	 * @param v1
	 * @param v2
	 * @return
	 * @see LinearAlgebra#distance(double[], double[])
	 */
	double distanceFast(double[] v1, double[] v2);
	
	/**
	 * divides a vector by a scalar
	 * @param v
	 * @param scalar
	 * @return
	 */
	double[] divide(double[] v, double scalar);
	
	/**
	 * the scalar product of 2 vectors
	 * @param v1
	 * @param v2
	 * @return
	 */
	double dotProduct(double[] v1, double[] v2);
	
	/**
	 * returns <code>true</code> iff this quadrilateral is left-handed
	 * @param ca1
	 * @param ca2
	 * @param ca3
	 * @param ca4
	 * @return true if these coordinates dictate a left-handed conformation
	 */
	boolean isLeftHandedTwist(double[] v1, double[] v2, double[] v3, double[] v4);
	
	/**
	 * multiplies each element of the given vector with a scalar
	 * @param v the vector
	 * @param scalar
	 */
	double[] multiply(double[] v, double scalar);
	
	double[] multiply(double[] v, double[][] m);
	
	/**
	 * the length of this vector
	 * @param v
	 * @return
	 */
	double norm(double[] v);
	
	/**
	 * normalizes a vector a.k.a. dividing each component by the length of the vector
	 * @param v the input vector
	 * @return a new <code>double[]</code> instance containing the normalized vector
	 */
	double[] normalize(double[] v);

	/**
	 * rototranslates a vector utilizing a given translation vector and rotation matrix
	 * @param vector the vector to be rototranslated
	 * @param translation the vector describing the wanted translation
	 * @param rotation the matrix describing the wanted rotation
	 * @return the manipulated coordinates of the input vector
	 */
	double[] transform(double[] vector, double[] translation, double[][] rotation);
	
	/**
	 * subtracts 2 vectors
	 * @param v1
	 * @param v2
	 * @return the difference vector
	 */
	double[] subtract(double[] v1, double[] v2);
	
	/**
	 * computes the vector product between 2 vectors
	 * @param v1
	 * @param v2
	 * @return
	 */
	double[] vectorProduct(double[] v1, double[] v2);

	/**
	 * Calculate the torsion angle, i.e. the angle between the normal vectors of
	 * the two plains a-b-c and b-c-d. See http://en.wikipedia.org/wiki/Dihedral_angle
	 * 
	 * @return the torsion angle in degrees, in range +-[0,180]. If either first
	 *         3 or last 3 atoms are colinear then torsion angle is not defined
	 *         and NaN is returned
	 */
	double torsionAngle(double[] v1, double[] v2, double[] v3, double[] v4);
}
