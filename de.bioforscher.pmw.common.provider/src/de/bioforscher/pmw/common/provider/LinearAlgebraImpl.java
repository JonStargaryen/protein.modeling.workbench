package de.bioforscher.pmw.common.provider;

import org.osgi.service.component.annotations.Component;
import org.osgi.service.component.annotations.Reference;

import de.bioforscher.pmw.api.LinearAlgebra;
import de.bioforscher.pmw.api.ModelConverter;

@Component(name = "de.bioforscher.pmw.linear.algebra")
public class LinearAlgebraImpl implements LinearAlgebra {
	@Reference
	private ModelConverter modelConverter;
	
	@Override
	public double[] add(double[] v1, double[] v2) {
		return new double[] { v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2] };
	}

	@Override
	public double angle(double[] v1, double[] v2) {
		double vDot = dotProduct(v1, v2) / (norm(v1) * norm(v2));
		if (vDot < -1.0) {
			vDot = -1.0;
		}
		if (vDot > 1.0) {
			vDot = 1.0;
		}
		return Math.toDegrees((double) Math.acos(vDot));
	}

	@Override
	public double distance(double[] v1, double[] v2) {
		return Math.sqrt(distanceFast(v1, v2));
	}
	@Override
	public double distance14(double[] v1, double[] v2, double[] v3, double[] v4) {
		double d14 = distance(v1, v4);
	    if(isLeftHandedTwist(v1, v2, v3, v4)) {
			d14 = -d14;
		}
	    return d14;
	}
	
	@Override
	public double distanceFast(double[] v1, double[] v2) {
		return (v1[0] - v2[0]) * (v1[0] - v2[0]) +
				(v1[1] - v2[1]) * (v1[1] - v2[1]) +
				(v1[2] - v2[2]) * (v1[2] - v2[2]);
	}

	@Override
	public double[] divide(double[] v, double scalar) {
		return new double[] { v[0] / scalar, v[1] / scalar, v[2] / scalar};
	}

	@Override
	public double dotProduct(double[] v1, double[] v2) {
		return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	}

	@Override
	public boolean isLeftHandedTwist(double[] v1, double[] v2, double[] v3, double[] v4) {
		double[] d21 = subtract(v2, v1);
		double[] d31 = subtract(v3, v1);
		double[] d41 = subtract(v4, v1);
		
		return dotProduct(d21, vectorProduct(d31, d41)) < 0.0;
	}

	@Override
	public double[] multiply(double[] v, double scalar) {
		return new double[] { v[0] * scalar, v[1] * scalar, v[2] * scalar};
	}

	@Override
	public double norm(double[] v) {
		return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	}

	@Override
	public double[] normalize(double[] v) {
		double length = norm(v);
		return new double[] { v[0] / length, v[1] / length, v[2] / length };
	}

	@Override
	public double[][] rotation(double[] v1, double[] v2, double[] v3) {
		double[] difference21 = normalize(subtract(v2, v1));
		double[] difference23 = normalize(subtract(v2, v3));
	    double[] difference13m = normalize(subtract(difference21, difference23));
	    double[] difference13p = normalize(add(difference21, difference23));
	    
	    return new double[][] { difference13m, difference13p,
	    	{ difference13m[1] * difference13p[2] - difference13m[2] * difference13p[1],
	    		difference13m[2] * difference13p[0] - difference13m[0] * difference13p[2],
	    		difference13m[0] * difference13p[1] - difference13m[1] * difference13p[0]
	    	}};
	}

	@Override
	public double[] rototranslate(double[] vector, double[] translation, double[][] rotation) {
	    double[] result = new double[3];
		double oldX = vector[0];
	    double oldY = vector[1];
	    double oldZ = vector[2];
	    result[0] = (oldX * rotation[0][0] + oldY * rotation[1][0] + oldZ * rotation[2][0]) + translation[0];
	    result[1] = (oldX * rotation[0][1] + oldY * rotation[1][1] + oldZ * rotation[2][1]) + translation[1];
	    result[2] = (oldX * rotation[0][2] + oldY * rotation[1][2] + oldZ * rotation[2][2]) + translation[2];
		return result;
	}

	@Override
	public double[] subtract(double[] v1, double[] v2) {
		return new double[] { v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2] };
	}

	@Override
	public double torsionAngle(double[] v1, double[] v2, double[] v3, double[] v4) {
		double[] ab = subtract(v1, v2);
		double[] cb = subtract(v3, v2);
		double[] bc = subtract(v2, v3);
		double[] dc = subtract(v4, v3);

		double[] abc = vectorProduct(ab, cb);
		double[] bcd = vectorProduct(bc, dc);

		double angl = angle(abc, bcd);

		/* calc the sign: */
		double[] vecprod = vectorProduct(abc, bcd);
		double val = dotProduct(cb, vecprod);
		if (val < 0.0)
			angl = -angl;

		return angl;
	}

	@Override
	public double[] vectorProduct(double[] v1, double[] v2) {
		return new double[] { v1[1] * v2[2] - v1[2] * v2[1],
				v1[2] * v2[0] - v1[0] * v2[2],
				v1[0] * v2[1] - v1[1] * v2[0]
		};
	}
}
