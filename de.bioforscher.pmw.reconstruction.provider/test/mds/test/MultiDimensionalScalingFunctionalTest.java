package mds.test;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

import de.bioforscher.pmw.reconstruction.ca.mds.MultiDimensionalScaling;
import junit.framework.TestCase;

public class MultiDimensionalScalingFunctionalTest extends TestCase {
	private static final String POCKET_ALIGN_CSV_PATH = "mds/test/rmsd_matrix_scheme1.csv";
	
	public void testMultiDimensionalScaling() {
		try {
			MultiDimensionalScaling mds = new MultiDimensionalScaling();
			double[][] dataPoints = parseDistancesFromPocketAlignCSV(getResourceAsFilepath(POCKET_ALIGN_CSV_PATH));
			List<String> dataLabels = parseLabels(getResourceAsFilepath(POCKET_ALIGN_CSV_PATH));
			List<double[]> embeddedDataPoints = mds.computeEmbedding(dataPoints, 3);
			for(int i = 0; i < dataLabels.size(); i++) {
				String s = Arrays.toString(embeddedDataPoints.get(i));
				System.out.println(dataLabels.get(i) + "," + s.replace("[", "").replace("]", "").replace(" ", ""));
			}
		} catch (IOException e) {
			fail("failed with IOException: " + e.getLocalizedMessage());
		}
	}
	
	private String getResourceAsFilepath(String filename) {
		ClassLoader ccl = Thread.currentThread().getContextClassLoader();
		Objects.requireNonNull(ccl);
		URL resource = ccl.getResource(filename);
		Objects.requireNonNull(resource);
		return resource.getPath();
	}
	
	private List<String> parseLabels(String filepath) throws IOException {
		String line = Files.readAllLines(new File(filepath).toPath()).get(0);		
		// omit first ,
		return Arrays.asList(line.substring(1, line.length()).split(",")).stream().collect(Collectors.toList());
	}
	
	private double[][] parseDistancesFromPocketAlignCSV(String filepath) throws IOException {
		List<String> lines = Files.readAllLines(new File(filepath).toPath());
		double[][] data = new double[lines.size() - 1][];
		
		// skip first line
		for(int i = 0; i < lines.size(); i++) {
			if(i == 0) {
				continue;
			}
			String[] tmpLine = lines.get(i).split(",");
			double[] tmpData = new double[lines.size() - 1];
			for(int j = 0; j < lines.size(); j++) {
				if(j == 0) {
					continue;
				}
				if(i == j) {
					tmpData[j - 1] = 0.0;
					continue;
				}
				tmpData[j - 1] = Double.valueOf(tmpLine[j]);
			}
			data[i - 1] = tmpData;
		}
		
		return data;
	}
}
