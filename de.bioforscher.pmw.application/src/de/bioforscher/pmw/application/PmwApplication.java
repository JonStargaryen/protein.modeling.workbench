package de.bioforscher.pmw.application;

import java.util.Base64;
import org.osgi.dto.DTO;
import org.osgi.service.component.annotations.Component;
import org.osgi.service.component.annotations.Reference;
import org.osgi.service.log.LogReaderService;
import org.osgi.service.log.LogService;
import org.slf4j.Logger;

import de.bioforscher.pmw.api.FeatureExtractor;
import de.bioforscher.pmw.api.ModelConverter;
import de.bioforscher.pmw.api.ModelPersistence;
import de.bioforscher.pmw.api.ReconstructionService;
import de.bioforscher.pmw.model.APIConstants;
import de.bioforscher.pmw.model.FeatureType;
import de.bioforscher.pmw.model.Project;
import de.bioforscher.pmw.model.ReconstructionLevel;
import osgi.enroute.configurer.api.RequireConfigurerExtender;
import osgi.enroute.google.angular.capabilities.RequireAngularWebResource;
import osgi.enroute.logger.api.RequireLoggerImplementation;
import osgi.enroute.rest.api.REST;
import osgi.enroute.rest.api.RESTRequest;
import osgi.enroute.twitter.bootstrap.capabilities.RequireBootstrapWebResource;
import osgi.enroute.webserver.capabilities.RequireWebServerExtender;

@RequireAngularWebResource(resource={"angular.js","angular-resource.js", "angular-route.js"}, priority=1000)
@RequireBootstrapWebResource(resource="css/bootstrap.css")
@RequireWebServerExtender
@RequireConfigurerExtender
@RequireLoggerImplementation
@Component(name="de.bioforscher.pmw")
public class PmwApplication implements REST {
	/** the flag indicating that a file was uploaded (with need to be base64-encoded to be processed by the REST-interface) - this value can also be used to decode the submitted string as e.g. the header information has be removed */
	private static final String FILE_ENCODING_FLAG = "base64,";
	private static final String FEATURE_CONTEXT = "feature";
	private static final String RECONSTRUCTION_CONTEXT = "reconstruction";
	private static final APIConstants SETTINGS = new APIConstants();
	/**
	 * TODO move these to config
	 */
	private static final boolean SUPPRESS_FRAMEWORK_MESSAGES = true;
	private static final int LOG_LEVEL = LogService.LOG_DEBUG;
	
	@Reference
	private FeatureExtractor featureExtractorService;
	@Reference
	private Logger logService;
	@Reference
	private ModelConverter modelConverterService;
	@Reference
	private ModelPersistence modelPersistenceService;
	@Reference
	private ReconstructionService reconstructionService;

	//TODO implement: some 'delta' function would be nice, so not the whole model has to be transfered but rather only the model's changes
	
	/*
	 * the interface to retrieve already persisted projects from the backend
	 */
	
	/**
	 * returns the requested project by querying the persistence service for a job with that ID
	 * @param id the UUID of the job requested
	 * @return the JSONified {@link ModelingProject}
	 * @throws Exception thrown upon not finding a corresponding job
	 */
	public Project getProject(RESTRequest request, String id) throws Exception {
		return this.modelPersistenceService.retrieve(id);
	}
	
	/*
	 * the interface to request calculations from the back-end
	 */
	
	public static class Computation extends DTO {
		public String projectId;
		public String context;
		public int value;
	}
	
	public interface ComputationRequest extends RESTRequest {
		public Computation _body();
	}
	
	public String postCalculation(ComputationRequest request) throws Exception {
		Project project = this.modelPersistenceService.retrieve(request._body().projectId);
		int value = request._body().value;
		String context = request._body().context;
		if(context.equals(FEATURE_CONTEXT)) {
			this.featureExtractorService.computeFeatures(project.proteins.get(0), FeatureType.values()[value]);
			this.modelPersistenceService.update(project);
			return FeatureType.values()[value].name();
		}
		if(context.equals(RECONSTRUCTION_CONTEXT)) {
			this.reconstructionService.reconstruct(project.proteins.get(0), ReconstructionLevel.values()[value]);
			this.modelPersistenceService.update(project);
			return ReconstructionLevel.values()[value].name();
		}
		
		throw new IllegalArgumentException("context not known");
	}
	
	/*
	 * the interface to create new ModelingJobs
	 */
	
	public interface CreateProjectRequest extends RESTRequest {
		public String _body();
	}
	
	/**
	 * the entry point to create {@link ModelingProject}s<br />
	 * this is done by querying the provided REST-interface<br />
	 * either, the user can:
	 * <ul>
	 * <li>provide a protein sequence as String and, thus, create a protein without</li>
	 * <li>upload a existing PDB file and, thus, parse all available information and investigate it</li>
	 * </ul>
	 * <br />
	 * the provided information will be wrapped into a {@link Protein} which itself will be wrapped into a modeling project
	 * @param input either a protein sequence or a uploaded, base64-encoded PDB file
	 * @return the UUID of the created project, so the front-end can retrieve it later on
	 * @throws Exception 
	 */
	public String postProject(CreateProjectRequest request) throws Exception {
		// create project with the specified inputs
		Project project;
		String input = request._body();
		if(input.contains(FILE_ENCODING_FLAG)) {
			// file uploaded - parse it into a protein object
			// strip encoding header as well as JSON-tail from the string
			input = input.substring(input.indexOf(FILE_ENCODING_FLAG) + FILE_ENCODING_FLAG.length(), input.length() - 2);
			byte[] uploadedFileContent = Base64.getDecoder().decode(input);
			
			// create protein
			project = this.modelConverterService.createModelingProject(uploadedFileContent);
		} else {
			// sequence only specified
			project = this.modelConverterService.createModelingProject(ModelConverter.DEFAULT_PROTEIN_NAME, ModelConverter.DEFAULT_PROTEIN_TITLE, input);
		}
		
//		// async approach
//		this.modelPersistenceService.create(project);
//		Executors.newSingleThreadExecutor().submit(() -> this.featureExtractorService.computeFeatures(project.proteins.get(0), FeatureType.values()));
		
		this.featureExtractorService.computeFeatures(project.proteins.get(0), FeatureType.values());
		this.modelPersistenceService.create(project);
		
		// return id to retrieve object later on
		return project._id;
	}
	
	/*
	 * interface to retrieve global application settings and constants
	 */
	
	public APIConstants getSettings(RESTRequest request) {
		return SETTINGS;
	}
	
	/**
	 * handles all PMW-internal logging needs
	 * @param reader
	 */
	@Reference
	void setLogReader(LogReaderService reader) {
		reader.addLogListener(e -> {
			String msg = e.getMessage();
			// suppress framework messages
			if(SUPPRESS_FRAMEWORK_MESSAGES && (msg.startsWith("ServiceEvent") || msg.startsWith("BundleEvent") || msg.startsWith("FrameworkEvent"))) {
				return;
			}
			
			// enforce correct log level
			int level = e.getLevel();
			if(LOG_LEVEL < level) {
				return;
			}
			
			switch (level) {
			case LogService.LOG_DEBUG:
				System.out.println("[DEBUG] " + msg);
				break;
			case LogService.LOG_INFO:
				System.out.println("[INFO] " + msg);
				break;
			case LogService.LOG_WARNING:
				System.out.println("[WARNING] " + msg);
				break;
			case LogService.LOG_ERROR:
				System.err.println("[ERROR] " + msg);
				break;
			}
		});
	}
}