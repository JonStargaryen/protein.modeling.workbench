package de.bioforscher.pmw.model.persistence.provider;

import org.osgi.service.component.annotations.Activate;
import org.osgi.service.component.annotations.Component;
import org.osgi.service.component.annotations.Reference;
import org.slf4j.Logger;

import aQute.open.store.api.DB;
import aQute.open.store.api.Store;
import de.bioforscher.pmw.api.ModelPersistence;
import de.bioforscher.pmw.model.Project;

@Component(name = "de.bioforscher.pmw.model.persistence")
public class ModelPersistenceImpl implements ModelPersistence {
	@Reference
	private Logger logger;
	@Reference
	private DB db;
	private Store<Project> projects;
	
	@Activate
	public void activate() throws Exception {
		this.logger.info("firing up mongodb persistence service");
		this.projects = this.db.getStore(Project.class, "projects");
	}
	
	@Override
	public void create(Project project) throws Exception {
		this.logger.info("creating persistence entry for " + project._id);
		this.projects.insert(project);
	}
	
	@Override
	public Project retrieve(String projectId) throws Exception {
		this.logger.info("retrieving persistence entry for " + projectId);
		Project project = new Project();
		project._id = projectId;
//		return this.projects.find( "select * from Projects where id=%s", projectId).first().get();
		return this.projects.find(project).first().get();
	}

	@Override
	public void update(Project project) throws Exception {
		this.logger.info("updating persistence entry " + project._id);
		this.projects.upsert(project);
	}
	
	@Override
	public void delete(String projectId) throws Exception {
		this.logger.info("deleting persistence entry " + projectId);
		this.projects.find(projectId).remove();
	}
}
