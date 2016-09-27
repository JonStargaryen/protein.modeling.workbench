package de.bioforscher.pmw.model.persistence.provider;

import java.util.List;

import org.osgi.service.component.annotations.Activate;
import org.osgi.service.component.annotations.Component;
import org.osgi.service.component.annotations.Reference;
import org.osgi.service.log.LogService;

import aQute.open.store.api.DB;
import aQute.open.store.api.Store;
import de.bioforscher.pmw.api.ModelPersistence;
import de.bioforscher.pmw.model.DefinedMotif;
import de.bioforscher.pmw.model.Fragment;
import de.bioforscher.pmw.model.Project;

@Component(name = "de.bioforscher.pmw.model.persistence")
public class ModelPersistenceImpl implements ModelPersistence {
	@Reference
	private DB db;
	@Reference
	private LogService logger;

	private Store<Project> projects;
	private Store<Fragment> fragments;
	
	@Activate
	public void activate() {
		this.logger.log(LogService.LOG_INFO, "firing up persistence provider and initing stores");
		try {
			this.projects = this.db.getStore(Project.class, "projects");
			this.fragments = this.db.getStore(Fragment.class, "fragments");
		} catch(Exception e) {
			this.logger.log(LogService.LOG_ERROR, "creating persistence provider failed", e);
		}
	}
	
	@Override
	public void createProject(Project project) throws Exception {
		this.logger.log(LogService.LOG_INFO, "creating persistence entry for project " + project._id);
		this.projects.insert(project);
	}
	
	@Override
	public Project retrieveProject(String projectId) throws Exception {
		this.logger.log(LogService.LOG_INFO, "retrieving persistence entry for project " + projectId);
		Project queryProject = new Project();
		queryProject._id = projectId;
//		return this.projects.find( "select * from Projects where id=%s", projectId).first().get();
		return this.projects.find(queryProject).first().get();
	}

	@Override
	public void updateProject(Project project) throws Exception {
		this.logger.log(LogService.LOG_INFO, "updating persistence entry project " + project._id);
		this.projects.upsert(project);
	}
	
	@Override
	public void deleteProject(String projectId) throws Exception {
		this.logger.log(LogService.LOG_INFO, "deleting persistence entry project " + projectId);
		this.projects.find(projectId).remove();
	}

	@Override
	public List<Fragment> retrieveFragments(DefinedMotif sequenceMotif) throws Exception {
		Fragment queryFragment = new Fragment();
		queryFragment.sequenceMotif = sequenceMotif;
		//TODO check this
		return this.fragments.find(queryFragment).collect();
	}

	@Override
	public void createFragment(Fragment fragment) throws Exception {
		this.logger.log(LogService.LOG_INFO, "creating persistence entry for fragment " + fragment);
		this.fragments.insert(fragment);
	}

	@Override
	public List<Fragment> retrieveFragments() throws Exception {
		this.logger.log(LogService.LOG_INFO, "retrieving all fragments");
		return this.fragments.all().collect();
	}
}
