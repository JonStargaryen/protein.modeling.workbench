package de.bioforscher.pmw.api;

import java.util.List;
import java.util.NoSuchElementException;

import de.bioforscher.pmw.model.DefinedMotif;
import de.bioforscher.pmw.model.Fragment;
import de.bioforscher.pmw.model.Project;

/**
 * The abstraction of the persistence layer - provides the classic CRUD-operations in the context of {@link Project}.
 * 
 * @author S
 *
 */
public interface ModelPersistence {
	/**
	 * persists the handed over {@link ModelingJob}<br />
	 * this is done by ensuring that the internally managed list remains unique (that means, existing references for this projectId are removed beforehand) and, subsequently, the new object is added to the list
	 * @param project
	 */
	void createProject(Project project) throws Exception;
	
	/**
	 * tries to retrieve a persisted {@link ModelingProject} by its projectId attribute - when none can be found, a {@link NoSuchElementException} is thrown
	 * @param projectId the id of the object to be retrieved
	 * @return a reference to the retrieved object
	 */
	Project retrieveProject(String projectId) throws Exception;
	
	/**
	 * updates an existing entry
	 * @param project
	 */
	void updateProject(Project project) throws Exception;
	
	/**
	 * removes a object identified by its projectId
	 * @param projectId the id of the object to be removed - when none can be found, a {@link NoSuchElementException} is thrown
	 */
	void deleteProject(String projectId) throws Exception;
	
	List<Fragment> retrieveFragments() throws Exception;
	
	/**
	 * TODO design: really what is the query? only a motif? the whole sequence? also topology information?
	 * @param sequenceMotif
	 * @return
	 */
	List<Fragment> retrieveFragments(DefinedMotif sequenceMotif) throws Exception;
	
	/**
	 * TODO design: is a fragment unique for one bin or do we store only one representative fragment - do we merge them outside or here when collisions occur
	 * @param fragment the fragment data to add
	 */
	void createFragment(Fragment fragment) throws Exception;
}
