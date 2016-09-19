package de.bioforscher.pmw.api;

import java.util.NoSuchElementException;

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
	void create(Project project) throws Exception;
	
	/**
	 * tries to retrieve a persisted {@link ModelingProject} by its projectId attribute - when none can be found, a {@link NoSuchElementException} is thrown
	 * @param projectId the id of the object to be retrieved
	 * @return a reference to the retrieved object
	 */
	Project retrieve(String projectId) throws Exception;
	
	/**
	 * updates an existing entry
	 * @param project
	 */
	void update(Project project) throws Exception;
	
	/**
	 * removes a object identified by its projectId
	 * @param projectId the id of the object to be removed - when none can be found, a {@link NoSuchElementException} is thrown
	 */
	void delete(String projectId) throws Exception;
}
