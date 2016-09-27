package de.bioforscher.pmw.model;

import java.util.ArrayList;
import java.util.List;

import org.osgi.dto.DTO;

/**
 * The root of the model. For sake of the persistence API, this was moved from EMF to DTOs. It feels kinda wrong however that all fields are directly accessible - however, is actually quite reasonable since it ensures that the representation of the data is the same in each part of the application (in the database, in the back-end as well as in its JSON-form).
 * @author S
 *
 */
public class Project extends DTO {
	public String _id;
	public long date;
	public List<Protein> proteins;
	// removed sequence - now the first protein is the general information container :x - maybe make that a dedicated field - let's see
//	public String sequence;
	public String name;
	
	public Project() {
		this.proteins = new ArrayList<>();
	}
}
