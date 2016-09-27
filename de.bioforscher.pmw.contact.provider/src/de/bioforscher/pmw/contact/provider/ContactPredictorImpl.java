package de.bioforscher.pmw.contact.provider;

import org.osgi.service.component.annotations.Component;
import org.osgi.service.component.annotations.Reference;
import org.osgi.service.log.LogService;

import de.bioforscher.pmw.api.ContactPredictor;

@Component(name = "de.bioforscher.pmw.contact")
public class ContactPredictorImpl implements ContactPredictor {
	@Reference
	private LogService logger;
}
