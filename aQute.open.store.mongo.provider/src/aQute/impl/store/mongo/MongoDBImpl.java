package aQute.impl.store.mongo;

import org.osgi.service.component.annotations.Activate;
import org.osgi.service.component.annotations.Component;
import org.osgi.service.component.annotations.Deactivate;
import org.osgi.service.component.annotations.Reference;
import org.osgi.service.log.LogService;

import com.mongodb.DB;
import com.mongodb.Mongo;
import com.mongodb.MongoClient;
import com.mongodb.MongoClientURI;

/**
 * This component is driven by a Managed Service Factory. It opens a Mongo DB,
 * gets a DB object and provides access to the stores. This component implements
 * the aQute.service.store service.
 */
@Component(name = "aQute.open.store.mongo")
public class MongoDBImpl implements aQute.open.store.api.DB {
	Mongo		mongo;
	LogService		log;
	DB	db;
	
	/**
	 * Activate method
	 * @throws Exception 
	 */
	@SuppressWarnings("deprecation")
	@Activate
	void activate() throws Exception {
		//TODO move to config again
		mongo = new MongoClient(new MongoClientURI("mongodb://localhost:27017"));
		db = mongo.getDB("pmw");
	}

	/**
	 * Close the db and unregister the collections
	 */
	@Deactivate
	void deactivate() {
		mongo.close();
	}

	public <T> MongoStoreImpl<T> getStore(Class<T> clazz, String name) throws Exception {
		return new MongoStoreImpl<T>(this, clazz, db.getCollection(name));
	}

	@Override
	public void drop() {
		throw new UnsupportedOperationException("operation 'drop' not implemented");
	}

	@Reference
	public void setLogService(LogService log) {
		this.log = log;
	}
}