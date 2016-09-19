# Protein Modeling Workbench
A Java-based membrane __Protein Modeling Workbench (PMW)__. The aim is predict the three-dimensional structure of membrane protein sequences and provide an interactive front-end which allows for the visualization of the results and the pieces of information which were utilized in the reconstruction process.

# Demo
Stay tuned for a demo of the current implementation.

# Developing
To be honest - getting this project running is not trivial (even though bndtools & OSGi enroute help a lot).

__IDE requirements__:
* eclipse
* bndtools plugin (https://dl.bintray.com/bndtools/bndtools/3.2.0)
* a bnd/cnf-workspace (this is also possible in IntelliJ IDEA, however compiling the project is significantly slower)
* to get started, you can follow [this tutorial](http://enroute.osgi.org/qs/050-start.html) and ensure that their example application is running

__pmw__:
* clone the repo into the bnd-workspace
* get a [MongoDB Server](https://www.mongodb.com/de) up and running

# Bundles
The application is implemented in a rudimentary service-oriented manner. This allows for modularity and dependency injection.

## API
Specifies the capabilities of all services (which can be found in subsequently listed bundles/services).

## Application
Ties back- and front-end and describes the application on a general level as well as provides the entry point.
#### back-end
exposes access to the database and implemented computations via a `REST` interface
#### front-end
provides an angularJS UI and visualizes the results

## Common
Implementation of `LinearAlgebra` functions utilized for computations by other services.

## Contact
Stub which potentially predict residue-residue, fragment-fragment, motif-motif and helix-helix interactions which are the foundation for the reconstruction cascade.

## Feature Extractor
Can be used to compute features describing residues in more detail (e.g. accessible surface area or secondary structure information).

## Model
The model of the application. Root element is a `Project` (1 job submitted by the user to the application). Projects contain proteins which contain chains which contain residues which contain atoms. Everything all the way potentially has additional information attached to it.

## Model Converter
Can be used to modify the model and work efficiently/conveniently with it. Keywords include:
#### create
new model project instances: either by sequence or PDB file
#### update
PDB representation of protein and residue objects - each stores its relevant information in a PDB-formatted `ATOM` record
#### convert
of 3-letter and 1-letter codes (the model stores 3-letter codes)
#### get
convenient access to residue and atom objects of a protein - backbone atoms have dedicated functions
#### remove
defined atoms or residues from the structure

## Model Persistance
Defines CRUD-operations for the model. Depends on a running MongoDB server.

## Reconstruction
Can reconstruct protein structures in small steps. Typically, predicted contacts are used to compose a distance matrix which describes the distances between all CA atoms of a structure. These information is the foundation of the reconstruction process consisting of multiple steps:
#### place alpha carbons
transform distance map to three-dimensional coordinates by multi-dimensional scaling (MDS) or algorithms such as FT-COMAR
#### place backbone atoms
based on the CA coordinates, the protein backbone can be reconstructed - this is done by Backbone Building from Quadrilaterals (BBQ)
#### place sidechain atoms
based on the backbone coordinates, all atom coordinates of the protein can be computed - the PULCHRA algorithm provides these capabilities
#### refine structure
algorithms such as simulated annealing (SA), can optimize the energy of the system - resulting in reasonably placed atoms
#### evaluate
in the final, the consistency of all atom coordinates is evaluated - if needed previous steps can be executed again, aiming at improving the protein structure quality
