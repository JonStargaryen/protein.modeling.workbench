'use strict';

(function () {
	var MODULE = angular.module('de.bioforscher.pmw',
			[ 'ngRoute', 'ngResource', 'nvd3', 'ngDropdowns' ]);

	MODULE.constant('design', {
	    defaultColor : '#454b52',
	    darkenedColor : '#a9bcc1'
	});
	
	// the mapping between 3letter and 1letter code of particular amino acids
	MODULE.constant('threeLetterMapping', {
		'ALA' : 'A', 'ARG' : 'R', 'ASN' : 'N', 'ASP' : 'D', 'CYS' : 'C', 'GLN' : 'Q', 'GLU' : 'E', 'GLY' : 'G', 'HIS' : 'H', 'ILE' : 'I', 'LEU' : 'L', 'LYS' : 'K', 'MET' : 'M', 'PHE' : 'F', 'PRO' : 'P', 'SER' : 'S', 'THR' : 'T', 'TRP' : 'W', 'TYR' : 'Y', 'VAL' : 'V'
	});
	
	// mapping of each amino acid to an integer used to color each amino acid the same way regardless of the sequence order
	MODULE.constant('groupMapping', {
		'A' : 1, 'R' : 2, 'N' : 3, 'D' : 4, 'C' : 5, 'Q' : 6, 'E' : 7, 'G' : 8, 'H' : 9, 'I' : 10, 'L' : 11, 'K' : 12, 'M' : 13, 'F' : 14, 'P' : 15, 'S' : 16, 'T' : 17, 'W' : 18, 'Y' : 19, 'V' : 20
	});
	
	MODULE.constant('renderModes',  [
	    "cartoon", "sline", "lines", "trace", "lineTrace", "tube", "spheres", "ballsAndSticks"
	]);
	
	/**
	 * navigation
	 */
	MODULE.config(['$locationProvider', '$routeProvider', function($locationProvider, $routeProvider) {
		// hide index.html# in browser url
		$locationProvider.html5Mode({
			enabled: true
		});
		
		// navigation rules
		$routeProvider.when('/home', {
			templateUrl: '/de.bioforscher.pmw/main/htm/home.htm'
		});
		$routeProvider.when('/about', {
			templateUrl: '/de.bioforscher.pmw/main/htm/about.htm'
		});
		$routeProvider.when('/project', {
			templateUrl: '/de.bioforscher.pmw/main/htm/project.htm'
		});
		$routeProvider.when('/protein', {
			templateUrl: '/de.bioforscher.pmw/main/htm/protein.htm'
		});
		$routeProvider.otherwise('/home');
	}]);
	
	/**
	 * on app start
	 */
	MODULE.run(['$rootScope', '$http', '$location', '$filter', 'renderModes', function ($rootScope, $http, $location, $filter, renderModes) {
		$rootScope.alerts = [];
		
		// load server-managed app constants
		$http.get('/rest/settings/').then(function (d) {
			$rootScope.renderModes = [];
			renderModes.forEach(function(entry) {
				$rootScope.renderModes.push({ text: $filter('splitAtUppercase')(entry), index: $rootScope.renderModes.length, raw: entry });
			});
			// convert Java enums to something easy to handle
			$rootScope.features = handleConstantSet(d.data.features);
			$rootScope.reconstructionLevels = handleConstantSet(d.data.reconstructionLevels);
			$rootScope.secondaryStructures = handleConstantSet(d.data.secondaryStructures);
			$rootScope.topologies = handleConstantSet(d.data.topologies);
			$rootScope.motifTypes = handleConstantSet(d.data.motifTypes);
			$rootScope.interactionTypes = handleConstantSet(d.data.interactionTypes);
		}, function (d) {
			 $rootScope.alerts.push({ type: 'danger', msg: 'loading global constants from server failed with ['+ d.status + '] '+ d.statusText });
		});
		
		// 'parses' 1 constant entry set and attaches it to a container
		handleConstantSet = function(set) {
			var container = [];
			set.forEach(function(entry) {
				// formatting rules are
				container.push({ text: $filter('replaceUnderscores')($filter('lowercase')(entry)), index: container.length, raw: entry });
			});
			return container;
		};
		
		// message/alert container
		$rootScope.closeAlert = function (index) {
			$rootScope.alerts.splice(index, 1);
		};
		$rootScope.page = function () {
			return $location.path();
		}
	}]);
	
	MODULE.controller('ProjectController', ['$scope', '$http', '$location', 'design', function($scope, $http, $location, design) {
		$scope.loading = true;
		$scope.options = {
			renderMembrane : true,
			renderContacts : true,
			renderPDB : false,
			renderMode : $scope.renderModes[0],
			coloringFeature : $scope.features[0]
		};
		
		/* post-construct-like: populate view with data */
		activate = function() {
			$scope.loading = true;
			$http.get('/rest/project/' + $location.search()['id']).then(function (d) {
				// deserialize the returned object to make it usable
				console.log(d.data);
				var project = angular.fromJson(d.data);
				$scope.projectId = project._id;
				$scope.projectName = project.name;
				$scope.projectDate = project.date;
				$scope.proteins = project.proteins;
				$scope.date = project.date;
//				$scope.sequence = project.sequence;
				// join all atom pdbRepresentations to one linked to the protein object
				composePdbRepresentation(project.proteins[0]);
				$scope.protein = project.proteins[0];
				visualizeProteinStructure();
				
				// kick off scrollbar of the sequence/pdb part
				$('.tse-scrollable').TrackpadScrollEmulator();
				
				$scope.loading = false;
			}, function (d) {
				$scope.alerts.push({ type: 'danger', msg: 'failed with [' + d.status + '] ' + d.statusText });
			});
		};
		
		// join all atom pdbRepresentations to one linked to the protein object
		composePdbRepresentation = function (protein) {
			//TODO this could be realized faster using other/sophisticated approaches
			var rep = "";
			protein.chains.forEach(function(chain) {
				chain.residues.forEach(function(residue) {
					residue.atoms.forEach(function(atom) {
						rep += atom.pdbRepresentation + "\n";
					});
				});
			});
			protein.pdbRepresentation = rep;
		};
		
		/* request feature computation */
		this.feature = function(value) {
			// check if this computation is reasonable
//			$scope.protein.availableFeatures.forEach(function(entry) {
//				if(entry.value === value) {
//					return;
//				}
//			});
			calculation('feature', value);
		};
		
		/* request reconstruction step */
		this.reconstruction = function(value) {
			// check if this computation is reasonable
			//TODO	
			
			calculation('reconstruction', value);
		}
		
		/* wraps feature/reconstruction requests and performs them by POSTing to the back-end */
		calculation = function(context, value) {
			if(!$scope.projectId || !context || value < 0) {
				return;
			}
			$http.post('/rest/calculation/', angular.toJson({
				projectId : $scope.projectId,
				context : context,
				value : value
			})).then(function (d) {
				activate();
				$scope.alerts.push({ type: 'success', msg: 'computed ' + d.data });
			}, function (d) {
				$scope.alerts.push({ type: 'danger', msg: 'failed with ['+ d.status + '] '+ d.statusText });
			});
		}
		
		/* PV related functions */
		visualizeProteinStructure = function (protein) {
			var options = {
				// div to be selected is not visible at that time
				width: 400,
				height: 400,
				antialias: true,
				quality : 'high',
				background:'#313a41'
			};

			if($scope.protein.pdbRepresentation === "") {
				// no atom coordinates - this is a modeling project
				return;
			}
			
			// ensure container is empty
			document.getElementById('protein-visualizer').innerHTML = '';
			
			$scope.viewer = pv.Viewer(document.getElementById('protein-visualizer'), options);
			$scope.structure = io.pdb($scope.protein.pdbRepresentation);
			mol.assignHelixSheet($scope.structure);
			$scope.viewer.options('fog', false);
			
			renderStructure();
		    renderMembrane();
		    renderContacts();
		    
			$scope.viewer.centerOn($scope.structure);
			$scope.viewer.autoZoom();
		};
		
		/* dedicated functions to render information (if appropriate and available) */
		renderStructure = function() {
			console.log('[PV] requested to draw using ' + $scope.options.renderMode.text);
			$scope.viewer.renderAs('protein', $scope.structure, $scope.options.renderMode.raw, { color: pv.color.uniform(design.darkenedColor) });
		}
		
		renderContacts = function() {
			if($scope.protein.contacts.length === 0 || !$scope.options.renderContacts) {
				return;
			}
			
		    // visualize contacts in protein structure
		    graph.links.forEach(function(link) {
		    	// sequential residues
		    	if(link.value > 1) {
		    		return;
		    	}
		    	var name1 = link.source.name.split("-")[0] + '.' + link.source.name.split("-")[2] + '.CA';
		    	var name2 = link.target.name.split("-")[0] + '.' + link.target.name.split("-")[2] + '.CA';
		    	var res1 = $scope.structure.atom(name1);
		    	var res2 = $scope.structure.atom(name2);
		    	var g = $scope.viewer.customMesh('contact');
		        g.addTube(res1.pos(), res2.pos(), 0.2, { cap : false, color : 'red' });
		    });
		}
		
		renderMembrane = function() {
		    if(!$scope.protein.membrane || !$scope.options.renderMembrane) {
		    	return;
		    }
		    
	    	$scope.protein.membrane.membraneMolecules.forEach(function(entry) {
	    		$scope.viewer.customMesh('membrane').addSphere([entry[0], entry[1], entry[2]], 0.75, { color : [1, 1, 1, 0.5] });
	    	});
		}
		
		/* monitor option changes and if appropriate manipulate the PV's visualization */
		$scope.$watch('options.renderMode', function(newVal, oldVal){
			if(newVal === oldVal)
				return;
//			console.log('render mode changed from ' + oldVal + ' to ' + newVal);
		    $scope.viewer.rm('protein');
		    renderStructure();
		});
		
		$scope.$watch('options.renderMembrane', function(newVal, oldVal){
			if(newVal === oldVal)
				return;
//			console.log('render membrane changed from ' + oldVal + ' to ' + newVal);
		    newVal ? $scope.viewer.show('membrane') : $scope.viewer.hide('membrane');
		    $scope.viewer.requestRedraw();
		});
		
		$scope.$watch('options.renderContacts', function(newVal, oldVal){
			if(newVal === oldVal)
				return;
//			console.log('render contacts changed from ' + oldVal + ' to ' + newVal);
		    newVal ? $scope.viewer.show('contact') : $scope.viewer.hide('contact');
		    $scope.viewer.requestRedraw();
		});
		
		// init post construct routine
		activate();
	}]);
	
//	MODULE.controller('ProjectController', ['$scope', '$http', '$location', 'design', 'groupMapping', function ($scope, $http, $location, design, groupMapping) {
//		$scope.project = null;
//		$scope.options = {
//			renderMembrane : true,
//			renderContacts : true,
//			featureIndex : 0,
//			renderModeIndex : 0
//		};
//		$scope.features;
//		$scope.reconstructionLevels;
//		$scope.renderModes;
//
//		var viewer = null;
//		var structure = null;
//		var graph = { "nodes":[], "links":[] };
//		var computationRunning = false;
//		
//		loadProjectByID = function (id) {
//			if (!id)
//				return;
//			
//			$http.get('/rest/project/' + id).then(
//				function (d) {
//					// deserialize the returned object (JSONObject.toString()) and make it usable
//					var project = angular.fromJson(d.data);
//					$scope.project = project;
//					var protein = $scope.project.proteins[0];
//					
//					protein.availableFeatureNames = [];
//					protein.availableFeatures.forEach(function(entry) {
//						protein.availableFeatureNames.push(entry.name);
//					});
//					// check whether there is information to visualize
//					if(protein.residueResidueInteractions.length > 0) {
//						buildGraph(protein);
//						visualizeProteinGraph(protein);
//						visualizeContactMatrix(protein);
//					}
//					// same for 3d coordinates
//					if(protein.reconstructionLevel.value > 0) {
//						visualizeProteinStructure(protein);
//					}
//				}, function (d) {
//					$scope.alerts.push({ type: 'danger', msg: 'failed with [' + d.status + '] ' + d.statusText });
//				});
//		};
//		
//		this.requestFeature = function(value) {
//			var alreadyPresent = false;
//			$scope.project.proteins[0].availableFeatures.forEach(function(entry) {
//				if(entry.value === value) {
//					alreadyPresent = true;
//				}
//			});
//			if(alreadyPresent) {
//				return;
//			}
//			requestCalculation('feature', value);
//		};
//		
//		this.requestReconstruction = function(value) {
//			// don't do anything when level is higher than the one requested
////			if(value <= $scope.project.proteins[0].reconstructionLevel) {
////				return;
////			}
//			requestCalculation('reconstruction', value);
//		}
//		
//		requestCalculation = function(context, value) {		
//			if($scope.project && context && value) {
//				$http.get('/rest/calculation/' + $scope.project.projectId + '/' + context + '/' + value).then( //$http.post('/rest/calculation/', angular.toJson(request)).then(
//					function (d) {
//						loadProjectByID($scope.project.projectId);
//						$scope.alerts.push({ type: 'success', msg: 'computed ' + d.data });
//					}, function (d) {
//						$scope.alerts.push({ type: 'danger', msg: 'failed with ['+ d.status + '] '+ d.statusText });
//					});
//			}
//		}
//		
//		buildGraph = function(protein) {	
//			//TODO: for performance: only update graph/move/simulate when it is actually visable
//			// add sequentially connected amino acids as well as their links
//			var residueCount = 0;
//			protein.chains.forEach(function(chain) {
//				var hasPreviousResidueInChain = false;
//				chain.residues.forEach(function(residue) {
//					graph.nodes.push({"name":chain.chainId + "-" + residue.aminoAcid + "-" + residue.residueNumber,"group":mapToGroup(residue.aminoAcid)});
//					// connect sequential neighbors
//					if(hasPreviousResidueInChain) {
//						graph.links.push({"source":residueCount, "target":residueCount-1, "value":3});
//					}
//					residueCount++;
//					hasPreviousResidueInChain = true;
//				});
//			});
//			
//			// add spatially connected residue links
//			protein.residueResidueInteractions.forEach(function(rri) {
//				graph.links.push({"source":rri.interactingResidues[0].residueId,"target":rri.interactingResidues[1].residueId,"value":1});
//			});
//		}
//		
//		mapToGroup = function(aminoAcid) {
//			return groupMapping[aminoAcid];
//		};
//		
//		visualizeProteinStructure = function (protein) {
//			var options = {
//				// div to be selected is not visible at that time
//				width: d3.select(".visualizer-tab:not(.ng-hide)")[0][0].clientWidth,
////				width: d3.select("#sequence-visualizer")[0][0].clientWidth,
//				height: 500,
//				antialias: true,
//				quality : 'high',
//				background:'#eee'
//			};
//
//			// ensure container is empty
//			document.getElementById('protein-visualizer').innerHTML = '';
//			
//			//TODO: this seems quite hacky
//			viewer = pv.Viewer(document.getElementById('protein-visualizer'), options);
//			structure = io.pdb(protein.pdbRepresentation);
//			mol.assignHelixSheet(structure);
//			viewer.options('fog', false);
//			
//			renderStructure();
//		    renderMembrane();
//		    renderContacts();
//		    
//			viewer.centerOn(structure);
//		    viewer.autoZoom();
//		};
//		
//		$scope.$watch('options.renderMembrane', function(newVal, oldVal){
//			if(newVal === oldVal)
//				return;
//		    newVal ? viewer.show('membrane') : viewer.hide('membrane');
//		    viewer.requestRedraw();
//		});
//		
//		$scope.$watch('options.renderContacts', function(newVal, oldVal){
//			if(newVal === oldVal)
//				return;
//		    newVal ? viewer.show('contact') : viewer.hide('contact');
//		    viewer.requestRedraw();
//		});
//		
//		$scope.$watch('options.renderModeIndex', function(newVal, oldVal){
//			if(newVal === oldVal)
//				return;
//		    viewer.rm('protein');
//		    renderStructure();
//		});
//		
//		renderStructure = function() {
//			viewer.renderAs('protein', structure, $scope.renderModes[$scope.options.renderModeIndex], { color: pv.color.uniform(design.darkenedColor) });
//		}
//		
//		renderContacts = function() {
//			if($scope.project.proteins[0].residueResidueInteractions.length === 0 || !$scope.options.renderContacts)
//				return;
//			
//		    // visualize contacts in protein structure
//		    graph.links.forEach(function(link) {
//		    	// sequential residues
//		    	if(link.value > 1)
//		    		return;
//		    	var name1 = link.source.name.split("-")[0] + '.' + link.source.name.split("-")[2] + '.CA';
//		    	var name2 = link.target.name.split("-")[0] + '.' + link.target.name.split("-")[2] + '.CA';
//		    	var res1 = structure.atom(name1);
//		    	var res2 = structure.atom(name2);
//		    	var g = viewer.customMesh('contact');
//		        g.addTube(res1.pos(), res2.pos(), 0.2, { cap : false, color : 'red' });
//		    });
//		}
//		
//		renderMembrane = function() {
//		    if(!$scope.project.proteins[0].membrane || !$scope.options.renderMembrane)
//		    	return;
//		    
//	    	$scope.project.proteins[0].membrane.membraneMolecules.forEach(function(entry) {
//	    		viewer.customMesh('membrane').addSphere([entry.x, entry.y, entry.z], 0.75, { color : [1.0, 0.75, 0.5, 0.35] });
//	    	});
//		}
//		
//		visualizeContactMatrix = function(protein) {
////			var width = d3.select("#sequence-visualizer")[0][0].clientWidth,
////				height = d3.select("#sequence-visualizer")[0][0].clientWidth;
//			var width = protein.size,
//				height = protein.size;
//
//			var x = d3.scale.ordinal().rangeBands([0, width]),
//		    	z = d3.scale.linear().domain([0, 4]).clamp(true),
//		    	color = d3.scale.category20();
//
//			var svg = d3.select("#contact-visualizer").append("svg")
//		    	.attr("width", width)
//		    	.attr("height", height)
//		    	.append("g");
//		
//			var matrix = [],
//		    	nodes = graph.nodes,
//		    	n = nodes.length;
//
//			// Compute index per node.
//			nodes.forEach(function(node, i) {
//				node.index = i;
//				node.count = 0;
//				matrix[i] = d3.range(n).map(function(j) { return {x: j, y: i, z: 0}; });
//			});
//
//			
//			// Convert links to matrix; count character occurrences.
//			graph.links.forEach(function(link) {
//				matrix[link.source.index][link.target.index].z += link.value;
//				matrix[link.target.index][link.source.index].z += link.value;
//				matrix[link.source.index][link.source.index].z += link.value;
//				matrix[link.target.index][link.target.index].z += link.value;
//				nodes[link.source.index].count += link.value;
//				nodes[link.target.index].count += link.value;
//			});
//
//			// Precompute the orders.
//			var orders = {
//				name: d3.range(n).sort(function(a, b) { return d3.ascending(nodes[a].index, nodes[b].index); }),
////				count: d3.range(n).sort(function(a, b) { return nodes[b].count - nodes[a].count; }),
////				group: d3.range(n).sort(function(a, b) { return nodes[b].group - nodes[a].group; })
//			};
//			
//			// The default sort order.
//			x.domain(orders.name);
//			
//			svg.append("rect")
//				.attr("class", "background")
//				.attr("width", width)
//				.attr("height", height);
//
//			var row = svg.selectAll(".row")
//				.data(matrix)
//				.enter().append("g")
//				.attr("class", "row")
//				.attr("transform", function(d, i) { return "translate(0," + x(i) + ")"; })
//				.each(row);
//			
//			row.append("line")
//				.attr("x2", width);
//			
//			var column = svg.selectAll(".column")
//				.data(matrix)
//				.enter().append("g")
//				.attr("class", "column")
//				.attr("transform", function(d, i) { return "translate(" + x(i) + ")rotate(-90)"; });
//				
//			column.append("line")
//				.attr("x1", -width);
//			
//			function row(row) {
//				var cell = d3.select(this).selectAll(".cell")
//					.data(row.filter(function(d) { return d.z; }))
//					.enter().append("rect")
//					.attr("class", "cell")
//					.attr("x", function(d) { return x(d.x); })
//					.attr("width", x.rangeBand())
//					.attr("height", x.rangeBand())
//					.style("fill-opacity", function(d) { return z(d.z); })
////					.style("fill", function(d) { return nodes[d.x].group == nodes[d.y].group ? c(nodes[d.x].group) : null; });
//					.style("fill", function(d) { return color(d.group); });
//			}
//		};
//		
//		visualizeProteinGraph = function (protein) {
////		var width = d3.select("#sequence-visualizer")[0][0].clientWidth, height = 500, radius = 2.5;
//			var width = d3.select(".visualizer-tab:not(.ng-hide)")[0][0].clientWidth, height = 500, radius = 2.5;
//			
//			var color = d3.scale.category20();
//			var force = d3.layout.force()
//		    	.gravity(.05)
//		    	.charge(-5)
//		    	.linkDistance(3.5)
//		    	.size([width, height]);
//			
//			var svg = d3.select("#graph-visualizer").append("svg")
//				.attr("width", width)
//				.attr("height", height);
//				
//			force.nodes(graph.nodes)
//				.links(graph.links)
//				.start();
//			
//			var link = svg.selectAll(".link")
//				.data(graph.links)
//			    .enter().append("line")
//				.attr("class", "link")
//				.style("stroke-width", function(d) { return Math.sqrt(d.value); });
//
//			var node = svg.selectAll(".node")
//		    	.data(graph.nodes)
//		    	.enter().append("circle")
//		    	.attr("r", radius - .75)
//		    	.style("fill", function(d) { return color(d.group); })
//		    	.style("stroke", function(d) { return d3.rgb(color(d.group)).darker(); })
//		    	.call(force.drag);
//	
//
//			node.append("title")
//				.text(function(d) { return d.name; });
//
//			force
//				.nodes(graph.nodes)
//				.links(graph.links)
//				.on("tick", tick)
//				.start();
//			  
//			function tick() {
//				node.attr("cx", function(d) { return d.x = Math.max(radius, Math.min(width - radius, d.x)); })
//					.attr("cy", function(d) { return d.y = Math.max(radius, Math.min(height - radius, d.y)); });
//
//				link.attr("x1", function(d) { return d.source.x; })
//					.attr("y1", function(d) { return d.source.y; })
//			        .attr("x2", function(d) { return d.target.x; })
//			        .attr("y2", function(d) { return d.target.y; });
//			}
//		}
//		
//		// load global settings from server - TODO: move this somewhere more appropiate
////		if(!$scope.features) {
//		$http.get('/rest/settings/').then(function (d) {
//			var settings = angular.fromJson(d.data);
//
//			// assign global 'constants' - depending on the back-end these values may change with version updates - TODO: there should be a better way to propagate them to the front-end
//			$scope.features = settings.supportedFeatures;
//			console.log($scope.features);
//			$scope.reconstructionLevels = settings.supportedReconstructionLevels;
//			console.log($scope.reconstructionLevels);
//			$scope.renderModes = settings.supportedRenderModes;
//			console.log($scope.renderModes);
//		}, function (d) {
//			$scope.alerts.push({ type: 'danger', msg: 'failed with ['+ d.status + '] '+ d.statusText });
//		});
////		}
//		
//		// as post construct: load the project by the search parameter
////		loadProjectByID($location.search()['id']);
//	}]);
	
	/**
	 * handles multiple tabs in the UI
	 */
	MODULE.controller('TabController', function () {
		this.tab = 1;
		this.selectTab = function (setTab){
			this.tab = setTab;
		};
		this.isSelected = function (checkTab){
			return this.tab === checkTab;
		};
	});
	
	/**
	 * handles the upload of user input on the initial project creation page
	 */
	MODULE.controller('SubmitController', ['$scope', '$http', '$location', function ($scope, $http, $location) {
		$scope.sequence = 'QITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGE';
		
		$scope.createProtein = function () {
			if($scope.sequence) {
				$http.post('/rest/project/', angular.toJson($scope.sequence)).then(
				function (d) {
					$location.path("/project").search({ id: d.data });
				}, function (d) {
					$scope.alerts.push({ type: 'danger', msg: 'failed with ['+ d.status + '] '+ d.statusText });
				});
			}
		};
	}]);
	
	/**
	 * filters for UI output of back-end enums
	 */
	// format PV render modes to something nice
	MODULE.filter('splitAtUppercase', function() {
		return function(input) {
			if(input) {
				return input
				// insert a space before all caps
			    .replace(/([A-Z])/g, ' $1')
			    .toLowerCase();
			}
		}
	});
	
	// convert Java enum underscores to whitespaces
	MODULE.filter('replaceUnderscores',function() {
	    return function(input) {
	        if (input) {
	            return input.replace(/_/g, ' ');    
	        }
	    }
	});
	
	// convert 3letter to 1letter code
	MODULE.filter('oneLetterCode', ['threeLetterMapping', function(threeLetterMapping) {
		return function(input) {
			if(input) {
				var f = threeLetterMapping[input];
				return f ? f : 'X';
			}
		}
	}]);
	
	MODULE.filter('helixRenaming', function() {
		return function(input) {
			if(input) {
				return input.replace("pi", "\u03C0-").replace("alpha ", "\u03B1-").replace("three10", "3-10-");
			}
		}
	});
	
	/**
	 * uploads files to the back-end as base64-encoded json query
	 */
	MODULE.directive('pdbUpload', ['$location', 'httpPostFactory', function ($location, httpPostFactory) {
	    return {
	        restrict: 'A',
	        scope: true,
	        link: function (scope, element, attr) {
	            element.bind('change', function () {
	            	var fr = new FileReader();
	            	fr.onload = function(event) {
	            		var data = {'file':event.target.result};
		                httpPostFactory('/rest/project/', data, function (callback) {
		                	$location.path("/protein").search({ id: callback });
		                });
	            	};
	            	var file = fr.readAsDataURL(element[0].files[0]);
	            });

	        }
	    };
	}]);
	
	MODULE.factory('httpPostFactory', function ($http) {
	    return function (file, data, callback) {
	        $http({
	            url: file,
	            method: "POST",
	            data: data,
	            headers: {'Content-Type': undefined}
	        }).success(function (response) {
	            callback(response);
	        });
	    };
	});
})();