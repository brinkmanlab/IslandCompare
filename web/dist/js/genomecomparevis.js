//Dependencies: jQuery, D3.js, d3.phylogram.js, newick.js, linearplot.css (for styling)

function MultiVis(targetNode){
    var self = this;
    const TOPPADDING = 24;
    const SEQUENCEHEIGHT = 20;
    const CONTAINERWIDTH = null;
    const TREECONTAINERWIDTH = 195;
    const TREETOPPADDING = -16+TOPPADDING;
    const LEFTPADDING = 185+TREECONTAINERWIDTH;
    const TEXTTOPPADDING = 19;
    const TEXTPADDING = 30;
    const GIFILTERFACTOR = 4000;
    const GENEFILTERFACTOR =400000;
    const SEQUENCEWIDTH=10;
    const GISIZE = SEQUENCEWIDTH;
    const GENESIZE = GISIZE/2 +1;
    const RIGHTPADDING = 40;
    const MAXIMUMINTERVALSIZE = 160;
    const XAXISHEIGHT = 100;
    const GIDEFAULTCOLOR = "#ff9933";

    this.container = d3.select(targetNode);
    this.backbone = new Backbone();
    this.sequences = this.backbone.getSequences();
    this.scale = null;
    this.sequenceOrder = null;
    this.isPrinterColors = false;
    this.verticalScrollVal = 0;
    this.newickData = null;
    this.newickRoot = null;
    this.trueBranchLengths = false;

    this.toggleTrueBranchLengths = function(){
        this.trueBranchLengths = !this.trueBranchLengths;
    };

    // Updates the vertical scale of the graph depending on numberSequences(int)
    // If the number of sequences are over 30 than no expanding of graph will occur,
    // Otherwise graph will expand to fill the space available to it
    this.updateVerticalScrollVal = function(numberSequences){
        var height = window.innerHeight
            || document.documentElement.clientHeight
            || document.body.clientHeight;

        if (numberSequences >= 30){
            this.verticalScrollVal = 0;
        }
        else{
            this.verticalScrollVal = ((height-TOPPADDING-XAXISHEIGHT)/numberSequences)-SEQUENCEHEIGHT;
            if (this.verticalScrollVal>MAXIMUMINTERVALSIZE){
                this.verticalScrollVal = MAXIMUMINTERVALSIZE;
            }
        }
    };

    // Returns the scale modified height of the graph
    this.getSequenceModHeight = function(){
        return SEQUENCEHEIGHT + this.verticalScrollVal;
    };

    // Returns the scale modified padding of the tree to get the alignment correct
    this.getTreeModPadding = function(){
        return TREETOPPADDING - (this.verticalScrollVal/2);
    };

    // TODO improve this (Not used in current implementation of alignervis)
    this.setSequenceOrderFromNames = function(arrayOrder){
        var order = [];
        for (var nameIndex=0;nameIndex<arrayOrder.length;nameIndex++){
            for (var seqIndex=0;seqIndex<this.sequences.length;seqIndex++){
                if (arrayOrder[nameIndex]===this.sequences[seqIndex].sequenceName){
                    order.push(this.sequences[seqIndex].sequenceId);
                }
            }
        }
        this.sequenceOrder=order;
    };

    // Returns the sequence order, if it is not initialized will create a default order from 0 -> sequence length -1
    this.getSequenceOrder = function(){
        if (this.sequenceOrder == null){
            this.sequenceOrder = [];
            for (var index = 0; index<this.sequences.length; index++){
                this.sequenceOrder.push(index);
            }
        }
        return this.sequenceOrder;
    };

    // Reorders the sequences (array) based on the current sequence order value
    this.reorderSequences = function(){
        var newOrder = [];
        var seqOrder = this.getSequenceOrder();
        var currentSeqs = this.sequences;
        for (var index=0;index<seqOrder.length;index++){
            newOrder.push(currentSeqs[seqOrder[index]]);
        }
        this.sequences = newOrder;
    };

    this.getGIFilterValue = function(){
        // Converted to show all gis as long as they contain more nucleotides than the gifilterfactor
        /*
         var windowSize = self.scale.domain()[1]-self.scale.domain()[0];
         return windowSize/self.getLargestSequenceSize()*GIFILTERFACTOR;
         */
        return GIFILTERFACTOR
    };

    this.getGeneFilterValue = function(){
        // Converted to a show all genes or none algorithm so genefiltervalue stays constant now
        /*
         var windowSize = self.scale.domain()[1]-self.scale.domain()[0];
         return windowSize/self.getLargestSequenceSize()*GENEFILTERFACTOR;
         */
        return GENEFILTERFACTOR;
    };

    // This sets the horizontal ranges of the sequences on the graph.
    // Where start is the base with the smallest position in the graph,
    // and end is the base with the largest position in the graph
    this.setScale = function(start,end){
        this.scale = d3.scale.linear()
            .domain([start,end])
            .range([0,this.visualizationWidth()])
            .clamp(false);
    };

    this.getSequence = function(index){
        return self.sequences[index];
    };

    this.getLargestSequenceSize = function(){
        var largestSize=0;
        for (var i = 0; i<this.sequences.length; i++){
            if (this.sequences[i].getSequenceSize()>largestSize){
                largestSize = this.sequences[i].getSequenceSize();
            }
        }
        return largestSize;
    };

    this.containerWidth = function() {
        if (CONTAINERWIDTH != null){
            return CONTAINERWIDTH;
        }
        else {
            return this.container.node().getBoundingClientRect().width;
        }
    };

    this.visualizationWidth = function(){
        if (CONTAINERWIDTH != null){
            return CONTAINERWIDTH-LEFTPADDING-RIGHTPADDING;
        }
        else{
            return this.container.node().getBoundingClientRect().width-LEFTPADDING;
        }
    };

    this.containerHeight = function() {
        //return this.sequences.length*this.getSequenceModHeight()-100;// The -100 fixes padding issues on islandviewer site, fix this later if required;
        return this.sequences.length*(this.getSequenceModHeight());
    };

    this.updateSequenceVisualization= function(sequenceIndex, newstart, newend){
        try {
            self.getSequence(sequenceIndex).updateScale(newstart, newend, this.visualizationWidth());
            self.transition();
        }catch(e){
            //Chart may not have been initialized yet (call parseandrender?)
        }
    };

    //Readjusts the graph for updated sequence domains, TODO (improve later, currently just re-renders graph)
    this.transition = function(){
        this.container.select("svg").remove();
        this.render();
    };

    //Resets the range of the graph and rerenders it
    this.resetAndRenderRange = function(){
        this.setScale(0,this.getLargestSequenceSize());
        this.transition();
    };

    //Resets the pointer to the root of the tree and rerenders it
    this.resetAndRenderGraph = function(){
        this.newickRoot = this.newickData;
        this.sequenceOrder = null;
        this.transition();
    };

    //Renders the graph
    this.render = function (){
        this.container.select("svg").remove();
        if (self.scale == null){
            self.setScale(0,this.getLargestSequenceSize());
        }
        // Gets the sequenceOrder of the graph
        var seqOrder = self.getSequenceOrder();

        // Modifes the spacing between sequences depending on the number of sequences to show
        self.updateVerticalScrollVal(seqOrder.length);

        //Add the SVG (Make sequence height 2 sequence higher than container height to fit svg TODO Refactor container height)
        var svg = this.container.append("svg")
            .attr("width",this.containerWidth())
            .attr("height",(this.containerHeight()+this.getSequenceModHeight()*2)+TOPPADDING);

        //Add the visualization container
        var visContainer = svg.append("g")
            .attr("class","visContainer");

        //Renders the tree if a root is set
        if (self.newickRoot != null) {
            //Add the tree visualization container
            var treeContainer = svg.append("g")
                .attr("class", "treeContainer")
                .attr("width", TREECONTAINERWIDTH)
                .attr("height", this.containerHeight())
                .attr("transform", "translate(" + 0 + "," + (this.getTreeModPadding()) + ")");

            // Add the tree and listeners to the visualization
            d3.phylogram.build('.treeContainer', self.newickRoot, {
                width: TREECONTAINERWIDTH,
                height: seqOrder.length * (this.getSequenceModHeight()),
                skipLabels: true, //removes the labels from the tree (can be shown to ensure mapping correctly)
                skipBranchLengthScaling: !this.trueBranchLengths,
                nodeCallback: function (d) {
                    // Sets the root of the tree to the clicked node
                    self.newickRoot = d;
                    //Set the order of the sequences (as linear plot and tree should match)
                    //Do a post-order traversal on tree looking for leaves.
                    var tree = [];
                    traverseTreeForOrder = function (node) {
                        if (node['children'] != null) {
                            for (var nodeIndex = 0; nodeIndex < node['children'].length; nodeIndex++) {
                                output = traverseTreeForOrder(node['children'][nodeIndex]);
                                if (output != null) {
                                    tree.push(output);
                                }
                            }
                        }
                        else {
                            var memory = (node['name'].replace(/'/g, "").split(/[\\/]/).pop());
                            memory = memory.split(".")[0];
                            return self.backbone.getSequenceIdFromName(memory);
                        }
                        return null;
                    };
                    traverseTreeForOrder(self.newickRoot);
                    //Set the sequenceOrder of the tree to the sequences obtained from traversal of the new treeRoot
                    self.sequenceOrder = tree;
                    //Re-renders the graph
                    self.transition();
                }
            });
        }

        //Holds the linear plot visualization except the scale to prevent clipping/overflow problems
        var sequenceHolder = visContainer.append("svg")
            .attr("width",this.visualizationWidth())
            .append("g")
            .attr("transform","translate("+ 0 +","+GISIZE+")");

        //Draw Homologous Region Lines
        var lines = [];

        for (var i=0; i<this.sequences.length-1; i++){
            var seqlines = sequenceHolder.append("g")
                .attr("class","all-homolous-regions");
            //Find a way to clean this line up
            var homologousRegions = (this.backbone.retrieveHomologousRegions(seqOrder[i],seqOrder[i+1]));

            //Homologous Region polygons (shaded regions)
            for (var j=0;j<homologousRegions.length;j++){
                var homolousRegion = seqlines.append("g")
                    .attr("class", "homologous-region")
                    .attr("transform","translate("+ 0 +","
                        +SEQUENCEWIDTH/2+")");

                //Build Shaded Polygon For Homologous Region
                var points = self.scale(this.sequences[i].scale(homologousRegions[j].start1))+","+(this.getSequenceModHeight()*i+SEQUENCEWIDTH/2)+" ";
                points += self.scale(this.sequences[i].scale(homologousRegions[j].end1))+","+(this.getSequenceModHeight()*i+SEQUENCEWIDTH/2)+" ";
                points += self.scale(this.sequences[i+1].scale(homologousRegions[j].end2))+","+(this.getSequenceModHeight()*(i+1)-SEQUENCEWIDTH/2)+" ";
                points += self.scale(this.sequences[i+1].scale(homologousRegions[j].start2))+","+(this.getSequenceModHeight()*(i+1)-SEQUENCEWIDTH/2)+" ";

                homolousRegion.append("polygon")
                    .attr("points",points)
                    .attr("stroke-width",1)
                    .append("title")
                    .text("["+homologousRegions[j].start1+","+homologousRegions[j].end1+"],"+
                        "["+homologousRegions[j].start2+","+homologousRegions[j].end2+"]");
            }
        }
        //Create the sequence list from the ordered sequence list
        var sequencedata = [];
        for (var index=0;index<seqOrder.length;index++){
            sequencedata.push(this.sequences[seqOrder[index]]);
        }

        //Create the sequences container on the svg
        var seq = sequenceHolder.selectAll("sequencesAxis")
            .data(sequencedata)
            .enter()
            .append("g")
            .attr("class", "sequences");

        //Add the sequences to the SVG
        seq.append("rect")
            .attr("x",0)
            .attr("y",function (d, i){
                return (i*self.getSequenceModHeight())+"px";
            })
            .attr("height", SEQUENCEWIDTH)
            .attr("rx",4)
            .attr("ry",4)
            .attr("width", function (d){
                return self.scale(d.scale(d.getSequenceSize()));
            });

        //Add the gis to the SVG
        var giFiltervalue = self.getGIFilterValue();
        var gis = seq.each(function(d, i){
            var genomicIslandcontainer = seq.append("g")
                .attr("class","genomicIslands")
                .attr("transform","translate("+ 0 +","
                    +0+")");
            for (var prediction_name in d.gi){
                var giClassContainer = genomicIslandcontainer.append("g")
                    .attr("class", prediction_name);

                for (var giIndex=0;giIndex< d.gi[prediction_name].length;giIndex++){
                    var startPosition = parseInt(d.gi[prediction_name][giIndex]['start']);
                    var endPosition = parseInt(d.gi[prediction_name][giIndex]['end']);

                    if ((endPosition - startPosition)>giFiltervalue) {
                        var rectpoints = self.scale(startPosition) + "," + (self.getSequenceModHeight() * i + GISIZE / 2) + " ";
                        rectpoints += self.scale(endPosition) + "," + (self.getSequenceModHeight() * i + GISIZE / 2) + " ";
                        rectpoints += self.scale(endPosition) + "," + (self.getSequenceModHeight() * i - GISIZE / 2) + " ";
                        rectpoints += self.scale(startPosition) + "," + (self.getSequenceModHeight() * i - GISIZE / 2) + " ";

                        var gi = giClassContainer.append("polygon")
                            .attr("points", rectpoints)
                            .attr("stroke-width", 1)
                            .attr("transform", "translate(" + 0 + "," + (GISIZE / 2) + ")");

                        // if color was given for this gi, then color the fill and stroke of this gi to the given color
                        if (d.gi[prediction_name][giIndex]['color'] != null) {
                            gi.attr("fill", d.gi[prediction_name][giIndex]['color'])
                                .attr("stroke", d.gi[prediction_name][giIndex]['color']);
                        }
                    }
                }
            }
        });

        //Add AMR genes to the SVG
        // TODO AMR legend
        var amrcontainer = seq.append("g")
            .attr("class", "amrs");
        var amrs = seq.each(function(d, i) {
            for (var amrIndex = 0; amrIndex < d.amr.length; amrIndex++){
                var amr = d.amr[amrIndex]
                var startPosition = parseInt(amr['start']);
                var endPosition = parseInt(amr['end']);
                var width = Math.abs(self.scale(endPosition) - self.scale(startPosition));
                var strand = (amr['strand'] == "+") ? true : false;

                // min ensures the angles on the trapezoids are at most 45 degrees
                // helps with visualization while zoomed in
                var rectpoints = self.scale(startPosition) + Math.min(+!strand * width / 3, GISIZE) + "," + (self.getSequenceModHeight() * i + GISIZE) + " ";
                rectpoints += self.scale(endPosition) - Math.min(+!strand * width / 3, GISIZE)+ "," + (self.getSequenceModHeight() * i + GISIZE) + " ";
                rectpoints += self.scale(endPosition) - Math.min(+strand * width / 3, GISIZE)  + "," + (self.getSequenceModHeight() * i) + " ";
                rectpoints += self.scale(startPosition) + Math.min(+strand * width / 3, GISIZE) + "," + (self.getSequenceModHeight() * i) + " ";

                amrcontainer.append("polygon")
                    .attr("points", rectpoints)
                    // Move + strand genes above the sequence, - strand below the sequence
                    .attr("transform", function() {
                        if (strand) {
                            return "translate(0," + ( -GISIZE ) + ")";
                        } else {
                            return "translate(0," + ( GISIZE ) + ")";
                        }
                    });
            }
        });

        //Add the brush for zooming and focusing
        var brush = d3.svg.brush()
            .x(self.scale)
            .on("brush", brushmove)
            .on("brushend", brushend);

        sequenceHolder.append("g")
            .attr("class", "brush")
            .call(brush)
            .selectAll('rect')
            .attr('height', this.containerHeight())
            .attr("transform","translate(0,"+(-0.7)*this.getSequenceModHeight()+")");

        function brushmove() {
            var extent = brush.extent();
        }

        function brushend() {
            var extent = brush.extent();
            self.setScale(extent[0],extent[1]);
            self.transition();
        }

        //Add the genes to the plot
        var geneFilterValue = self.getGeneFilterValue();
        if((self.scale.domain()[1]-self.scale.domain()[0])<geneFilterValue) {
            var geneContainer = sequenceHolder.append("g")
                .attr("class", "genes")
                .attr("transform", "translate(0," + (GENESIZE / 4 + GENESIZE/2) + ")");
            seq.each(function(d, i){
                $.ajax({
                    url: genomesGenesUrl + d.sequenceId + "?start=" + Math.round(self.scale.domain()[0]) + "&end=" + Math.round(self.scale.domain()[1]),
                    async: false,
                    success: function(data){
                        for (var geneIndex = 0; geneIndex < data['genes'].length; geneIndex++){
                            if (data['genes'][geneIndex]['strand'] == 1) {
                                var rectpoints = self.scale((data['genes'][geneIndex]['start'])) + "," + (self.getSequenceModHeight() * i) + " ";
                                rectpoints += self.scale((data['genes'][geneIndex]['end'])) + "," + (self.getSequenceModHeight() * i) + " ";
                                rectpoints += self.scale((data['genes'][geneIndex]['end'])) + "," + (self.getSequenceModHeight() * i - GENESIZE / 2 - 1) + " ";
                                rectpoints += self.scale((data['genes'][geneIndex]['start'])) + "," + (self.getSequenceModHeight() * i - GENESIZE / 2 - 1) + " ";
                            }

                            else if (data['genes'][geneIndex]['strand'] == -1) {
                                var rectpoints = self.scale((data['genes'][geneIndex]['start'])) + "," + (self.getSequenceModHeight() * i + GENESIZE / 2 - 1) + " ";
                                rectpoints += self.scale((data['genes'][geneIndex]['end'])) + "," + (self.getSequenceModHeight() * i + GENESIZE / 2 - 1) + " ";
                                rectpoints += self.scale((data['genes'][geneIndex]['end'])) + "," + (self.getSequenceModHeight() * i + GENESIZE - 1) + " ";
                                rectpoints += self.scale((data['genes'][geneIndex]['start'])) + "," + (self.getSequenceModHeight() * i + GENESIZE - 1) + " ";
                            }
                            else {
                                console.log("An error occured while trying to draw: " + d.genes[geneIndex['name']]);
                                continue;
                            }

                            var genename = data['genes'][geneIndex]['name'];

                            geneContainer.append("polygon")
                                .attr("points", rectpoints)
                                .attr("stroke-width", 1)
                                .append("title")
                                .text(function (d, i) {
                                    return genename;
                                });
                        }
                    },
                    headers: {
                        "Authorization": 'Token' + ' ' + getAuthenticationCookie()
                    }
                })
            });
        }

        //Adds the xAvis TODO Need a different implementation for IslandViewer
        var xAxis = d3.svg.axis()
            .scale(self.scale)
            .orient("bottom")
            .tickSize(-10)
            .tickFormat(d3.format("s"));

        visContainer.append("g")
            .attr("class","xAxis")
            .attr("transform", "translate(0," + ((seqOrder.length*self.getSequenceModHeight()+XAXISHEIGHT/6) + ")"))
            .call(xAxis)
            .append("rect")
            .attr("width",this.visualizationWidth())
            .attr("height",2)
            .attr("transform","translate(0,"+0+")");

        // Add the genome names to the visualization
        var textContainer = svg.append("g")
            .attr("class","sequenceLabels");

        var text = textContainer.selectAll("text")
            .data(sequencedata)
            .enter()
            .append("text");

        //Add SVG Text Element Attributes
        var textLabels = text.attr("y", function(d,i){ return (i)*self.getSequenceModHeight()})
            .text(function(d){return d.shownName});

        textContainer.attr("transform","translate("+(TREECONTAINERWIDTH+TEXTPADDING)+","+(TEXTTOPPADDING+TOPPADDING)+")");

        //Aligns the viscontainer to the right to make room for other containers
        visContainer.attr("transform","translate("+LEFTPADDING+","+((GISIZE/2)+TOPPADDING)+")");

        //Change all colors gray if togglePrinterColors = true
        //Could be moved to when sequences are rendered but will end up with messier code
        if (this.isPrinterColors){
            this.setPrinterColors();
        }
    };

    this.togglePrinterColors = function(){
        this.isPrinterColors=!this.isPrinterColors;
    };

    this.setPrinterColors = function(){
        $(".sequences").attr("class","sequences print");
        $(".genomicIslands").attr("class","genomicIslands print");
        $(".genes").attr("class","genes print");
        $(".node").attr("class","nodes print");
    };

    return this;
}

function Backbone() {
    var self = this;
    this.sequences = [];
    this.backbone = [[]];

    // Retrieves the sequenceId from sequences given the name of a sequence.
    // Will return -1 if no sequence is found
    this.getSequenceIdFromName = function (name) {
        for (var index = 0; index < this.sequences.length; index++) {
            if (this.sequences[index]['sequenceName'] == (name)) {
                return this.sequences[index]['sequenceId'];
            }
        }
        return -1;
    };

    this.addSequence = function (sequenceId, sequenceSize, sequenceName, givenName) {
        sequenceName = sequenceName || null;
        givenName = givenName || null;
        seq = new Sequence(sequenceId, sequenceSize, sequenceName, givenName);
        this.sequences.push(seq);
        return seq
    };

    this.getSequences = function () {
        return self.sequences;
    };

    this.addHomologousRegion = function (seqId1, seqId2, start1, end1, start2, end2) {
        if (this.backbone[seqId1] === undefined) {
            this.backbone[seqId1] = [];
        }

        if (this.backbone[seqId1][seqId2] === undefined) {
            this.backbone[seqId1][seqId2] = [];
        }

        if (this.backbone[seqId2] === undefined) {
            this.backbone[seqId2] = [];
        }

        if (this.backbone[seqId2][seqId1] === undefined) {
            this.backbone[seqId2][seqId1] = [];
        }
        this.backbone[seqId1][seqId2].push(new HomologousRegion(start1, end1, start2, end2));
        this.backbone[seqId2][seqId1].push(new HomologousRegion(start2, end2, start1, end1));
    };

    this.retrieveHomologousRegions = function (seqId1, seqId2) {
        try {
            var homologousRegions = this.backbone[seqId1][seqId2];
        } catch (e) {
            return [];
        }
        if (homologousRegions === undefined) {
            return [];
        }
        else {
            return homologousRegions;
        }
    };

    // retrieve json from server and render
    this.retrieveJsonAndRender = function (multiVis, url, geneUrl) {
        var backbonereference = this;
        $.ajax({
            url: url,
            success: function (data) {
                for (var genomeIndex = 0; genomeIndex < data['genomes'].length; genomeIndex++) {
                    var currentseq = backbonereference.addSequence(genomeIndex, data['genomes'][genomeIndex]['length'],
                        data['genomes'][genomeIndex]['name'], data['genomes'][genomeIndex]['givenName']);
                    // Add GIs to appropriate sequence
                    for (var giIndex = 0; giIndex < data['genomes'][genomeIndex]['gis'].length; giIndex++) {
                        currentseq.addGI(data['genomes'][genomeIndex]['gis'][giIndex]);
                    }
                    // Add genes to appropriate sequence
                    for (var geneIndex = 0; geneIndex < data['genomes'][genomeIndex]['genes'].length; geneIndex++) {
                        currentseq.addGene(data['genomes'][genomeIndex]['genes'][geneIndex]);
                    }
                    // at this scale, individual scaling for sequences may not be usable...so used fixed scale
                    currentseq.updateScale(0, multiVis.getLargestSequenceSize(), multiVis.getLargestSequenceSize());
                }

                for (var sequenceIndex = 0; sequenceIndex < Object.keys(data['backbone']).length; sequenceIndex++) {
                    for (var regionIndex = 0; regionIndex < data['backbone'][Object.keys(data['backbone'])[sequenceIndex]].length; regionIndex++) {
                        backbonereference.addHomologousRegion(parseInt(Object.keys(data['backbone'])[sequenceIndex]), parseInt(Object.keys(data['backbone'])[sequenceIndex]) + 1,
                            data['backbone'][Object.keys(data['backbone'])[sequenceIndex]][regionIndex][0][0],
                            data['backbone'][Object.keys(data['backbone'])[sequenceIndex]][regionIndex][0][1],
                            data['backbone'][Object.keys(data['backbone'])[sequenceIndex]][regionIndex][1][0],
                            data['backbone'][Object.keys(data['backbone'])[sequenceIndex]][regionIndex][1][1]);
                    }
                }
                multiVis.newickData = Newick.parse(data['newick']);
                multiVis.newickRoot = multiVis.newickData;
                // if a gene end point is given then retrieve it async
                if (geneUrl != null) {
                    backbonereference.retrieveGeneDataAsync(geneUrl);
                }
                multiVis.render();
            }
        })
    };

    //load gene data
    this.retrieveGeneDataAsync = function (url) {
        var backbonereference = this;
        $.ajax({
            url: url,
            success: function (data) {
                for (var genomeIndex = 0; genomeIndex < backbonereference.sequences.length; genomeIndex++) {
                    var genomeName = backbonereference.sequences[genomeIndex].sequenceName;
                    for (var geneIndex = 0; geneIndex < data[genomeName].length; geneIndex++) {
                        backbonereference.sequences[genomeIndex].addGene(data[genomeName][geneIndex]);
                    }
                }
            }
        });
    };
}

// Object that holds homologous regions
function HomologousRegion(start1,end1,start2,end2){
    this.start1 = start1;
    this.end1 = end1;
    this.start2 = start2;
    this.end2 = end2;

    // Negative numbers are inversions so remove negative and swap start and end points
    if (this.start1 < 0){
        memory = this.start1;
        this.start1=this.end1*-1;
        this.end1=memory*-1;
    }

    // Negative numbers are inversions so remove negative and swap start and end points
    if (this.start2 < 0){
        memory = this.start2;
        this.start2=this.end2*-1;
        this.end2 = memory*-1;
    }

    return this;
}

//Object to hold sequences
function Sequence(sequenceId, sequenceSize, sequenceName, givenName){
    this.sequenceId = sequenceId;
    this.sequenceName = sequenceName;
    this.sequenceSize = sequenceSize;
    this.shownName = givenName || sequenceName;
    this.genes = [];
    this.gi = {};
    this.amr = [];
    this.scale = null;

    this.updateScale = function (start,end, containerwidth){
        this.scale = d3.scale.linear().domain([start,end]).range([0,containerwidth]);
    };

    this.getSequenceSize = function(){
        return this.sequenceSize;
    };

    this.addGI = function(giClass, giDict){
        if (this.gi.hasOwnProperty(giClass)) {
            this.gi.push(giDict);
        }
        else{
            this.gi[giClass] = giDict;
        }
    };

    this.addAMR = function(amrList) {
        this.amr.push.apply(this.amr, amrList);
    };

    this.addGene = function(geneDict) {
        this.genes.push(geneDict);
    };

    return this;
}