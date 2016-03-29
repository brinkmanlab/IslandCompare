//Notes: seems appropriate to move Multivis.sequences to Backbone
//Dependencies: jQuery, D3.js

function MultiVis(targetNode){
    var self = this;
    const SEQUENCEHEIGHT = 20;
    const CONTAINERWIDTH = null;
    const TREECONTAINERWIDTH = 195;
    const TREETOPPADDING = 2;
    const LEFTPADDING = 85+TREECONTAINERWIDTH;
    const GISIZE = 8;
    const GENESIZE = 5;
    const GIFILTERFACTOR = 8000;
    const GENEFILTERFACTOR =400000;
    const SEQUENCEWIDTH=8;
    const RIGHTPADDING = 40;

    this.container = d3.select(targetNode);
    this.backbone = new Backbone();
    this.sequences = this.backbone.getSequences();
    this.scale = null;
    this.treeData = null;
    this.sequenceOrder = null;
    this.isPrinterColors = false;
    this.verticalScrollVal = 0;
    this.treeRoot = null;

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
            this.verticalScrollVal = (height/numberSequences)-SEQUENCEHEIGHT;
        }
    };

    // Returns the scale modifed height of the graph
    this.getSequenceModHeight = function(){
        return SEQUENCEHEIGHT + this.verticalScrollVal;
    };

    // Returns the scale modifed padding of the tree to get the alignment correct
    this.getTreeModPadding = function(){
        return TREETOPPADDING - (this.verticalScrollVal/2);
    };

    // TODO improve this
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

    // This sets the horizontal ranges of the graph. Where start is the base with the smallest position in the graph,
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
        this.treeRoot = this.treeData;
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
            .attr("height",(this.containerHeight()+this.getSequenceModHeight()*2));

        //Add the visualization container
        var visContainer = svg.append("g")
            .attr("class","visContainer");

        //Add the tree visualization container
        var treeContainer = svg.append("g")
            .attr("class","treeContainer")
            .attr("width",TREECONTAINERWIDTH)
            .attr("height",this.containerHeight())
            .attr("transform", "translate(" + 0 + "," + (this.getTreeModPadding()) + ")");

        //Add the tree
        var cluster = d3.layout.cluster()
            .separation(function(a, b) { return (a.parent == b.parent ? 1 : 1 ) })
            .size([seqOrder.length*(this.getSequenceModHeight()), TREECONTAINERWIDTH/1.5]);

        function elbow(d, i) {
            return "M" + d.source.y + "," + d.source.x
                + "V" + d.target.x + "H" + d.target.y;
        }

        var i = 0;
        root = this.treeRoot;
        update(root);

        function update(source) {
            // Compute the new tree layout.
            var nodes = cluster.nodes(root).reverse(),
                links = cluster.links(nodes);

            // Declare the nodes
            var node = treeContainer.selectAll("g.node")
                .data(nodes, function(d) { return d.id || (d.id = ++i); });

            // Enter the nodes.
            var nodeEnter = node.enter().append("g")
                .attr("class", "node")
                .attr("transform", function(d) {
                    return "translate(" + d.y + "," + d.x + ")"; });

            nodeEnter.append("circle")
                .attr("r", 5)
                .style("fill", function(d) { return "lightsteelblue"})
                .on('click',function(d){
                    //Sets the root of the tree to the clicked node
                    self.treeRoot = d;
                    //Set the order of the sequences (as linear plot and tree should match)
                    //Do a post-order traversal on tree looking for leaves.
                    var tree = [];
                    traverseTreeForOrder = function(node){
                        if (node['children']!=null){
                            for (var nodeIndex=0;nodeIndex<node['children'].length;nodeIndex++){
                                output = traverseTreeForOrder(node['children'][nodeIndex]);
                                if (output != null) {
                                    tree.push(output);
                                }
                            }
                        }
                        else{
                            return self.backbone.getSequenceIdFromName(node['name']);
                        }
                        return null;
                    };
                    traverseTreeForOrder(self.treeRoot);
                    //Set the sequenceOrder of the tree to the sequences obtained from traversal of the new treeRoot
                    self.sequenceOrder = tree;
                    self.transition();
                });

            // Adds the genome ids to the tree, (used to test if matching sequences correctly)
            /*
            nodeEnter.append("text")
                .attr("x", function(d) {
                    return d.children || d._children ? -13 : 13; })
                .attr("dy", ".35em")
                .attr("text-anchor", function(d) {
                    return d.children || d._children ? "end" : "start"; })
                .text(function(d) {
                    if (d.name==="None"){
                        return ''
                    }
                    else {
                        return d.name;
                    }
                })
                .style("fill-opacity", 1);
            */

            // Declare the links¦
            var link = treeContainer.selectAll("path.link")
                .data(links, function(d) { return d.target.id; });

            // Enter the links.
            link.enter().insert("path", "g")
                .attr("class", "link")
                .attr("d", elbow);
        }
        //Holds the linear plot visualization except the scale to prevent clipping/overflow problems
        var sequenceHolder = visContainer.append("svg")
            .attr("width",this.visualizationWidth())
            .append("g")
            .attr("transform","translate("+ 0 +","
                +GISIZE/2+")");

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
                var points = self.scale(this.sequences[i].scale(homologousRegions[j].start1))+","+this.getSequenceModHeight()*i+" ";
                points += self.scale(this.sequences[i].scale(homologousRegions[j].end1))+","+this.getSequenceModHeight()*i+" ";
                points += self.scale(this.sequences[i+1].scale(homologousRegions[j].end2))+","+this.getSequenceModHeight()*(i+1)+" ";
                points += self.scale(this.sequences[i+1].scale(homologousRegions[j].start2))+","+this.getSequenceModHeight()*(i+1)+" ";

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
            for (var giIndex=0;giIndex< d.gi.length;giIndex++){
                if ((d.gi[giIndex]['end']-d.gi[giIndex]['start'])>giFiltervalue) {
                    var rectpoints = self.scale((d.gi[giIndex]['start'])) + "," + (self.getSequenceModHeight() * i + GISIZE / 2) + " ";
                    rectpoints += self.scale((d.gi[giIndex]['end'])) + "," + (self.getSequenceModHeight() * i + GISIZE / 2) + " ";
                    rectpoints += self.scale((d.gi[giIndex]['end'])) + "," + (self.getSequenceModHeight() * i - GISIZE / 2) + " ";
                    rectpoints += self.scale((d.gi[giIndex]['start'])) + "," + (self.getSequenceModHeight() * i - GISIZE / 2) + " ";

                    genomicIslandcontainer.append("polygon")
                        .attr("points", rectpoints)
                        .attr("stroke-width", 1)
                        .attr("transform","translate("+0+","+(GISIZE/2)+")");
                }
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
            var genes = seq.each(function (d, i) {
                var geneContainer = sequenceHolder.append("g")
                    .attr("class", "genes")
                    .attr("transform", "translate(0," + (GENESIZE / 4 + GENESIZE/2) + ")");
                for (var geneIndex = 0; geneIndex < d.genes.length; geneIndex++) {
                    //Only show genes if window is smaller than geneFilterValue
                    var rectpoints = self.scale((d.genes[geneIndex]['start'])) + "," + (self.getSequenceModHeight() * i + GENESIZE / 2) + " ";
                    rectpoints += self.scale((d.genes[geneIndex]['end'])) + "," + (self.getSequenceModHeight() * i + GENESIZE / 2) + " ";
                    rectpoints += self.scale((d.genes[geneIndex]['end'])) + "," + (self.getSequenceModHeight() * i - GENESIZE / 2) + " ";
                    rectpoints += self.scale((d.genes[geneIndex]['start'])) + "," + (self.getSequenceModHeight() * i - GENESIZE / 2) + " ";

                    var genename = d.genes[geneIndex]['name'];

                    geneContainer.append("polygon")
                        .attr("points", rectpoints)
                        .attr("stroke-width", 1)
                        .append("title")
                        .text(function (d, i) {
                            return genename;
                        });
                }
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
            .attr("transform", "translate(0," + (self.containerHeight() + ((this.sequences.length-3.6)/this.sequences.length)*self.getSequenceModHeight()) + ")")
            .call(xAxis)
            .append("rect")
            .attr("width",this.visualizationWidth())
            .attr("height",2)
            .attr("transform","translate(0,"+0+")");

        //Add the SVG Text Element to the svgContainer (Which shows the genome id)
        //Used to test if tree and appropriate sequence is mapping correctly

        /*
        var textContainer = svg.append("g")
            .attr("class","sequenceLabels");

        var text = textContainer.selectAll("text")
            .data(this.sequences)
            .enter()
            .append("text");

        //Add SVG Text Element Attributes
        var textLabels = text.attr("y", function(d,i){ return i*self.getSequenceModHeight()})
            .text(function(d){return d.sequenceName});

        textContainer.attr("transform","translate("+TREECONTAINERWIDTH+",25)");
        */

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

        textContainer.attr("transform","translate("+(TREECONTAINERWIDTH-50)+","+17+")");

        //Aligns the viscontainer to the right to make room for other containers
        visContainer.attr("transform","translate("+LEFTPADDING+","+(GISIZE/2)+")");

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
    };

    return this;
}

function Backbone(){
    var self = this;
    this.sequences = [];
    this.backbone = [[]];

    // Retrieves the sequenceId from sequences given the name of a sequence.
    // Will return -1 if no sequence is found
    this.getSequenceIdFromName=function(name){
        for (var index=0;index<this.sequences.length;index++){
            if (this.sequences[index]['sequenceName']==(name)){
                return this.sequences[index]['sequenceId'];
            }
        }
        return -1;
    };

    this.addSequence = function (sequenceId, sequenceSize, sequenceName, givenName) {
        sequenceName = sequenceName || null;
        givenName = givenName || null;
        seq = new Sequence(sequenceId,sequenceSize,sequenceName, givenName);
        this.sequences.push(seq);
        return seq
    };

    this.getSequences = function(){
        return self.sequences;
    };

    this.addHomologousRegion = function(seqId1,seqId2,start1,end1,start2,end2){
        if (this.backbone[seqId1] ===undefined){
            this.backbone[seqId1] = [];
        }

        if (this.backbone[seqId1][seqId2]===undefined){
            this.backbone[seqId1][seqId2] = [];
        }

        if (this.backbone[seqId2] ===undefined){
            this.backbone[seqId2] = [];
        }

        if (this.backbone[seqId2][seqId1]===undefined){
            this.backbone[seqId2][seqId1] = [];
        }
        this.backbone[seqId1][seqId2].push(new HomologousRegion(start1,end1,start2,end2));
        this.backbone[seqId2][seqId1].push(new HomologousRegion(start2,end2,start1,end1));
    };

    this.retrieveHomologousRegions = function(seqId1,seqId2){
        try {
            var homologousRegions = this.backbone[seqId1][seqId2];
        } catch(e){
            return [];
        }
        if (homologousRegions === undefined){
            return [];
        }
        else{
            return homologousRegions;
        }
    };

    // retrieve json from server and render
    this.retrieveJsonAndRender=function(url,multiVis){
        var backbonereference = this;
        $.ajax({
            url:url,
            success: function(data){
                for (var genomeIndex=0;genomeIndex<data['genomes'].length;genomeIndex++){
                    var currentseq = backbonereference.addSequence(genomeIndex,data['genomes'][genomeIndex]['length'],
                        data['genomes'][genomeIndex]['name'], data['genomes'][genomeIndex]['givenName']);
                    // Add GIs to appropriate sequence
                    for (var giIndex=0;giIndex<data['genomes'][genomeIndex]['gis'].length;giIndex++){
                        currentseq.addGI(data['genomes'][genomeIndex]['gis'][giIndex]);
                    }
                    // Add genes to appropriate sequence
                    for (var geneIndex=0;geneIndex<data['genomes'][genomeIndex]['genes'].length;geneIndex++){
                        currentseq.addGene(data['genomes'][genomeIndex]['genes'][geneIndex]);
                    }
                    // at this scale, individual scaling for sequences may not be usable...so used fixed scale
                    currentseq.updateScale(0,multiVis.getLargestSequenceSize(), multiVis.getLargestSequenceSize());
                }

                for (var sequenceIndex=0;sequenceIndex<Object.keys(data['backbone']).length;sequenceIndex++){
                    for (var regionIndex=0;regionIndex<data['backbone'][Object.keys(data['backbone'])[sequenceIndex]].length;regionIndex++){
                        backbonereference.addHomologousRegion(parseInt(Object.keys(data['backbone'])[sequenceIndex]),parseInt(Object.keys(data['backbone'])[sequenceIndex])+1,
                            data['backbone'][Object.keys(data['backbone'])[sequenceIndex]][regionIndex][0][0],
                            data['backbone'][Object.keys(data['backbone'])[sequenceIndex]][regionIndex][0][1],
                            data['backbone'][Object.keys(data['backbone'])[sequenceIndex]][regionIndex][1][0],
                            data['backbone'][Object.keys(data['backbone'])[sequenceIndex]][regionIndex][1][1]);
                    }
                }
                multiVis.treeData = data['tree'];
                multiVis.treeRoot = multiVis.treeData;
                multiVis.render();
            }
        })
    };

    //Parses and then renders a backbone file and other json objects in the target multivis object
    //TODO: This function is broken. It is also less efficient than retrieveJSONandRender. Though it is used in islandviewer, fix this
    this.parseAndRenderBackbone= function(backboneFile,multiVis,genomeData,treeData,isFixedScale){
        var backbonereference = this;
        genomeData = genomeData || null;
        multiVis.treeData = treeData || null;
        d3.tsv(backboneFile, function(data){
            var numberSequences = (Object.keys(data[0]).length)/2;

            var choicelist = [];
            for (var seq1 = 0; seq1 < numberSequences-1; seq1++){
                for (var seq2 = 1; seq2< numberSequences; seq2++){
                    choicelist.push([seq1,seq2]);
                }
            }
            var largestBase = [];
            for (var j=0; j<numberSequences; j++){
                largestBase[j] = 0;
            }

            for (var row=1; row<data.length; row++){
                for (var choice=0; choice<choicelist.length; choice++) {
                    for (var k=0; k<largestBase.length; k++){
                        if (Number(data[row]["seq"+k+"_rightend"]) > largestBase[k]){
                            largestBase[k] = Number(data[row]["seq"+k+"_rightend"]);
                        }
                    }

                    //Dont Load "Matches" that do not contain a homologous region
                    if (data[row]["seq" + choicelist[choice][0] + "_rightend"] == 0 || data[row]["seq" + choicelist[choice][1] + "_rightend"] == 0){
                        continue;
                    }

                    backbonereference.addHomologousRegion( choicelist[choice][0],  choicelist[choice][1],
                        data[row]["seq" + choicelist[choice][0] + "_leftend"],
                        data[row]["seq" + choicelist[choice][0] + "_rightend"],
                        data[row]["seq" + choicelist[choice][1] + "_leftend"],
                        data[row]["seq" + choicelist[choice][1] + "_rightend"]);
                }
            }
            for (var i=0; i<numberSequences; i++){
                if (genomeData != null){
                    var currentseq = backbonereference.addSequence(i, largestBase[i],genomeData[i]['name']);
                    for (var arrayIndex=0;arrayIndex<genomeData[i]['gis'].length;arrayIndex++) {
                        currentseq.addGI(genomeData[i]['gis'][arrayIndex]);
                    }
                    for (var geneIndex=0;geneIndex<genomeData[i]['genes'].length;geneIndex++){
                        currentseq.addGene(genomeData[i]['genes'][geneIndex]);
                    }
                }
                else {
                    var currentseq = backbonereference.addSequence(i, largestBase[i]);
                }
                // If isFixedScale is true, then scaling is not done for individual sequences
                if (isFixedScale) {
                    currentseq.updateScale(0, multiVis.getLargestSequenceSize(),multiVis.getLargestSequenceSize());
                }
                else {
                    currentseq.updateScale(0, largestBase[i], multiVis.visualizationWidth());
                }
            }

            //Set the order of the sequences (as linear plot and tree should match)
            //Do a post-order traversal on tree looking for leaves.
            var treeOrder = [];
            traverseTreeForOrder = function(node){
                if (node['children']!=null){
                    for (var nodeIndex=0;nodeIndex<node['children'].length;nodeIndex++){
                        output = traverseTreeForOrder(node['children'][nodeIndex]);
                        if (output != null) {
                            treeOrder.push(output);
                        }
                    }
                }
                else{
                    return node['name'];
                }
                return null;
            };
            traverseTreeForOrder(treeData);
            multiVis.setSequenceOrderFromNames(treeOrder);
            multiVis.reorderSequences();
            multiVis.render();
        });
    }
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
    this.gi = [];
    this.scale = null;

    this.updateScale = function (start,end, containerwidth){
        this.scale = d3.scale.linear().domain([start,end]).range([0,containerwidth]);
    };

    this.getSequenceSize = function(){
        return this.sequenceSize;
    };

    this.addGI = function(giDict){
        this.gi.push(giDict);
    };

    this.addGene = function(geneDict) {
        this.genes.push(geneDict);
    };

    return this;
}