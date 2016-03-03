//Notes: seems appropriate to move Multivis.sequences to Backbone

function MultiVis(targetNode){
    var self = this;
    const SEQUENCEHEIGHT = 90;
    const CONTAINERWIDTH = 1115;
    const TREECONTAINERWIDTH = 140;
    const LEFTPADDING = 85+TREECONTAINERWIDTH;
    const GISIZE = 30;
    const GENESIZE = 17;
    const GIFILTERFACTOR = 0;
    const GENEFILTERFACTOR =60000;
    const SEQUENCEWIDTH=8;
    const RIGHTPADDING = 40;

    this.container = d3.select(targetNode);
    this.backbone = new Backbone();
    this.sequences = this.backbone.getSequences();
    this.scale = null;
    this.treeData = null;
    this.sequenceOrder = null;
    this.isPrinterColors = false;

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
        var windowSize = self.scale.domain()[1]-self.scale.domain()[0];
        return windowSize/self.getLargestSequenceSize()*GIFILTERFACTOR;
    };

    this.getGeneFilterValue = function(){
        var windowSize = self.scale.domain()[1]-self.scale.domain()[0];
        return windowSize/self.getLargestSequenceSize()*GENEFILTERFACTOR;
    };

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
        //return this.sequences.length*SEQUENCEHEIGHT-100;// The -100 fixes padding issues on islandviewer site, fix this later if required;
        return this.sequences.length*SEQUENCEHEIGHT;
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

    this.resetAndRenderRange = function(){
        this.setScale(0,this.getLargestSequenceSize());
        this.transition();
    };

    this.render = function (){
        this.container.select("svg").remove();
        if (self.scale == null){
            self.setScale(0,this.getLargestSequenceSize());
        }

        //Add the SVG
        var svg = this.container.append("svg")
            .attr("width",this.containerWidth())
            .attr("height",this.containerHeight());

        //Add the visualization container
        var visContainer = svg.append("g")
            .attr("class","visContainer");

        //Add the tree visualization container
        var treeContainer = svg.append("g")
            .attr("class","treeContainer")
            .attr("width",TREECONTAINERWIDTH)
            .attr("height",this.containerHeight())
            .attr("transform", "translate(" + 0 + "," + (-GISIZE/2+2.5) + ")");

        //Add the tree
        var cluster = d3.layout.cluster()
            .separation(function(a, b) { return (a.parent == b.parent ? 1 : 1 ) })
            .size([this.containerHeight(), TREECONTAINERWIDTH/1.5]);

        function elbow(d, i) {
            return "M" + d.source.y + "," + d.source.x
                + "V" + d.target.x + "H" + d.target.y;
        }

        var i = 0;
        root = this.treeData;
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
        var seqOrder = self.getSequenceOrder();

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
                var points = self.scale(this.sequences[i].scale(homologousRegions[j].start1))+","+SEQUENCEHEIGHT*i+" ";
                points += self.scale(this.sequences[i].scale(homologousRegions[j].end1))+","+SEQUENCEHEIGHT*i+" ";
                points += self.scale(this.sequences[i+1].scale(homologousRegions[j].end2))+","+SEQUENCEHEIGHT*(i+1)+" ";
                points += self.scale(this.sequences[i+1].scale(homologousRegions[j].start2))+","+SEQUENCEHEIGHT*(i+1)+" ";

                homolousRegion.append("polygon")
                    .attr("points",points)
                    .attr("stroke-width",1)
                    .append("title")
                    .text("["+homologousRegions[j].start1+","+homologousRegions[j].end1+"],"+
                        "["+homologousRegions[j].start2+","+homologousRegions[j].end2+"]");
            }
        }

        //Create the sequences container on the svg
        var seq = sequenceHolder.selectAll("sequencesAxis")
            .data(this.sequences)
            .enter()
            .append("g")
            .attr("class", "sequences");

        //Add the sequences to the SVG
        seq.append("rect")
            .attr("x",0)
            .attr("y",function (d, i){
                return (i*SEQUENCEHEIGHT)+"px";
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
            //TODO Do something with positive or reverse strand
            for (var giIndex=0;giIndex< d.gi.length;giIndex++){
                if ((d.gi[giIndex]['end']-d.gi[giIndex]['start'])>giFiltervalue) {
                    var rectpoints = self.scale((d.gi[giIndex]['start'])) + "," + (SEQUENCEHEIGHT * i + GISIZE / 2) + " ";
                    rectpoints += self.scale((d.gi[giIndex]['end'])) + "," + (SEQUENCEHEIGHT * i + GISIZE / 2) + " ";
                    rectpoints += self.scale((d.gi[giIndex]['end'])) + "," + (SEQUENCEHEIGHT * i - GISIZE / 2) + " ";
                    rectpoints += self.scale((d.gi[giIndex]['start'])) + "," + (SEQUENCEHEIGHT * i - GISIZE / 2) + " ";

                    genomicIslandcontainer.append("polygon")
                        .attr("points", rectpoints)
                        .attr("stroke-width", 1)
                        .attr("transform","translate("+0+","+1+")");
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
            .attr("transform","translate(0,"+(-0.7)*SEQUENCEHEIGHT+")");

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
        var genes = seq.each(function(d, i){
            var geneContainer = sequenceHolder.append("g")
                .attr("class","genes")
                .attr("transform","translate(0,"+GENESIZE/4+")");
            for (var geneIndex=0;geneIndex< d.genes.length;geneIndex++){
                if((d.genes[geneIndex]['end']- d.genes[geneIndex]['start'])>geneFilterValue) {
                    var rectpoints = self.scale((d.genes[geneIndex]['start'])) + "," + (SEQUENCEHEIGHT * i + GENESIZE / 2) + " ";
                    rectpoints += self.scale((d.genes[geneIndex]['end'])) + "," + (SEQUENCEHEIGHT * i + GENESIZE / 2) + " ";
                    rectpoints += self.scale((d.genes[geneIndex]['end'])) + "," + (SEQUENCEHEIGHT * i - GENESIZE / 2) + " ";
                    rectpoints += self.scale((d.genes[geneIndex]['start'])) + "," + (SEQUENCEHEIGHT * i - GENESIZE / 2) + " ";

                    var genename = d.genes[geneIndex]['name'];

                    geneContainer.append("polygon")
                        .attr("points", rectpoints)
                        .attr("stroke-width", 1)
                        .append("title")
                        .text(function(d, i) { return genename; });
                }
            }
        });

        //Adds the xAvis TODO Need a different implementation for IslandViewer
        var xAxis = d3.svg.axis()
            .scale(self.scale)
            .orient("bottom")
            .tickSize(-10)
            .tickFormat(d3.format("s"));

        visContainer.append("g")
            .attr("class","xAxis")
            .attr("transform", "translate(0," + (SEQUENCEHEIGHT*(self.sequences.length-0.65)+(GISIZE/2)) + ")")
            .call(xAxis)
            .append("rect")
            .attr("width",this.visualizationWidth())
            .attr("height",2)
            .attr("transform","translate(0,"+0+")");

        //Add the SVG Text Element to the svgContainer
        //Used to test if tree and appropriate sequence is mapping correctly

        /*
        var textContainer = svg.append("g")
            .attr("class","sequenceLabels");

        var text = textContainer.selectAll("text")
            .data(this.sequences)
            .enter()
            .append("text");

        //Add SVG Text Element Attributes
        var textLabels = text.attr("y", function(d,i){ return i*SEQUENCEHEIGHT})
            .text(function(d){return d.sequenceName});

        textContainer.attr("transform","translate("+TREECONTAINERWIDTH+",25)");
        */

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

    this.addSequence = function (sequenceId, sequenceSize, sequenceName) {
        sequenceName = sequenceName || null;
        if (sequenceName != null){
            seq = new Sequence(sequenceId,sequenceSize,sequenceName)
        }
        else {
            seq = new Sequence(sequenceId, sequenceSize);
        }
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

    //Parses and then renders a backbone file in the target multivis object
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

                    //Dont Load "Matches" that do not contain a homologous region TODO theres a bug here
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
function Sequence(sequenceId, sequenceSize, sequenceName){
    this.sequenceId = sequenceId;
    this.sequenceName = sequenceName;
    this.sequenceSize = sequenceSize;
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