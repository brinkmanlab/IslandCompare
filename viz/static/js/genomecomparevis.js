//Dependencies: jQuery, D3.js, d3.phylogram.js, newick.js, linearplot.css (for styling)

function MultiVis(targetNode){
    var self = this;
    const TOPPADDING = 24;
    const SEQUENCEHEIGHT = 20;
    const CONTAINERWIDTH = null;
    const TREECONTAINERWIDTH = 195;
    const TREETOPPADDING = -16+TOPPADDING*2;
    const LEFTPADDING = 185+TREECONTAINERWIDTH;
    const TEXTTOPPADDING = 19+TOPPADDING;
    const TEXTPADDING = 30;
    const GIFILTERFACTOR = 4000;
    const GENEFILTERFACTOR =400000;
    const MINZOOM = GENEFILTERFACTOR / 100;
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
    this.cluster = null;

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
    this.getSequenceOrder = function() {
        if (this.sequenceOrder == null) {
            let treeOrder = this.traverseTreeForOrder(this.newickRoot);
            this.sequenceOrder = this.backbone.getIndicesFromIds(treeOrder);
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
        //Minimum zoom
        var diff = Math.abs(start - end);
        if (diff < MINZOOM) {
            //Center on region
            start = start + (diff / 2) - (MINZOOM / 2);
            end = start + MINZOOM;
        }
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
        this.container.select("svg").remove(); //Redundant?
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

    // Post-order tree traversal to get the order of the sequences
    this.traverseTreeForOrder = function (node) {
        var tree = [];
        if (node.branchset != null) {
            for (var nodeIndex = 0; nodeIndex < node.branchset.length; nodeIndex++) {
                var output = this.traverseTreeForOrder(node.branchset[nodeIndex]);
                Array.prototype.push.apply(tree, output);
            }
        }
        else {
            return [node['name']];
        }
        return tree;
    };

    //Renders the graph
    this.render = function (){
        this.container.select("svg").remove();

        // Gets the sequenceOrder of the graph
        var seqOrder = self.getSequenceOrder();

        // Modifes the spacing between sequences depending on the number of sequences to show
        self.updateVerticalScrollVal(seqOrder.length);

        //Add the SVG (Make sequence height 2 sequence higher than container height to fit svg TODO Refactor container height)
        var svg = this.container.append("svg")
            .attr("width",this.containerWidth())
            .attr("height",(this.containerHeight()+this.getSequenceModHeight()*2)+TOPPADDING);

        //Scale is set after adding svg as adding svg changes container width
        if (self.scale == null){
            self.setScale(0,this.getLargestSequenceSize());
        }

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
                    var seqIds = self.traverseTreeForOrder(self.newickRoot);
                    self.sequenceOrder = self.backbone.getIndicesFromIds(seqIds);
                    //Re-renders the graph
                    self.transition();
                }
            });
        }

        //Holds the linear plot visualization except the scale to prevent clipping/overflow problems
        var sequenceHolder = visContainer.append("svg")
            .attr("width",this.visualizationWidth())
            .append("g")
            .attr("transform","translate("+ 0 +","+(GISIZE / 2 + TOPPADDING)+")");

        //Draw Homologous Region Lines
        var lines = [];
        for (var i = 0; i < seqOrder.length - 1; i++){
            var seqlines = sequenceHolder.append("g")
                .attr("class","all-homolous-regions");

            var seqId = this.backbone.sequences[seqOrder[i]].sequenceId;
            var homologousRegions = this.backbone.retrieveHomologousRegions(seqId);

            //Homologous Region polygons (shaded regions)
            for (var j = 0; j < homologousRegions.length; j++){
                var homolousRegion = seqlines.append("g")
                    .attr("class", "homologous-region")
                    .attr("transform","translate("+ 0 +","+ SEQUENCEWIDTH/2 +")");

                //Build Shaded Polygon For Homologous Region
                var points = self.scale(homologousRegions[j].start1)+","+(this.getSequenceModHeight()*i+SEQUENCEWIDTH/2)+" ";
                points += self.scale(homologousRegions[j].end1)+","+(this.getSequenceModHeight()*i+SEQUENCEWIDTH/2)+" ";
                points += self.scale(homologousRegions[j].end2)+","+(this.getSequenceModHeight()*(i+1)-SEQUENCEWIDTH/2)+" ";
                points += self.scale(homologousRegions[j].start2)+","+(this.getSequenceModHeight()*(i+1)-SEQUENCEWIDTH/2)+" ";

                homolousRegion.append("polygon")
                    .attr("points",points)
                    .attr("stroke-width",1)
                    .append("title")
                    .text("["+homologousRegions[j].start1+","+homologousRegions[j].end1+"],"+
                        "["+homologousRegions[j].start2+","+homologousRegions[j].end2+"]");
            }
        }

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
            // Handle accidental zoom
            if (Math.abs(extent[0] - extent[1]) < 10) return;
            self.setScale(extent[0],extent[1]);
            self.transition();
        }

        //Create the sequence list from the ordered sequence list
        var sequenceData = [];
        for (var index=0;index<seqOrder.length;index++){
            sequenceData.push(this.sequences[seqOrder[index]]);
        }

        //Create the sequences container on the svg
        var seq = sequenceHolder.selectAll("sequencesAxis")
            .data(sequenceData)
            .enter()
            .append("g")
            .attr("class", "sequences")
            .attr("sequence", function(sequence) {
                return sequence.sequenceId;
            });

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
                var width = self.scale(d.scale(d.getSequenceSize()));
                if (width >= 0) {
                    return width;
                }
            });

        //Add GIs to each sequence
        seq.append("g").attr("class","genomicIslands").each(function(seqData, i){
            var methods = d3.select(this).selectAll("g")
                .data(Object.keys(seqData.gi))
                .enter().append("g")
                .attr("class", function(method) {
                    return method; //This assigns the tool used to generate [sigi, islandpath, merged]
                });
            methods.each(function(method) {
                d3.select(this).selectAll("polygon")
                    .data(seqData.gi[method].filter(function(gi) {
                        // Filter out the GIs that aren't within the current scale or are too small
                        var withinScale = gi.start < self.scale.domain()[1] && gi.end > self.scale.domain()[0];
                        var withinLength = gi.end - gi.start > self.getGIFilterValue();
                        return withinScale && withinLength;
                    }))
                    .enter().append("polygon")
                    .attr("points", function(gi) {
                        var startPosition = self.scale(parseInt(gi.start));
                        var endPosition = self.scale(parseInt(gi.end));
                        return startPosition + "," + (self.getSequenceModHeight() * i + GISIZE) + " " +
                                 endPosition + "," + (self.getSequenceModHeight() * i + GISIZE) + " " +
                                 endPosition + "," + (self.getSequenceModHeight() * i) + " " +
                               startPosition + "," + (self.getSequenceModHeight() * i) + " ";
                    })
                    .attr("start", function(gi) { return gi.start; })
                    .attr("end",function(gi) { return gi.end; })
                    .attr("fill", function(gi) {
                        if (gi.color != null) { return gi.color; }
                    })
                    .attr("stroke", function(gi) {
                        if (gi.color != null) { return gi.color; }
                    })
                    .attr("class", function(gi) {
                        let res = "";
                        if (gi.cluster != null) { res += "cluster-" + gi.cluster; }
                        if (res !== "") res += " ";
                        if (gi.parent) { res += "merged_component"; }
                        return res;
                    })
                    .on("mouseenter", function(gi) {
                        if(gi.cluster != null) { self.highlightCluster(".cluster-" + gi.cluster); }
                        return false;
                    })
                    .on("mouseleave", function(gi) {
                        if (gi.cluster != null) { self.unhighlightCluster(".cluster-" + gi.cluster); }
                        return false;
                    })
                    .on("click", function(gi) {
                        if (gi.cluster != null) {
                            try {
                                self.openClusterPage(self, gi.cluster, gi.color);
                            } catch(err) {
                                alert("Problem encountered when constructing cluster view: \n\n" + err);
                            }
                        }
                    })
                    .append("title").text(function(gi) {
                        if(gi.cluster != null) { return "Click to open GI cluster " + gi.cluster + " view"; }
                    });
            });
        });
        // Hide GIs if a cluster is selected
        if (self.cluster) {
            $("svg .genomicIslands polygon").not(".cluster-" + self.cluster).hide();
        }

        //Add AMR genes to each sequence
        seq.append("g").attr("class", "amrs").each(function(d, i) {
            // AMRs are given more detail at a smaller scale
            var details = ((self.scale.domain()[1] - self.scale.domain()[0]) < self.getGeneFilterValue() * 3);
            d3.select(this).selectAll("polygon")
                .data(d.amr.filter(function(d) {
                    // Filter out the AMRs that aren't within the current scale
                    return (d.start < self.scale.domain()[1] && d.end > self.scale.domain()[0]);
                }))
                .enter().append("polygon")
                .attr("points", function(d) {
                    var startPosition = self.scale(parseInt(d.start));
                    var endPosition = self.scale(parseInt(d.end));
                    var plus_strand = (d.strand === "+");
                    var polypoints = "";
                    // If scale is close, AMRs will be rendered as trapezoids
                    if (details) {
                        var width = Math.abs(endPosition - startPosition);
                        // min ensures the angles on the trapezoids are at most 45 degrees
                        polypoints += startPosition + Math.min(+!plus_strand * width / 3, GISIZE) + "," + (self.getSequenceModHeight() * i) + " ";
                        polypoints += endPosition - Math.min(+!plus_strand * width / 3, GISIZE) + "," + (self.getSequenceModHeight() * i) + " ";
                        polypoints += endPosition - Math.min(+plus_strand * width / 3, GISIZE) + "," + (self.getSequenceModHeight() * i - GISIZE / 2) + " ";
                        polypoints += startPosition + Math.min(+plus_strand * width / 3, GISIZE) + "," + (self.getSequenceModHeight() * i - GISIZE / 2);
                    // Else AMRs will be rectangles to improve appearance from a wider view
                    } else {
                        polypoints += startPosition + "," + (self.getSequenceModHeight() * i) + " ";
                        polypoints += startPosition + 1 + "," + (self.getSequenceModHeight() * i) + " ";
                        polypoints += startPosition + 1 + "," + (self.getSequenceModHeight() * i - GISIZE / 2) + " ";
                        polypoints += startPosition + "," + (self.getSequenceModHeight() * i - GISIZE / 2);
                    }
                    return polypoints;
                })
                .attr("transform", function(d) {
                    // Positive strand AMRs are in place, negative strand AMRs need to be moved to the other side
                    if (d.strand === "-") {
                        return "translate(0," + (3 * GISIZE / 2) + ")";
                    }
                });
        });

        //Add the genes to the plot
        var geneFilterValue = self.getGeneFilterValue();
        if((self.scale.domain()[1] - self.scale.domain()[0]) < geneFilterValue) {
            var geneContainer = sequenceHolder.append("g")
                .attr("class", "genes")
                .attr("transform", "translate(0," + (GENESIZE / 4 + GENESIZE/2) + ")");
            var seqCount = 0;
            var domain_start = Math.round(self.scale.domain()[0]);
            var domain_end = Math.round(self.scale.domain()[1]);
            seq.each(function(d, i){
                var genes = self.backbone.getSequences().find(function(s){return d.sequenceId == s.sequenceId;}).genes;
                seqCount++;
                for (var gene of genes.filter(function(gene){ return gene['start'] >= domain_start && gene['end'] <= domain_end; })){
                    if (gene['strand'] == '+' || gene['strand'] == '.') {
                        var rectpoints = self.scale((gene['start'])) + "," + (self.getSequenceModHeight() * i) + " ";
                        rectpoints += self.scale((gene['end'])) + "," + (self.getSequenceModHeight() * i) + " ";
                        rectpoints += self.scale((gene['end'])) + "," + (self.getSequenceModHeight() * i - GENESIZE / 2 - 1) + " ";
                        rectpoints += self.scale((gene['start'])) + "," + (self.getSequenceModHeight() * i - GENESIZE / 2 - 1) + " ";
                    }

                    else if (gene['strand'] == '-') {
                        var rectpoints = self.scale((gene['start'])) + "," + (self.getSequenceModHeight() * i + GENESIZE / 2 - 1) + " ";
                        rectpoints += self.scale((gene['end'])) + "," + (self.getSequenceModHeight() * i + GENESIZE / 2 - 1) + " ";
                        rectpoints += self.scale((gene['end'])) + "," + (self.getSequenceModHeight() * i + GENESIZE - 1) + " ";
                        rectpoints += self.scale((gene['start'])) + "," + (self.getSequenceModHeight() * i + GENESIZE - 1) + " ";
                    }

                    geneContainer.append("polygon")
                        .attr("points", rectpoints)
                        .attr("stroke-width", 1)
                        .attr("class", gene['type'])
                        .style("opacity", 0.0)
                        .append("title")
                            .text(function (d, i) {
                                var hover = "";
                                if(gene['gene']) {
                                    hover = hover.concat("Gene name: " + gene['gene'] + "\n");
                                }
                                if(gene['locus_tag']) {
                                    hover = hover.concat("Locus tag: " + gene['locus_tag'] + "\n");
                                }
                                if(gene['product']) {
                                    hover = hover.concat("Product: " + gene['product'])
                                }
                                if (hover == "") hover = "No provided annotations";
                                return hover;
                            });
                }
                // Wait until all genes have been added so they can be made visible simultaneously
                if (seqCount === seqOrder.length) {
                    $(".loadingGenes").hide();
                    $(".genesLegend-toggle").show();
                    geneContainer.selectAll("polygon")
                        .style("opacity", 1);
                }
            });
        } else {
            $(".genesLegend-toggle").hide();
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
            .data(sequenceData)
            .enter()
            .append("text");

        //Add SVG Text Element Attributes
        var textLabels = text.attr("y", function(d,i){ return (i)*self.getSequenceModHeight()})
            .text(function(d){return d.shownName});

        // Shorten long genbank file names that extend into the visContainer
        text.each(function() {
            if(this.getBoundingClientRect().width > LEFTPADDING - TREECONTAINERWIDTH - TEXTPADDING - 5) {
                d3.select(this)
                    .attr("textLength", LEFTPADDING - TREECONTAINERWIDTH - TEXTPADDING - 5)
                    .attr("lengthAdjust", "spacingAndGlyphs");
            }
        });

        textContainer.attr("transform","translate("+(TREECONTAINERWIDTH+TEXTPADDING)+","+(TEXTTOPPADDING+TOPPADDING)+")");

        //Aligns the viscontainer to the right to make room for other containers
        visContainer.attr("transform","translate("+LEFTPADDING+","+((GISIZE/2)+TOPPADDING)+")");

        //Change all colors gray if togglePrinterColors = true
        //Could be moved to when sequences are rendered but will end up with messier code
        if (this.isPrinterColors){
            this.setPrinterColors();
        }

        // Colour GIs by specified method
        this.predictorToggle();
        this.GIColourToggle();
    };

    // Resets the scale when the sidebar is toggled
    $(".sidebar-toggle").on("click", function() {
        setTimeout(function() {
            self.setScale(self.scale.domain()[0], self.scale.domain()[1]);
            self.transition();
        }, 500);
    });

    this.togglePrinterColors = function(){
        this.isPrinterColors=!this.isPrinterColors;
    };

    this.setPrinterColors = function(){
        $(".sequences").attr("class","sequences print");
        $(".genomicIslands").attr("class","genomicIslands print");
        $(".genes").attr("class","genes print");
        $(".node").attr("class","nodes print");
        $(".amrs").attr("class", "amrs print");
    };

    this.predictorToggle = function() {
        if ($(".GIColourBody :checkbox[name=sigi]").is(":checked")) {
            $(".sigi").show();
        } else {
            $(".sigi").hide();
        }
        if ($(".GIColourBody :checkbox[name=islandpath]").is(":checked")) {
            $(".islandpath").show();
        } else {
            $(".islandpath").hide();
        }
    };

    this.GIColourToggle = function() {
        if ($("input[name=GIColour]:checked").val() === "similarity") {
            $(".GIColourBody :checkbox").attr("disabled", true);
            $(".predictorLegend").hide();
            $(".merged").show();
            $(".merged_component").hide();
        } else {
            $(".GIColourBody :checkbox").removeAttr("disabled");
            $(".predictorLegend").show();
            $(".merged").hide();
            $(".merged_component").show();
        }
    };

    this.saveImage = function(mode) {
        // https://stackoverflow.com/a/38085847
        var svg = this.container.select("svg"),
            img = new Image(),
            serializer = new XMLSerializer(),
            width = svg.node().getBoundingClientRect().width,
            height = svg.node().getBoundingClientRect().height;

        // get the css from the stylesheet
        var style = "\nsvg {background-color: white;}\n";
        for (var i=0; i<document.styleSheets.length; i++) {
            var sheet = document.styleSheets[i];
            if (sheet.href) {
                var sheetName = sheet.href.split('/').pop();
                if (sheetName === 'linearplot.css') {
                    var rules = sheet.rules;
                    if (rules) {
                        for (var j = 0; j < rules.length; j++) {
                            style += (rules[j].cssText + '\n');
                        }
                    }
                }
            }
        }
        // add the css to the svg
        svg.insert("defs", ":first-child")
            .append("style")
            .attr('type', 'text/css')
            .html(style);

        var svgStr = serializer.serializeToString(svg.node());

        if(mode === "png") {
            var canvas = document.createElement("canvas");
            canvas.width = width;
            canvas.height = height;
            var context = canvas.getContext("2d");
            img.onload = function() {
                context.drawImage(img, 0, 0);
                var png = canvas.toDataURL("image/png");
                var link = document.createElement("a");
                link.setAttribute("href", png);
                link.setAttribute("download", "result.png");
                link.click();
            };
        }
        img.src = 'data:image/svg+xml;base64,'+window.btoa(svgStr);
        if(mode === "svg") {
            var link = document.createElement("a");
            link.setAttribute("href", img.src);
            link.setAttribute("download", "result.svg");
            link.click();
        }
    };

    this.highlightCluster = function(clusterClass) {
        $(clusterClass).attr("filter", "url(#shadow)");
        //$(".sigi, .islandpath").hide();
        var otherIslands = $("svg .genomicIslands polygon").not(clusterClass);
        otherIslands.attr("fill-opacity", "0.4");
        otherIslands.attr("stroke-opacity", "0");
        var genes = $(".CDS, .tRNA, .rRNA, .gene");
        genes.attr("fill-opacity", "0.4");
        genes.attr("stroke-opacity", "0");
    };

    this.unhighlightCluster = function(clusterClass) {
        $(clusterClass).removeAttr("filter");
        //$(".sigi, .islandpath").show();
        var otherIslands = $("svg .genomicIslands polygon").not(clusterClass);
        otherIslands.removeAttr("fill-opacity");
        otherIslands.removeAttr("stroke-opacity");
        var genes = $(".CDS, .tRNA, .rRNA, .gene");
        genes.removeAttr("fill-opacity");
        genes.removeAttr("stroke-opacity");
    };

    this.toggleClusterView = function(cluster, color) {
        $("svg .genomicIslands polygon").not(".cluster-" + cluster).toggle();
        $(".GIColour").children().toggle();
        $(".clusterButton").toggle();
        $("#clusterLegend").toggle();
        if (self.cluster == null) {
            self.cluster = cluster;
            $("#clusterLegend circle").attr("fill", color);
            $("#clusterLegendSpan").text("Cluster " + cluster + " Genomic Island");
        } else {
            self.cluster = null;
        }
    };

    this.zoomCluster = function() {
        var clusterClass = ".cluster-" + self.cluster;
        var start_array = [];
        var end_array = [];
        $(clusterClass).each(function() {
            start_array.push(parseInt($(this).attr("start")));
            end_array.push(parseInt($(this).attr("end")));
        });
        var start = Math.min.apply(null, start_array);
        var end = Math.max.apply(null, end_array);
        var buffer = parseFloat(end - start) / 80;
        self.setScale(start-buffer, end+buffer);
        self.transition();
        self.unhighlightCluster(clusterClass);
    };

    this.openClusterPage = function(multiVis, clusterNum, color) {
        var clusterDict = {};
        clusterDict.cluster = clusterNum;
        clusterDict.color = color;
        clusterDict.sequences = [];
        // Get all GIs in cluster
        var cluster = $(".cluster-" + clusterNum);
        // For each cluster record ID, start, end
        cluster.each(function() {
            var seqID = $(this).parents(".sequences").attr("sequence");
            var start = Number($(this).attr("start"));
            var end = Number($(this).attr("end"));
            var index = clusterDict.sequences.findIndex(i => i.id === seqID);
            if (index === -1) {
                clusterDict.sequences.push({"id": seqID, "islands": []});
                index = clusterDict.sequences.length - 1;
            }
            clusterDict.sequences[index].islands.push([start, end, seqID]);
        });
        // Assign sequence names
        for (var i = 0; i < multiVis.sequences.length; i++) {
            var currentSeq = multiVis.sequences[i];
            var index = clusterDict.sequences.findIndex(i => i.id === currentSeq.sequenceId);
            if (index !== -1) {
                clusterDict.sequences[index].name =  currentSeq.sequenceName;
                clusterDict.sequences[index].amr = currentSeq.amr;
                clusterDict.sequences[index].genes = currentSeq.genes;
            }
        }

        clusterDict.alignment = multiVis.backbone.backbone;

        var clusterPage = open("static/cluster.html");
        clusterPage.onload = ()=>{
            clusterPage.createClusterVisualization(clusterDict);
        };

    };

    return this;
}

function Backbone() {
    var self = this;
    this.sequences = [];
    this.backbone = {};

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

    this.getIndicesFromIds = function(ids) {
        let indices = [];
        for (id of ids) {
            for (var seq_i = 0; seq_i < this.sequences.length; seq_i++) {
                if (id === this.sequences[seq_i].sequenceId) {
                    indices.push(seq_i);
                    break;
                }
            }
        }
        return indices;
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

    this.addHomologousRegion = function (seqId, start1, end1, start2, end2) {
        if (this.backbone[seqId] === undefined) {
            this.backbone[seqId] = [];
        }
        this.backbone[seqId].push(new HomologousRegion(start1, end1, start2, end2));
    };

    this.retrieveHomologousRegions = function (seqId) {
        try {
            var homologousRegions = this.backbone[seqId];
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
        if (!this.gi.hasOwnProperty(giClass)) {
            this.gi[giClass] = [];
        }

        this.gi[giClass].push(giDict);
    };

    this.addAMR = function(amrList) {
        this.amr.push.apply(this.amr, amrList);
    };

    this.addGene = function(geneDict) {
        this.genes.push(geneDict);
    };

    return this;
}