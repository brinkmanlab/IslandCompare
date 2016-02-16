//Notes: seems appropriate to move Multivis.sequences to Backbone

function MultiVis(targetNode){
    var self = this;
    const SEQUENCEHEIGHT = 150;
    const CONTAINERWIDTH = 1115;
    const LEFTPADDING = 85;
    // const NUMBERAXISTICKS = 5; BROKEN FOR CURRENT IMPLEMENTATION

    this.container = d3.select(targetNode);
    this.backbone = new Backbone();
    this.sequences = this.backbone.getSequences();

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
            return CONTAINERWIDTH-LEFTPADDING;
        }
        else{
            return this.container.node().getBoundingClientRect().width-LEFTPADDING;
        }
    };

    this.containerHeight = function() {
        return this.sequences.length*SEQUENCEHEIGHT-100;// The -100 fixes padding issues on islandviewer site, fix this later if required;
    };

    this.updateSequenceVisualization= function(sequenceIndex, newstart, newend){
        try {
            self.getSequence(sequenceIndex).updateScale(newstart, newend, this.visualizationWidth());
            self.transition();
        }catch(e){
            //Chart may not have been initialized yet (call parseandrender?)
        }
    };

    //Readjusts the graph for updated sequence domains, (improve later, currently just re-renders graph)
    this.transition = function(){
        this.container.select("svg").remove();
        this.render();
    };

    this.render = function (){
        this.container.select("svg").remove();
        var scale = d3.scale.linear().domain([0,this.getLargestSequenceSize()]).range([0,this.visualizationWidth()]);

        //Add the SVG
        var svg = this.container.append("svg")
            .attr("width",this.containerWidth())
            .attr("height",this.containerHeight());

        var visContainer = svg.append("g")
            .attr("class","visContainer");

        //Draw Homologous Region Lines
        var lines = [];

        for (var i=0; i<this.sequences.length-1; i++){
            var seqlines = visContainer.append("g")
                .attr("class","all-homolous-regions");
            var homologousRegions = (this.backbone.retrieveHomologousRegions(i,i+1));

            //Homologous Region polygons (shaded regions)
            for (var j=0;j<homologousRegions.length;j++){
                var homolousRegion = seqlines.append("g")
                    .attr("class", "homologous-region");

                //Build Shaded Polygon For Homologous Region
                var points = this.sequences[i].scale(homologousRegions[j].start1)+","+SEQUENCEHEIGHT*i+" ";
                points += this.sequences[i].scale(homologousRegions[j].end1)+","+SEQUENCEHEIGHT*i+" ";
                points += this.sequences[i+1].scale(homologousRegions[j].end2)+","+SEQUENCEHEIGHT*(i+1)+" ";
                points += this.sequences[i+1].scale(homologousRegions[j].start2)+","+SEQUENCEHEIGHT*(i+1)+" ";

                homolousRegion.append("polygon")
                    .attr("points",points)
                    .attr("stroke","#808080")
                    .attr("stroke-width",1)
                    .attr("fill","#C0C0C0")
                    .append("title")
                    .text("["+homologousRegions[j].start1+","+homologousRegions[j].end1+"],"+"["+homologousRegions[j].start2+","+homologousRegions[j].end2+"]");
            }
        }

        /*
         //BROKEN FOR CURRENT IMPLEMENTATION
         //Prepare the xAxis Fix this for new scales
         var xAxis = d3.visContainer.axis().scale(scale)
         .orient("bottom")
         .ticks(NUMBERAXISTICKS)
         .tickFormat(d3.format("s"));

         //Add the xAxis to the SVG
         var sequenceaxis = visContainer.selectAll("sequencesAxis")
         .data(this.sequences)
         .enter()
         .append("g")
         .attr("class", "x axis")
         .call(xAxis);


         //Modify the attributes of the axis on the SVG
         sequenceaxis.attr("x",0)
         .attr("transform",function (d, i){
         return "translate(0,"+(i*SEQUENCEHEIGHT)+")";
         });
         */

        //Create the sequences container on the svg
        var seq = visContainer.selectAll("sequencesAxis")
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
            .attr("height", 2)
            .attr("width", function (d){
                return d.scale(d.getSequenceSize());
            });

        //Add the gis to the SVG
        var gis = seq.each(function(d, i){
            var genomicIslandcontainer = seq.append("g")
                .attr("class","genomicIslands");
            for (var giIndex=0;giIndex< d.gi.length;giIndex++){
                genomicIslandcontainer.append("rect")
                    .attr("class","gi")
                    //TODO Do something with positive or reverse strand
                    .attr("width", Math.abs(d.scale(d.gi[giIndex]['end']- d.gi[giIndex]['start'])))
                    .attr("height", 8)
                    .attr("style","fill:rgb(0,0,255)")
                    .attr("transform","translate("+ d.scale(d.gi[giIndex]['start']) +","
                        +i*SEQUENCEHEIGHT+")");
            }
        });

        //Add the SVG Text Element to the svgContainer
        var textContainer = svg.append("g")
            .attr("class","sequenceLabels");

        var text = textContainer.selectAll("text")
            .data(this.sequences)
            .enter()
            .append("text");

        //Add SVG Text Element Attributes
        var textLabels = text.attr("y", function(d,i){ return i*SEQUENCEHEIGHT})
            .text(function(d){return d.sequenceName});

        visContainer.attr("transform","translate("+LEFTPADDING+",20)");
        textContainer.attr("transform","translate(0,25)");
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
    this.parseAndRenderBackbone= function(backboneFile,multiVis,genomeData,isFixedScale){
        var backbonereference = this;
        genomeData = genomeData || null;
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
                }
                else {
                    var currentseq = backbonereference.addSequence(i, largestBase[i]);
                }
                currentseq.updateScale(0,largestBase[i],multiVis.visualizationWidth());
            }

            //If fixedscale variable is set, than scaling is not done for individual sequences
            //TODO refactor this and above statement, will reduce calculations done
            if (!isFixedScale) {
                for (var j = 0; j < numberSequences; j++) {
                    backbonereference.sequences[j].updateScale(0, multiVis.getLargestSequenceSize(), multiVis.visualizationWidth());
                }
            }
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

    return this;
}




