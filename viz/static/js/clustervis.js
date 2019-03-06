function createClusterVisualization() {
    var SEQSTART = 115; // Islands will begin at this horizontal distance from the left - makes room for position.
    var SVGPADDING = 30; // Accounts for padding of each div
    var SEQPADDING = 25; // Vertical distance between each island in a genome
    var TEXTOFFSET = 15; // Adjustment to align position text with island
    var RECTHEIGHT = 10; // Height of the island rect elements

    var clusterDict = window.clusterDict;
    if (clusterDict === undefined) {
        clusterDict = JSON.parse(sessionStorage.getItem("clusterDict"));
    } else {
        // Store the clusterDict so it can be retrieved if the page is refreshed
        sessionStorage.setItem("clusterDict", JSON.stringify(clusterDict));
    }

    var cluster = clusterDict.cluster;
    $("#header").html("Cluster " + cluster);
    var color = clusterDict.color;
    var visWidth = $("#visualizationBody").width();
    var scale = getScale(clusterDict.sequences, visWidth - SVGPADDING - SEQSTART);

    // D3.js //
    var seq = d3.select("#visualizationBody")
        .selectAll("div")
        .data(clusterDict.sequences);

    var div = seq.enter().append("div");
    div.attr("class", function(seq) {
        return seq.name;
    });
    // Add genome name
    div.append("p").append("b")
        .html(function(seq) {
            return seq.name;
        });

    var svg = div.append("svg")
        .attr("height", function(seq) {
            return seq.islands.length * SEQPADDING;
        })
        .attr("width", visWidth - SVGPADDING);

    var rectGroup = svg.selectAll("g")
        .data(function(seq) {
            return seq.islands;
        });

    rectGroup.enter().append("g");
    // Add start and Stop text
    rectGroup.append("text")
        .text(function(d) {
            return d[0] + " - " + d[1];
        })
        .attr("y", function(d, i) {
            return (i * SEQPADDING) + TEXTOFFSET;
        });
    // Add the rect representing the genomic island
    rectGroup.append("rect")
        .attr("y", function(d, i) {
            return (i * SEQPADDING) + RECTHEIGHT / 2;
        })
        .attr("height", RECTHEIGHT)
        .attr("width", function(d) {
            return scale(Math.abs(d[1] - d[0]));
        })
        .attr("rx", 5)
        .attr("ry", 5)
        .attr("fill", color)
        .attr("transform", "translate(" + SEQSTART +",0)");
    // Add the genes
    rectGroup.append("g")
        .attr("class", "genes")
        .attr("transform", function(d, i) {
            return "translate(" + SEQSTART + "," + i * SEQPADDING + ")";
        }).each(function(island) {
            var geneGroup = d3.select(this);
            // For each island, an asynchronous call is made to fetch the genes within its boundaries.
            // This prevents the page from freezing up while it gets all the genes.

            // The genes are added to the island once they have been retreived from the database.
            geneGroup.selectAll("polygon")
                .data(data.genes)
                .enter().append("polygon")
                .attr("points", function(gene) {
                    var geneStart = gene.start - island[0];
                    // Cut off genes that extend beyond the island boundary
                    var geneEnd = Math.min(gene.end - island[0], island[1] - island[0]);
                    return scale(geneStart) + "," + RECTHEIGHT / 2 + "," +
                            scale(geneEnd) + "," + RECTHEIGHT / 2 + "," +
                            scale(geneEnd) + "," + RECTHEIGHT + "," +
                            scale(geneStart) + "," + RECTHEIGHT;
                })
                .attr("class", function(gene) {
                    return gene.type; // gene, rRNA, or tRNA
                })
                .attr("transform", function(gene) {
                    // Move negative strand genes to the bottom half of the island
                    if (gene.strand === -1) {
                        return "translate(0," + RECTHEIGHT / 2 +")";
                    }
                })
                .append("title")
                    .text(function(gene) {
                        var hover = "";
                        if(gene.gene) {
                            hover = hover.concat("Gene name: " + gene.gene + "\n");
                        }
                        if(gene.locus_tag) {
                            hover = hover.concat("Locus tag: " + gene.locus_tag + "\n");
                        }
                        if(gene.product) {
                            hover = hover.concat("Product: " + gene.product)
                        }
                        return hover;
                    });
        });
    // Add the AMR genes
    rectGroup.append("g").attr("class", "amrs")
        .attr("transform", function(d, i) {
            return "translate(" + SEQSTART + "," + i * SEQPADDING + ")"
        })
        .each(function(island, i) {
            var start = island[0];
            var end = island[1];
            var seq = island[2];
            var index = clusterDict.sequences.findIndex(i => i.id === seq);
            var amrs = clusterDict.sequences[index].amr;
            for (var i = 0; i < amrs.length; i++) {
                if (amrs[i].start < end && amrs[i].end > start) {
                    var amr_start = scale(amrs[i].start - start);
                    var amr_end = scale(amrs[i].end - start);
                    var neg_strand = (amrs[i].strand === "-");
                    d3.select(this).append("polygon")
                        .attr("points", function() {
                            // Returns a trapezoid shape.
                            // Shape depends on whether it is on negative strand and
                            // whether it is within the boundaries of the GI
                            var points = "";
                            if (neg_strand) {
                                if (amrs[i].start < start) {
                                    points += "0,0,0," + RECTHEIGHT / 2 + ",";
                                } else {
                                    points += amr_start + ",0," +
                                        (amr_start + RECTHEIGHT / 2) + "," + RECTHEIGHT / 2 + ",";
                                }
                                if (amrs[i].end > end) {
                                    points += scale(end - start) + "," + RECTHEIGHT / 2 +
                                        "," + scale(end - start) + ",0";
                                } else {
                                    points += (amr_end - RECTHEIGHT / 2) + "," + RECTHEIGHT / 2 + "," +
                                            amr_end + ",0";
                                }
                            } else {
                                if (amrs[i].start < start) {
                                    points += "0,0,0," + RECTHEIGHT / 2 + ",";
                                } else {
                                    points += (amr_start + RECTHEIGHT / 2) + ",0," + amr_start + "," + RECTHEIGHT / 2 + ",";
                                }
                                if (amrs[i].end > end) {
                                    points += scale(end - start) + "," + RECTHEIGHT / 2 +
                                        "," + scale(end - start) + ",0";
                                } else {
                                    points += amr_end + "," + RECTHEIGHT / 2 + "," +
                                        (amr_end - RECTHEIGHT / 2) + ",0";
                                }
                            }
                            return points;
                        })
                        .attr("transform", function() {
                            // Move negative strand AMR genes to the underside of the GI
                            if (neg_strand) {
                                return "translate(0, "+ 3 / 2 * RECTHEIGHT +")"
                            }
                        });
                }
            }
    });
}

function getScale(sequences, containerWidth) {
    var lengths = [];
    for (var seqIndex = 0; seqIndex < sequences.length; seqIndex++) {
        var islands = sequences[seqIndex].islands;
        for (var islandIndex = 0; islandIndex < islands.length; islandIndex++) {
            lengths.push(Math.abs(islands[islandIndex][1] - islands[islandIndex][0]));
        }
    }
    var max = Math.max.apply(null, lengths);
    return d3.scale.linear()
        .domain([0,max])
        .range([0,containerWidth])
        .clamp(true);
}