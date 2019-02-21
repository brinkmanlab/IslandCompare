<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/5.7.0/d3.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/PapaParse/4.6.3/papaparse.min.js"></script>
<script src="/plugins/visualizations/islandcompare/static/js/d3.phylogram.js"></script>
<script src="/plugins/visualizations/islandcompare/static/js/newick.js"></script>
<script src="/plugins/visualizations/islandcompare/static/js/genomecomparevis.js"></script>
<link rel="stylesheet" href="/plugins/visualizations/islandcompare/static/css/linearplot.css" />

<section class="content">
    <div class="row">
        <div class="col-xs-12">
            <div class="box">
                <div class="box-header with-border">
                    <div class="controlRow">
                        <div class="controlBody">
                            <span id="#controlTitle">Controls: </span>
                            <button id="resetTree">Reset Tree</button>
                            <button id="resetRange">Reset Zoom</button>
                            <button id="toggleBranchLength">Toggle Branches</button>
                            <button id="printerColors">Toggle Grayscale</button>
                        </div>
                        <div class="legendBody">
                            <span id="legendTitle">Legend: </span>
                            <svg id="genesLegend" height="10" width="10">
                                <circle cx="5" cy="5" r="4" stroke-width="1"></circle>
                            </svg>
                            <span class="legendText">Genes</span>
                            <svg id="homologosRegionLegend" height="10" width="10">
                                <circle cx="5" cy="5" r="4" stroke-width="1"></circle>
                            </svg>
                            <span class="legendText">Homologous Region</span>
                            <svg id="genomeLegend" height="10" width="10">
                                <circle cx="5" cy="5" r="4" stroke-width="1"></circle>
                            </svg>
                            <span class="legendText">Genome</span>
                            <svg id="amrLegend" height="10" width="10">
                                <polygon points="3,0,7,0,10,10,0,10"><title>AntiMicrobial Resistance Genes</title></polygon>
                            </svg>
                            <span title="AntiMicrobial Resistance Genes" class="legendText">AMR Genes</span>

                            <input type="checkbox" name="sigi" value="sigi" id="sigi-box" checked>Sigi-HMM
                            <input type="checkbox" name="islandpath" value="islandpath" id="islandpath-box" checked>IslandPath
                        </div>
                        <!--
                        --Currently Not used
                        <div class="controlSide">
                            <button id="toggleControls">(<<)</button>
                        </div>
                        -->
                    </div>
                </div>
                <!-- /.box-header -->
                <div class="box-body">
                    <div id="visualization-body"></div>
                </div>
                <!-- /.box-body -->
            </div>
            <!-- /.box -->
        </div>
        <!-- /.col -->
    </div>
    <!-- /.row -->
</section>
<!-- /.content -->

<script>
    $(document).ready(function() {
        let Galaxy = getGalaxyInstance();
        //Galaxy.data.dialog(function(url){alert(url)})
        var container = new MultiVis("#visualization-body");
        var parser = Papa.parse((new URLSearchParams(location.search)).get('dataset_id'), {
            download: true,
            delimiter: "\t",
            worker: true,
            skipEmptyLines: true,
            beforeFirstChunk: function(chunk) {
                // Add the newick data to generate a tree
                container.newickData = Newick.parse(chunk.match(/^##newick: ([^\n]+)/)[1]);
                // Set the root of the tree
                container.newickRoot = container.newickData;
                var reg = /^##sequence-region ([^ ]+) [^ ]+ ([^ ]+)/;
                var seq;
                while (seq = reg.exec(chunk)) {
                    container.backbone.addSequence(seq[1], seq[2], seq[1], seq[1]);
                }
                return chunk.replace(/^#.*\n/g, '')
            },
            step: function(row) {
                switch (row[2]) {
                    case 'genomic_island':
                        //Add GI
                        container.backbone.getSequences().search(function(seq){seq.sequenceId == row[1]}).addGI(row[2], )
                        break;
                    case 'match':
                        //Add alignment
                        target = /Parent=([^ ]+) ([^ ]+) ([^ ;\n]+)(?: ([\+\-\.]))?/.exec(row[9]);
                        container.backbone.addHomologousRegion(row[1], target[1], row[4], row[5], target[2], target[3]); //TODO strand
                        break;
                    case 'gene':
                        //Add AMR
                        container.backbone.getSequences().search(function(seq){seq.sequenceId == row[1]}).addAMR([{start: row[4], end: row[5], strand:row[7]}]);
                        break;
                }
            },
            complete: function() {
            // at this scale, individual scaling for sequences may not be usable...so used fixed scale
            currentSeq.updateScale(0, container.getLargestSequenceSize(), container.getLargestSequenceSize());

            // Add homolous regions to Visualization
            // TODO Fix ids for homologous region, doesnt use genome ids but rather goes from 0 -> N
            for (var i = 0; i < container.newickData["branchset"].length - 1; i++){
            for (var j = 0; j < data["alignment"][i].length; j++){
                container.backbone.addHomologousRegion(
                    i,
                    i + 1,
                    data["alignment"][i][j][0],
                    data["alignment"][i][j][1],
                    data["alignment"][i][j][2],
                    data["alignment"][i][j][3]
                )
            }
        }

        // render the graph
        container.render();
            }
        });


        $.ajax({
            url: url,
            beforeSend: function (request) {
                request.setRequestHeader("Authorization", 'Token'+' '+getAuthenticationCookie());
            },
            success: function(data) {
                // Add the newick data to generate a tree
                container.newickData = Newick.parse(data["newick"]);
                // Set the root of the tree
                container.newickRoot = container.newickData;

                // Add Genomes to Visualization
                for (var genomeIndex in container.newickData["branchset"]){
                    var genomeId = container.newickData["branchset"][genomeIndex]["name"];
                    var currentSeq = container.backbone.addSequence(
                        genomeId,
                        data["genomes"][genomeId]["length"],
                        data["genomes"][genomeId]["name"],
                        data["genomes"][genomeId]["name"]
                    );

                    currentSeq.addAMR(data["genomes"][genomeId]["amr_genes"]);

                    if ("user" in data["genomes"][genomeId]["genomic_islands"]){
                        currentSeq.addGI("user", data["genomes"][genomeId]["genomic_islands"]["user"])
                    }
                    else {
                        currentSeq.addGI("islandpath", data["genomes"][genomeId]["genomic_islands"]["islandpath"]);
                        currentSeq.addGI("sigi", data["genomes"][genomeId]["genomic_islands"]["sigi"]);
                        currentSeq.addGI("merged", data["genomes"][genomeId]["genomic_islands"]["merged"]);
                    }

                    // at this scale, individual scaling for sequences may not be usable...so used fixed scale
                    currentSeq.updateScale(0, container.getLargestSequenceSize(), container.getLargestSequenceSize());
                }

                // Add homolous regions to Visualization
                // TODO Fix ids for homologous region, doesnt use genome ids but rather goes from 0 -> N
                for (var i = 0; i < container.newickData["branchset"].length - 1; i++){
                    for (var j = 0; j < data["alignment"][i].length; j++){
                        container.backbone.addHomologousRegion(
                            i,
                            i + 1,
                            data["alignment"][i][j][0],
                            data["alignment"][i][j][1],
                            data["alignment"][i][j][2],
                            data["alignment"][i][j][3]
                        )
                    }
                }

                // render the graph
                container.render();
            }
        });

        var isPrinterColors = false;

        // Add listeners to button controls
        //Button resets any tree zooming that was done
        $("#resetTree").on("click", this, function () {
            container.resetAndRenderGraph();
        });
        //Button resets the range of the linear plot to (0,maxSize)
        $("#resetRange").on("click", this, function () {
            container.resetAndRenderRange();
        });
        //Button toggles if show true branch lengths or not
        $("#toggleBranchLength").on("click", this, function () {
            container.toggleTrueBranchLengths();
            container.transition();
        });
        //Button toggles printer colors on the linear plot
        $("#printerColors").on("click", this, function () {
            container.togglePrinterColors();
            container.render();
            // Update the Legend Colors
            isPrinterColors = !isPrinterColors;
            if (isPrinterColors) {
                $("#genomicIslandLegend").attr("class", "#genomicIslandsLegend print");
                $("#genesLegend").attr("class", "#genesLegend print");
                $("#genomeLegend").attr("class", "#genomeLegend print");
                $("#amrLegend").attr("class", "#amrLegend print");
            }
            else {
                $("#genomicIslandLegend").attr("class", "#genomicIslandsLegend");
                $("#genesLegend").attr("class", "#genesLegend");
                $("#genomeLegend").attr("class", "#genomeLegend");
                $("#amrLegend").attr("class", "#amrLegend");
            }
        });

        $("#sigi-box").change(function() {
            $(".sigi").toggle();
        });

        $("#islandpath-box").change(function() {
            $(".islandpath").toggle();
        });
    });
</script>