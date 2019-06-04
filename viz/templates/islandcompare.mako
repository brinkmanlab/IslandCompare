<html>
<head>
    <script src="static/js/d3.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/PapaParse/4.6.3/papaparse.min.js"></script>
    <script src="static/js/d3.phylogram.js"></script>
    <script src="static/js/newick.js"></script>
    <script src="static/js/genomecomparevis.js"></script>

    <!-- Tell the browser to be responsive to screen width -->
    <meta content="width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no" name="viewport">
    <!-- Bootstrap 3.3.6 -->
    <link rel="stylesheet" href="static/bootstrap/css/bootstrap.min.css">
    <!-- Font Awesome -->
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.5.0/css/font-awesome.min.css">
    <!-- Ionicons -->
    <!--link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/ionicons/2.0.1/css/ionicons.min.css"-->
    <!-- Theme style -->
    <!--link rel="stylesheet" href="static/css/AdminLTE.min.css"-->
    <!-- AdminLTE Skins -->
    <!--link rel="stylesheet" href="static/css/skins/skin-blue-light.min.css"-->
    <!-- Bootstrap Diag -->
    <!--link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap3-dialog/1.34.7/css/bootstrap-dialog.min.css"-->
    <!-- jQuery 2.2.3 -->
    <script src="static/plugins/jQuery/jquery-2.2.3.min.js"></script>
    <!-- Bootstrap 3.3.6 -->
    <script src="static/bootstrap/js/bootstrap.min.js"></script>
    <!-- slimscroll -->
    <!--script src="static/plugins/slimScroll/jquery.slimscroll.min.js"></script-->
    <!-- FastClick -->
    <!--script src="static/plugins/fastclick/fastclick.js"></script-->
    <!-- AdminLTE App -->
    <!--script src="static/js/app.min.js"></script-->
    <!-- Bootstrap Diag -->
    <script src="https://cdnjs.cloudflare.com/ajax/libs/bootstrap3-dialog/1.34.7/js/bootstrap-dialog.min.js"></script>
    <link rel="stylesheet" href="static/css/linearplot.css" />
</head>
<body>
<section class="content">
    <div class="box">
        <div class="box-header with-border">
            <div class="controlRow row">
                <div class="col-xs-5">
                    <div id="btn-div" class="clearfix btn-group btn-group-sm">
                        <button id="toggleBranchLength" class="btn btn-primary btn-flat">Toggle Branches</button>
                        <button id="resetTree" class="btn btn-primary btn-flat" title="Resets the tree to show all clades. Use reset zoom to also reset to the full view.">Reset Tree</button>
                        <button id="resetGIs" class="btn btn-primary btn-flat clusterButton" style="display: none;">
                            Reset GIs
                        </button>
                        <div class="btn-group btn-group-sm">
                            <button class="btn btn-primary btn-flat dropdown-toggle" data-toggle="dropdown">
                                Save image <span class="caret"></span>
                            </button>
                            <ul class="dropdown-menu">
                                <li><a class="dropdown-item" id="saveSvg" href="#">svg</a></li>
                                <li><a class="dropdown-item" id="savePng" href="#">png</a></li>
                            </ul>
                        </div>
                        <button id="printerColors" class="btn btn-primary btn-flat" style="margin-left: 0px;">Toggle Colour</button>
                        <button id="resetRange" class="btn btn-primary btn-flat">Reset Zoom</button>
                        <button id="zoomCluster" class="btn btn-primary btn-flat clusterButton" style="display: none;">
                            Zoom to Cluster
                        </button>
                    </div>
                </div>
                <div class="col-xs-1" id="GITitle">
                    <span class="clearfix">Colour GIs by: </span>
                </div>
                <div class="col-xs-1 GIColour">
                    <div class="GIColourBody">
                        <div class="clearfix">
                            <label>
                                <input type="radio" name="GIColour" value="similarity" checked/>Similarity
                            </label>
                            <label id="predictorButton">
                                <input type="radio" name="GIColour" value="predictor" />Predictor
                            </label>
                        </div>
                        <div class="clearfix">
                            <span id="sigiToggle">
                                <label>
                                    <input type="checkbox" name="sigi" value="sigi" id="sigi-box" checked disabled>Sigi-HMM
                                </label>
                            </span>
                            <span id="islandpathToggle">
                                <label>
                                    <input type="checkbox" name="islandpath" value="islandpath" id="islandpath-box" checked disabled>IslandPath
                                </label>
                            </span>
                        </div>
                    </div>
                </div>
                <div class="col-xs-1" id="legendTitle">
                    <span class="clearfix">Legend: </span>
                </div>
                <div class="col-xs-2 legend">
                    <div class="legendBody clearfix">
                        <span class="legendText">
                            <svg id="homologousRegionLegend" height="10" width="10">
                                <circle cx="5" cy="5" r="4" stroke-width="1"></circle>
                            </svg>
                            Homologous Region
                        </span>

                        <span title="AntiMicrobial Resistance Gene" class="legendText">
                            <svg id="amrLegend" height="10" width="10">
                                <circle cx="5" cy="5" r="4" stroke-width="1"><title>AntiMicrobial Resistance Gene</title></circle>
                            </svg>
                            AMR Gene
                        </span>
                        <span class="loadingGenes" hidden>
                            <i class='fa fa-spinner fa-spin '></i> Loading Genes
                        </span>
                        <span hidden class="genesLegend-toggle">
                            <span class="legendText">
                                <svg id="genesLegend" height="10" width="10">
                                    <circle cx="5" cy="5" r="4" stroke-width="1"></circle>
                                </svg>
                                Gene
                            </span>

                            <span class="legendText">
                                <svg id="tRNALegend" height="10" width="10">
                                    <circle cx="5" cy="5" r="4" stroke-width="1"></circle>
                                </svg>
                                tRNA
                            </span>

                            <span class="legendText">
                                <svg id="rRNALegend" height="10" width="10">
                                <circle cx="5" cy="5" r="4" stroke-width="1"></circle>
                                </svg>
                                rRNA
                            </span>
                        </span>
                    </div>
                    <div class="predictorLegend clearfix" hidden>
                        <span id="sigiLegend" class="legendText">
                            <svg height="10" width="10">
                                <circle cx="5" cy="5" r="4" stroke-width="1"></circle>
                            </svg>
                            Sigi-HMM Genomic Island
                        </span>
                        <span id="islandpathLegend" class="legendText">
                            <svg height="10" width="10">
                                <circle cx="5" cy="5" r="4" stroke-width="1"></circle>
                            </svg>
                            IslandPath Genomic Island
                        </span>

                    </div>
                    <div id="clusterLegend" class="clearfix" hidden>
                        <svg height="10" width="10">
                            <circle cx="5" cy="5" r="4" stroke-width="1"></circle>
                        </svg>
                        <span id="clusterLegendSpan" class="legendText"></span>
                    </div>
                </div>
                <!--
                --Currently Not used
                <div class="controlSide">
                    <button id="toggleControls">(<<)</button>
                </div>
                -->
                <div id="tips" class="col-xs-2"><span class="clearfix"><span style="font-weight: bold;">Tips:</span> Select a node of a tree to show only that clade. Click and drag over a genomic region to zoom in. Select a genomic island to view its cluster (note that the default is islands are coloured by cluster). Hover over genes to see the annotation.</span></div>
            </div>
        </div>
        <!-- /.box-header -->
        <div class="box-body">
            <div id="visualization-body">
                <h1 id="loading" class="text-center"><i class='fa fa-spinner fa-spin'></i> Loading results</h1>
            </div>
            <svg>
                <filter id="shadow" x="-5" y="-5" width="200" height="200">
                    <feGaussianBlur in="SourceAlpha" stdDeviation="5" />
                    <feColorMatrix result="matrixOut" in="offOut" type="matrix"
                                   values="0   0   0   0   0
                                0   0   0   0   0
                                0   0   0   0   0
                                0   0   0   2   0" />
                    <feMerge>
                        <feMergeNode/>
                        <feMergeNode in="SourceGraphic"/>
                    </feMerge>
                </filter>
            </svg>
        </div>
        <!-- /.box-body -->
    </div>
    <!-- /.box -->
</section>

<script>
    $(document).ready(function() {
        //let Galaxy = getGalaxyInstance();
        //Galaxy.data.dialog(function(url){alert(url)})
        var container = new MultiVis("#visualization-body");

        var treeOrder;
        var parser = Papa.parse("/galaxy/datasets/" + (new URLSearchParams(location.search)).get('dataset_id') + "/display", {
            download: true,
            delimiter: "\t",
            worker: false,
            skipEmptyLines: true,
            beforeFirstChunk: function(chunk) {
                // Add the newick data to generate a tree
                let nwk = chunk.match(/^##newick: ([^\n]+)/m);
                if (nwk && nwk.length > 1) {
                    container.newickData = Newick.parse(nwk[1]);
                } else {
                    console.log("No newick found in first chunk.", chunk);
                }
                // Set the root of the tree
                container.newickRoot = container.newickData;
                var reg = /^##sequence-region ([^ ]+) [^ ]+ ([^ \n]+)/mg;
                var seq;
                while (seq = reg.exec(chunk)) {
                    container.backbone.addSequence(seq[1], Number(seq[2]), Number(seq[1]), seq[1]);
                }

                // Traverse tree to find sequence id order
                treeOrder = container.traverseTreeForOrder(container.newickRoot);

                container.sequenceOrder = container.backbone.getIndicesFromIds(treeOrder);

                return chunk.replace(/^#.*\n/mg, '');
            },
            step: function(stream, parser) {
                var row = stream.data[0];
                row[3] = Number(row[3]);
                row[4] = Number(row[4]);
                switch (row[2]) {
                    case 'genomic_island':
                        //Add GI
                        var program;
                        switch(row[1]) {
                            case 'Colombo/SigiHMM':
                                program = 'sigi';
                                break;
                            case 'islandpath':
                                program = 'islandpath';
                                break;
                            default:
                                program = 'merged';
                                break;
                        }
                        let cluster = /cluster=([^;\n]+)/.exec(row[8]);
                        cluster = cluster ? cluster[1] : null;
                        let color = /color=([^;\n]+)/.exec(row[8]);
                        color = color ? color[1] : null;
                        let parent = /Parent=([^;\n]+)/.exec(row[8]);
                        parent = parent ? parent[1] : null;
                        container.backbone.getSequences().find(function(seq){return seq.sequenceId === row[0];}).addGI(program, {start:row[3], end:row[4], cluster: cluster, color: color, parent: parent});
                        break;
                    case 'match':
                        //Add alignment
                        let target = /Target=([^ ]+) ([^ ]+) ([^ ;\n]+)(?: ([\+\-\.]))?/.exec(row[8]);
                        target[2] = Number(target[2]);
                        target[3] = Number(target[3]);
                        let dist = treeOrder.indexOf(target[1]) - treeOrder.indexOf(row[0]);
                        if (dist == 1 || dist == -1) {
                            if (row[6] == "-") {
                                row[3] *= -1;
                                row[4] *= -1;
                            }
                            if (target[4] == "-") {
                                target[2] *= -1;
                                target[3] *= -1;
                            }
                            if (dist > 0)
                                container.backbone.addHomologousRegion(row[0], row[3], row[4], target[2], target[3]);
                            else
                                container.backbone.addHomologousRegion(target[1], target[2], target[3], row[3], row[4]);
                        }
                        break;
                    case 'gene':
                    case 'pseudogene':
                        if (row[1] == "RGI-CARD") {
                            //Add AMR
                            container.backbone.getSequences().find(function (seq) {
                                return seq.sequenceId == row[0]
                            }).addAMR([{start: row[3], end: row[4], strand: row[6]}]);
                        } else {
                            //Add gene
                            var name = /Name=([^;\n]+)/.exec(row[8]);
                            name = name ? name[1] : "";
                            var gene_type = /Type=([^;\n]+)/.exec(row[8]);
                            gene_type = gene_type ? gene_type[1] : "";
                            var locus_tag = /locus_tag=([^;\n]+)/.exec(row[8]);
                            locus_tag = locus_tag ? locus_tag[1] : "";
                            var product = /product=([^;\n]+)/.exec(row[8]);
                            product = product ? product[1] : "";
                            container.backbone.getSequences().find(function (seq) {
                                return seq.sequenceId == row[0]
                            }).addGene({start: row[3], end: row[4], strand: row[6], gene: name, locus_tag: locus_tag, product: product, type: gene_type})
                        }
                        break;
                }
            },
            complete: function() {
                // at this scale, individual scaling for sequences may not be usable...so used fixed scale
                for (currentSeq of container.backbone.getSequences())
                    currentSeq.updateScale(0, container.getLargestSequenceSize(), container.getLargestSequenceSize());

                $("#loading").hide();

                // render the graph
                container.render();
                // Hook resize event to redraw
                //this.container.on('resize', function(){
                //    self.transition();
                //});
                $(window).on('resize', function(){
                    container.resetAndRenderRange();
                });
            }
        });

        var isPrinterColors = false;

        // Add listeners to button controls
        $(document).ready(function() {
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
                    $("#tRNALegend").attr("class", "#tRNALegend print");
                    $("#rRNALegend").attr("class", "#rRNALegend print");
                    $("#genomeLegend").attr("class", "#genomeLegend print");
                    $("#amrLegend").attr("class", "#amrLegend print");
                }
                else {
                    $("#genomicIslandLegend").attr("class", "#genomicIslandsLegend");
                    $("#genesLegend").attr("class", "#genesLegend");
                    $("#tRNALegend").attr("class", "#tRNALegend");
                    $("#rRNALegend").attr("class", "#rRNALegend");
                    $("#genomeLegend").attr("class", "#genomeLegend");
                    $("#amrLegend").attr("class", "#amrLegend");
                }
            });
            $("#saveSvg").on("click", this, function() {
                container.saveImage("svg");
            });
            $("#savePng").on("click", this, function() {
                container.saveImage("png");
            });
            $("#resetGIs").on("click", this, function() {
                container.toggleClusterView(container.cluster, null);
            });
            $("#zoomCluster").on("click", this, function() {
                container.zoomCluster();
            });
            // Changes to sigi/islandpath checkboxes will hide/show corresponding GIs
            $(".GIColourBody :checkbox").change(container.predictorToggle);
            // Changes to GI colour method will hide/show the merged GI set and enable/disable the predictor checkboxes
            $("input[name=GIColour]").change(container.GIColourToggle);
            // Set all buttons to have the same width
            var buttons = $("#btn-div").children();
            //buttons.width(buttons.width());
        });
    });
</script>
</body>
</html>