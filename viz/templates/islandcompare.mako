<html>
<head>
    <script src="/plugins/visualizations/islandcompare/static/js/d3.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/PapaParse/4.6.3/papaparse.min.js"></script>
    <script src="/plugins/visualizations/islandcompare/static/js/d3.phylogram.js"></script>
    <script src="/plugins/visualizations/islandcompare/static/js/newick.js"></script>
    <script src="/plugins/visualizations/islandcompare/static/js/genomecomparevis.js"></script>
    <link rel="stylesheet" href="/plugins/visualizations/islandcompare/static/css/linearplot.css" />
    <!-- Tell the browser to be responsive to screen width -->
    <meta content="width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no" name="viewport">
    <!-- Bootstrap 3.3.6 -->
    <link rel="stylesheet" href="/plugins/visualizations/islandcompare/static/bootstrap/css/bootstrap.min.css">
    <!-- Font Awesome -->
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.5.0/css/font-awesome.min.css">
    <!-- Ionicons -->
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/ionicons/2.0.1/css/ionicons.min.css">
    <!-- Theme style -->
    <link rel="stylesheet" href="/plugins/visualizations/islandcompare/static/css/AdminLTE.min.css">
    <!-- AdminLTE Skins -->
    <link rel="stylesheet" href="/plugins/visualizations/islandcompare/static/css/skins/skin-blue-light.min.css">
    <!-- Bootstrap Diag -->
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap3-dialog/1.34.7/css/bootstrap-dialog.min.css">
    <!-- jQuery 2.2.3 -->
    <script src="/plugins/visualizations/islandcompare/static/plugins/jQuery/jquery-2.2.3.min.js"></script>
    <!-- Bootstrap 3.3.6 -->
    <script src="/plugins/visualizations/islandcompare/static/bootstrap/js/bootstrap.min.js"></script>
    <!-- slimscroll -->
    <script src="/plugins/visualizations/islandcompare/static/plugins/slimScroll/jquery.slimscroll.min.js"></script>
    <!-- FastClick -->
    <script src="/plugins/visualizations/islandcompare/static/plugins/fastclick/fastclick.js"></script>
    <!-- AdminLTE App -->
    <script src="/plugins/visualizations/islandcompare/static/js/app.min.js"></script>
    <!-- Bootstrap Diag -->
    <script src="https://cdnjs.cloudflare.com/ajax/libs/bootstrap3-dialog/1.34.7/js/bootstrap-dialog.min.js"></script>

</head>
<body>
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
        //let Galaxy = getGalaxyInstance();
        //Galaxy.data.dialog(function(url){alert(url)})
        var container = new MultiVis("#visualization-body");

        var treeOrder;
        var parser = Papa.parse("/api/histories/" + hda.history.id + "/contents/"+ (new URLSearchParams(location.search)).get('dataset_id') + "/display", {
            download: true,
            delimiter: "\t",
            worker: false,
            skipEmptyLines: true,
            beforeFirstChunk: function(chunk) {
                // Add the newick data to generate a tree
                container.newickData = Newick.parse(chunk.match(/^##newick: ([^\n]+)/m)[1]);
                // Set the root of the tree
                container.newickRoot = container.newickData;
                var reg = /^##sequence-region ([^ ]+) [^ ]+ ([^ \n]+)/mg;
                var seq;
                while (seq = reg.exec(chunk)) {
                    container.backbone.addSequence(seq[1], seq[2], seq[1], seq[1]);
                }

                // Traverse tree to find sequence id order
                treeOrder = container.traverseTreeForOrder(container.newickRoot);

                return chunk.replace(/^#.*\n/mg, '')
            },
            step: function(stream, parser) {
                var row = stream.data[0];
                switch (row[2]) {
                    case 'genomic_island':
                        //Add GI
                        container.backbone.getSequences().find(function(seq){return seq.sequenceId == row[0];}).addGI(row[1], {start:row[3], end:row[4]})
                        break;
                    case 'match':
                        //Add alignment
                        var target = /Target=([^ ]+) ([^ ]+) ([^ ;\n]+)(?: ([\+\-\.]))?/.exec(row[8]);
                        if (row[6] == "-") {
                            row[3] *= -1;
                            row[4] *= -1;
                        }
                        if (target[4] == "-") {
                            target[2] *= -1;
                            target[3] *= -1;
                        }
                        container.backbone.addHomologousRegion(row[0], target[1], row[3], row[4], target[2], target[3]);
                        break;
                    case 'gene':
                        //Add AMR
                        container.backbone.getSequences().find(function(seq){return seq.sequenceId == row[0]}).addAMR([{start: row[4], end: row[5], strand:row[7]}]);
                        break;
                }
            },
            complete: function() {
                // at this scale, individual scaling for sequences may not be usable...so used fixed scale
                for (currentSeq in container.backbone.getSequences())
                    currentSeq.updateScale(0, container.getLargestSequenceSize(), container.getLargestSequenceSize());

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
</body>
</html>