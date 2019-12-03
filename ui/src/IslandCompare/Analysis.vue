<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <b-container class="Analysis" fluid>
        <b-row align-h="between">
            <b-col xl="5">
                <b-tabs v-model="current_tab" class="analysis-tabs">
                    <b-tab title="Recent Jobs">
                        <!-- TODO https://bootstrap-vue.js.org/docs/components/tabs#add-custom-content-to-tab-title -->
                        <!-- Shows running and completed jobs -->
                        <Jobs v-bind:invocationsPromise="getInvocations(workflowPromise)" ref="jobs" @click.native.prevent="$router.push('/history')">
                        </Jobs>
                    </b-tab>
                    <b-tab class="help" title="Instructions">
                        <p>Please see <b-link to="/about" >About</b-link> and <b-link to="/faq">FAQ</b-link> for more information.</p>
                        <p>Check out these <b-link :to="`visualize?src=${origin}/demo/listeria_sample_analysis.gff3`">example Listeria</b-link> or <b-link :to="`visualize?src=${origin}/demo/pseudomonas_sample_analysis.gff3`">example Pseudomonas</b-link> analyses.</p>

                        <p><b-button @click="start_tour('tour')" class="tutorial-start" variant="info" size="sm">Tutorial</b-button></p>

                        <ol>
                            <li>
                                <em>Upload your data</em>
                                <p>Drag and drop your files into the 'Input datasets' box on the right. Alternatively, click the upload button <i class="icon-file-upload"></i> and select your datasets to upload. <em>Permitted dataset formats are Genbank or EMBL.</em></p>
                            </li>
                            <li><em>Select your data by clicking them</em><p>Hold Ctrl (⌘ for mac) to select multiple. Hold Shift to select a range.</p></li>
                            <li><em>Upload a phylogenetic tree in newick format (Optional)</em><p>The tree will decide the order of alignment in the visualization. The file must have a '.newick' extension and refer to the genomes by accession. If no tree is provided, one will be automatically computed.</p></li>
                            <li><em>Select a reference genome (Optional)</em><p>If you are analysing draft genomes, select a reference genome that is closest to your datasets species. The drafts will be stitched in the order of alignment to the reference. If no reference is selected the draft contigs will be stitched in the order they appear in the uploaded dataset. The dataset label will be used for the stitched sequence identifier. If you specify a newick file, it must include this identifier.</p></li>
                            <li><em>Specify a label to identify the analysis, and click submit</em><p>The pending job will appear in the Recent Jobs tab. Once complete a "Visualize" button will appear in the Job History page along with the option to download the analysis.</p></li>
                        </ol>
                        <p><em>Be sure to bookmark this page to return to your work. The above URL is unique to you, and will be retained for three months following your last activity.</em></p>
                    </b-tab>
                </b-tabs>
            </b-col>
            <b-col xl="6">
                <JobRunner v-bind:historyPromise="historyPromise"
                           v-bind:workflowPromise="workflowPromise"
                           v-bind:selection_validator="selection=>selection.length<2?'You must select more than one dataset for comparison':null"
                           @galaxy-workflow-invocation="current_tab=0"
                />
            </b-col>
        </b-row>
    </b-container>
</template>

<script>
    import { getConfiguredWorkflow, getUploadHistory, getInvocations } from "@/app";
    // TODO async load these components
    import JobRunner from '@/components/JobRunner'
    import Jobs from '@/components/Jobs'

    import {updateRoute, gidPromise} from "@/auth";

    export default {
        name: "Analysis",
        components: {
            JobRunner,
            Jobs,
        },
        data() { return {
            current_tab: 1,
            workflowPromise: getConfiguredWorkflow(),
            historyPromise: getUploadHistory(),
            origin: window.location.origin,
        }},
        props: {
            tour: {
                type: String,
                default: '',
            }
        },
        methods: {
            getInvocations: getInvocations,
            start_tour(tour) {
                if (this.workflow && tour) this.$tours[tour].start();
            }
        },
        computed: {
            jobCount() {
                if (this.workflow)
                    return this.workflow.get_invocations_query().count();
                return 0;
            },
        },
        asyncComputed: {
            workflow() { return this.workflowPromise },
            //history() { return this.historyPromise },
            uuid() { return gidPromise },
        },
        watch: {
            jobCount(count) {
                // Switch current tab to recent jobs if there are jobs
                if (count > 0) {
                    this.current_tab = 0;
                }
            },
            uuid() {
                updateRoute(this.$router, this.$route);
            },
            workflow() {
                this.start_tour(this.tour);
            },
            tour() {
                if (this.tour) {
                    updateRoute(this.$router, this.$route);
                    if (this.workflow) this.$tours[this.tour].start();
                }
            }
        },
        activated() {
            updateRoute(this.$router, this.$route);
        },
        deactivated() {
            if (this.tour) this.$tours[this.tour].stop();
        }
    }
</script>

<style scoped>
    .analysis-tabs >>> .nav-tabs {
        font-size: 0.8em;
    }

    .analysis-tabs >>> .tab-content {
        border: 1px solid #dee2e6;
        border-top: 0;
        border-bottom-left-radius: 0.25rem;
        border-bottom-right-radius: 0.25rem;
        min-height: 70vh;
    }

    .JobRunner >>> .galaxy-workflow-parameter-dataset .b-table-sticky-header {
        min-height: 4rem;
    }

    .JobRunner >>> .galaxy-workflow-parameters .Input-datasets .b-table-sticky-header {
        height: 10rem;
    }

    .help {
        /*max-width: 30em;*/
        padding: 1em;
        font-size: 0.8em;
    }

    .help ol {
        padding-left: inherit;
        padding-top: 0;
    }

    .help em {
        font-weight: bold;
    }

    .Jobs {
        padding: 0.5vw;
        border-collapse: separate;
    }

    .Jobs >>> thead {
        display: none;
    }

    .Jobs >>> .galaxy-workflow-invocation > :not(.galaxy-history) {
        display: none;
        /*visibility: collapse;*/
    }

    .Jobs >>> .galaxy-workflow-invocation > .galaxy-history > * {
        display: none;
    }

    .Jobs >>> .galaxy-workflow-invocation .galaxy-history-label, .Jobs >>> .galaxy-workflow-invocation .galaxy-workflow-invocation-progress:not(.new), .Jobs >>> .galaxy-workflow-invocation-state.new {
        display: table-cell;
        font-size: 0.7em;
    }

    .Jobs >>> .galaxy-workflow-invocation .galaxy-workflow-invocation-progress {
        width: 100%;
    }

    .Jobs >>> .galaxy-workflow-invocation .galaxy-history-label {
        padding-left: 1em;
        padding-right: 1em;
        white-space: nowrap;
    }
</style>