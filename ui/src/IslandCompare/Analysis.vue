<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <b-container class="Analysis" fluid>
        <b-row align-h="between">
            <b-col xl="5">
                <b-tabs v-model="current_tab">
                    <b-tab title="Recent Jobs">
                        <!-- TODO https://bootstrap-vue.js.org/docs/components/tabs#add-custom-content-to-tab-title -->
                        <!-- Shows running and completed jobs -->
                        <Jobs v-bind:invocationsPromise="getInvocations(workflowPromise)" ref="jobs" @click.native.prevent="$router.push('/history')">
                        </Jobs>
                    </b-tab>
                    <b-tab class="help" title="Instructions">
                        <p>Check out these <b-link :to="`visualize?src=${origin}/demo/listeria_sample_analysis.gff3`">example Listeria</b-link> or <b-link :to="`visualize?src=${origin}/demo/pseudomonas_sample_analysis.gff3`">example Pseudomonas</b-link> analyses.</p>
                        <p>Upload your data to be analysed by dragging and dropping it into the box to the right. Alternatively, click the upload button <i class="icon-file-upload"></i> and select your datasets to upload.</p>
                        <p><em>Permitted dataset formats are Genbank or EMBL.</em></p>
                        <p>Once your datasets have been uploaded, select them by clicking in the box to the right. Hold Ctrl (⌘ for mac) to select multiple. Hold Shift to select a range.</p>
                        <p>Optionally, you can upload a newick formatted tree which will decide the order of alignment in the visualization.</p>
                        <p>If you are analysing draft genomes, select a reference genome that is closest to your datasets species. The drafts will be stitched in the order of alignment to the reference. If no reference is selected the draft contigs will be stitched in the order they appear in the uploaded dataset. The dataset label will be used for the stitched sequence identifier. If you specify a newick file, it must include this identifier.</p>
                        <p>Now that you have selected your data to compare, specify a label to identify the analysis, and click submit.</p>
                        <p>The pending job will appear in the Recent Jobs tab. Once complete a "Visualize" button will appear in the Job History page along with the option to download the analysis.</p>
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

    import {getUUID, gidPromise} from "@/auth";

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
        methods: {
            getInvocations: getInvocations,
            updateUUID() {
                // Force uuid into url when navigating to this page
                const uuid = getUUID();
                if (!('uuid' in this.$route.query) && uuid) {
                    this.$router.replace({query: {uuid: uuid}});
                }
            }
        },
        computed: {
            jobCount() {
                if (this.workflow)
                    return WorkflowInvocation.query().has('history').with('history', q=>q.where('deleted', false).where('$isDeleting', false)).where('workflow_id', this.workflow.id).count();
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
                this.updateUUID();
            },
        },
        mounted() {
            this.updateUUID();
        },
    }
</script>

<style scoped>
    .help {
        /*max-width: 30em;*/
        padding: 1em;
        font-size: 0.8em;
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

    .Jobs >>> .galaxy-workflow-invocation .galaxy-history-label, .Jobs >>> .galaxy-workflow-invocation .galaxy-workflow-invocation-progress {
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