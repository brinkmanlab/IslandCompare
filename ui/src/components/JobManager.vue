<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <div class="JobManager">
        <b-tabs v-model="current_tab">
            <b-tab title="Recent Jobs">
                <!-- TODO https://bootstrap-vue.js.org/docs/components/tabs#add-custom-content-to-tab-title -->
                <!-- Shows running and completed jobs -->
                <Jobs v-bind:workflow="workflow" ref="jobs">
                </Jobs>
            </b-tab>
            <b-tab class="help" title="Instructions">
                <p>Upload your data to be analysed by dragging and dropping it into the box to the right. Alternatively, click the "Upload datasets" button and select your datasets to upload.</p>
                <p><em>Permitted dataset formats are Genbank or EMBL.</em></p>
                <p>Once your datasets have been uploaded, select them by clicking in the box to the right. Hold Ctrl (⌘ for mac) to select multiple. Hold Shift to select a range.</p>
                <p>Now that you have selected your data to compare, make any necessary changes to the analysis parameters, and click submit.</p>
                <p>The pending job will appear below. Once complete a "Visualize" button will appear along with the option to download the analysis.</p>
                <p><em>Be sure to bookmark this page to return to your work. The above URL is unique to you.</em></p>
            </b-tab>
        </b-tabs>
        <!-- Displays history tagged with 'user_data' and invokes workflow from selected datasets -->
        <JobRunner v-bind:workflow="workflow"
                   v-bind:history="history"
                   v-bind:workflow_params="{
                      minimum_island_size: 5000,
                      minimum_homologous_region: 50,
                      minimum_cluster_size: 2,
                   }"
                   v-bind:selection_validator="selection_validator"
                   v-bind:upload_callback="upload_callback"
        >
            <template v-slot:workflow_params="params">
                <slot name="workflow_params" v-bind="params"></slot>
            </template>
        </JobRunner>
    </div>
</template>

<script>
    // This component is responsible for managing loading galaxies state and preparing the "user_data" history and loading the target workflow.
    // This role may need to be broken up in the future.
    import { getConfiguredWorkflow, getUploadHistory } from "@/app";
    import { galaxy_load } from "@/store";
    // TODO async load these components
    import JobRunner from './JobRunner.vue'
    import Jobs from './Jobs'

    let galaxy = null;
    galaxy_load.then(module=>galaxy = module); // This should always happen before anything uses it. (hopefully)

    export default {
        name: "JobManager",
        components: {
            JobRunner,
            Jobs,
        },
        props: {
            permitted_file_extensions: {
                type: Array,
                default: ()=>[],
            },
            selection_validator: {
                type: Function,
                default: selection => selection.length === 0 ? "Invalid dataset selection" : null,
            }
        },
        data() { return {
            current_tab: 1,
        }},
        methods: {
            upload_callback(file) {
                //This should never get called before loading the history
                let ext = file.name.match(/[^.]+$/);
                if (this.permitted_file_extensions.length === 0 || (ext && this.permitted_file_extensions.includes(ext[0]))) return file;
                let tmp_id = file.name+Math.floor(Math.random()*10**16).toString();
                galaxy.history_contents.HistoryDatasetAssociation.insert({data: {id: tmp_id, file: file, name: "Incorrect file format: " + file.name, hid: -1, history_id: this.history.id}});
                return null;
            },
        },
        asyncComputed: {
            workflow: getConfiguredWorkflow,
            history: getUploadHistory,
        },
        watch: {
            workflow(workflow, oldVal) {
                // Switch current tab to recent jobs if there are jobs
                if (oldVal === null
                    && workflow
                    && galaxy.workflows.WorkflowInvocation.query().has('history').with('history', q=>q.where('deleted', false)).where('workflow_id', workflow.id).count() > 0
                ) {
                    this.current_tab = 0;
                }
            },
        }
    }
</script>

<style scoped>
    .JobManager {
        display: grid;
        grid-template-areas:
                "help jobrunner";
        grid-template-columns: minmax(10em, 30em) 1fr;
        grid-template-rows: minmax(30em, auto) minmax(30em, auto);
        grid-gap: 1em;
    }

    .JobRunner {
        justify-self: stretch;
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
    }

    .help {
        max-width: 30em;
        padding: 1em;
    }

    .help em {
        color: red;
    }
</style>