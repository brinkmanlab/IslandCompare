<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <div class="JobManager">
        <section class="help">
            <p>Upload your data to be analysed by dragging and dropping it into the box to the right. Alternatively, click the "Upload datasets" button and select your datasets to upload.</p>
            <p><em>Permitted dataset formats are Genbank or EMBL.</em></p>
            <p>Once your datasets have been uploaded, select them by clicking in the box to the right. Hold Ctrl (⌘ for mac) to select multiple. Hold Shift to select a range.</p>
            <p>Now that you have selected your data to compare, make any necessary changes to the analysis parameters, and click submit.</p>
            <p>The pending job will appear below. Once complete a "Visualize" button will appear along with the option to download the analysis.</p>
            <p><em>Be sure to bookmark this page to return to your work. The above URL is unique to you.</em></p>
        </section>
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
        <!-- Shows running and completed jobs -->
        <Jobs v-bind:workflow="workflow" ref="jobs">
            <template v-slot:header="">
                <h2>History</h2>
            </template>
            <template v-slot:functions="slot">
                <slot name="invocation_functions" v-bind="slot"></slot>
            </template>
        </Jobs>
    </div>
</template>

<script>
    // This component is responsible for managing loading galaxies state and preparing the "user_data" history and loading the target workflow.
    // This role may need to be broken up in the future.
    import { galaxy_load } from "@/store";

    // TODO async load these components
    import JobRunner from './JobRunner.vue'
    import Jobs from './Jobs'
    let galaxy = null;
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
            workflow_name: "IslandCompare",
            buildORM: galaxy_load.then(module=>{
                //Lazy load galaxy ORM as it is BIG
                galaxy = module;
                galaxy.register(this.$store);
                this.fetchedHistories = galaxy.histories.History.$fetch();
                this.fetchedWorkflows = galaxy.workflows.StoredWorkflow.$fetch({query: {show_published: true}});
            }),
            fetchedHistories: null,
            fetchedWorkflows: null,
            fetchedInvocations: null,
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
            async workflow() {
                // Load the workflow and all its components
                let galaxy = await galaxy_load;
                await this.buildORM;
                await this.fetchedWorkflows;
                let workflow = galaxy.workflows.StoredWorkflow.query().where('name', this.workflow_name).first();
                if (!workflow) {
                    let err = "IslandCompare workflow could not be found";
                    throw err;
                }
                if (!this.fetchedInvocations) this.fetchedInvocations = galaxy.workflows.WorkflowInvocation.$fetch({params: {url: workflow.url}, query: {view: "element", step_details: false}});
                await this.fetchedInvocations;
                await this.fetchedHistories;

                let result = galaxy.workflows.StoredWorkflow.query().with('invocations', invocations => {
                    invocations.whereHas('history', history => {
                        history.where('deleted', false).where('tags', tags => tags.includes(workflow.id));
                    }).with('workflow').with('history');
                }).find(workflow.id);
                return result || workflow; //TODO clean up logic above to allow for empty workflows
            },
            async history() {
                // Load the user_data history and all its datasets
                let galaxy = await galaxy_load;
                await this.buildORM;
                await this.fetchedHistories;
                let history = galaxy.histories.History.query().where('tags', tags=>tags.includes('user_data')).first();
                if (!history) {
                    let response = await galaxy.histories.History.$create({
                        data: {
                            name: "Uploaded data", //TODO replace with app name
                        }
                    });
                    history = galaxy.histories.History.find(response.id);
                    history.tags.push('user_data');
                    history.upload();
                } else {
                    await galaxy.history_contents.HistoryDatasetAssociation.$fetch({
                        params: {
                            url: history.contents_url,
                        }
                    });
                    await galaxy.history_contents.HistoryDatasetCollectionAssociation.$fetch({
                        params: {
                            url: history.contents_url,
                        }
                    });
                }
                return history;
            }
        },
    }
</script>

<style scoped>
    .JobManager {
        display: grid;
        grid-template-areas:
                "help jobrunner"
                "jobs jobs";
        grid-template-columns: minmax(10em, 30em) 1fr;
        grid-template-rows: minmax(30em, auto) minmax(30em, auto);
        grid-gap: 1em;
    }

    .JobManager .JobRunner {
        justify-self: stretch;
    }

    .JobManager .Jobs {
        grid-area: jobs;
        border-top: solid 2px gray;
        margin-top: 1em;
        padding-top: 1em;
    }

    .JobManager .help {
        max-width: 30em;
        padding: 1em;
    }

    .JobManager .help em {
        color: red;
    }
</style>