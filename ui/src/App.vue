<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <div id="app">
        <section class="help">
            <p>Upload your data to be analysed by dragging and dropping it into the box to the right. Alternatively, click the "Upload datasets" button and select your datasets to upload.</p>
            <p><em>Permitted dataset formats are Genbank or EMBL.</em></p>
            <p>Once your datasets have been uploaded, select them by clicking in the box to the right. Hold Ctrl to select multiple. Hold Shift to select a range.</p>
            <p>Now that you have selected your data to compare, make any necessary changes to the analysis parameters, and click submit.</p>
            <p>The pending job will appear below. Once complete a "Visualize" button will appear along with the option to download the analysis.</p>
            <p><em>Be sure to bookmark this page to return to your work. The above URL is unique to you.</em></p>
        </section>
        <!-- Displays history tagged with 'user_data' and invokes workflow from selected datasets -->
        <JobRunner v-if="workflow && history"
                   v-bind:workflow="workflow"
                   v-bind:history="history"
                   v-bind:workflow_params="{
                      minimum_island_size: 5000,
                      minimum_homologous_region: 50,
                      minimum_cluster_size: 2,
                   }"
                   v-bind:selection_validator="selection=>selection.length<2?'You must select more than one dataset for comparison':null"
                   @toast="toast"
        >
            <template v-slot:workflow_params="params">
                <label>Minimum island size<input type="number" min="0" required v-model.number.lazy="params.minimum_island_size" /></label>
                <label>Minimum homologous region<input type="number" min="0" required v-model.number.lazy="params.minimum_homologous_region" /></label>
                <label>Minimum cluster size<input type="number" min="1" required v-model.number.lazy="params.minimum_cluster_size" /></label>
            </template>
        </JobRunner>
        <!-- Shows running and completed jobs -->
        <Jobs v-if="workflow" v-bind:workflow="workflow" ref="jobs">
            <template v-slot:header="">
                <h2>History</h2>
            </template>
            <template v-slot:functions="slot">
                <template v-if="slot.history.state === 'ok' && slot.outputs['IslandCompare Result']">
                    <a v-bind:href="`/plugins/visualizations/islandcompare/show?dataset_id=${slot.outputs['IslandCompare Result'].id}` | auth">Visualize</a>
                    <a v-bind:href="`/datasets/${slot.outputs['IslandCompare Result'].id}/display?to_ext=gff3` | auth">Download</a>
                </template>
            </template>
        </Jobs>
        <Toast ref="toast" />
    </div>
</template>

<script>
    import * as galaxy from '@/galaxy'
    import JobRunner from './components/JobRunner.vue'
    import Jobs from './components/Jobs'
    import Toast from './components/Toast'
    export default {
        name: 'app',
        components: {
            JobRunner,
            Jobs,
            Toast,
        },
        data() { return {
            workflow_name: "IslandCompare",
            fetchedHistories: galaxy.histories.History.$fetch(),
        }},
        methods: {
            toast(e) {
                this.$refs['toast'].show(e);
                console.log(e); // eslint-disable-line no-console
            },
        },
        asyncComputed: {
            async workflow(){
                await galaxy.workflows.StoredWorkflow.$fetch({query: { show_published: true }});
                //await galaxy.workflows.StoredWorkflow.$get({params: {id: "f597429621d6eb2b"}}); //eslint-disable-line
                let workflow = galaxy.workflows.StoredWorkflow.query().where('name', this.workflow_name).first();
                console.log(workflow); //eslint-disable-line
                if (!workflow) {
                    let err = "IslandCompare workflow could not be found";
                    this.toast(err);
                    throw err;
                }
                await galaxy.workflows.WorkflowInvocation.$fetch({params: {url: workflow.url}, /*query: {view: "element"}*/});
                await this.fetchedHistories;

                let result = galaxy.workflows.StoredWorkflow.query().with('invocations', invocations=>{
                    invocations.whereHas('history', history=>{
                        history.where('deleted', false).where('tags', tags=>tags.includes(workflow.id));
                    }).with('workflow').with('history');
                }).find(workflow.id);
                return result || workflow; //TODO clean up logic above to allow for empty workflows
            },
            async history() {
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
        mounted() {
            //TODO wait for galaxy server
        },
        errorCaptured(err, vm, info) { //eslint-disable-line
            this.toast(err);
            return false;
        }
    }
</script>

<style>
    #app {
        font-family: 'Avenir', Helvetica, Arial, sans-serif;
        -webkit-font-smoothing: antialiased;
        -moz-osx-font-smoothing: grayscale;
        display: grid;
        grid-template-areas:
                "help jobrunner"
                "jobs jobs";
        grid-template-columns: minmax(10em, 30em) 1fr;
        grid-template-rows: minmax(30em, auto);
        grid-gap: 1em;
        margin: 10%;
    }

    #app .JobRunner {
        justify-self: stretch;
    }

    #app .Jobs {
        grid-area: jobs;
        border-top: solid 2px gray;
        margin-top: 1em;
        padding-top: 1em;
    }

    #app .help {
        max-width: 30em;
    }

    #app .help em {
        color: red;
    }

    #app .Toast {
        position: fixed;
        top: 20%;
        left: 50%;
        transform: translate(-50%, -50%);
    }
</style>
