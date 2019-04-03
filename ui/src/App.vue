<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <div id="app">
        <section class="help">
            <p>Upload your data to be analysed by dragging and dropping it into the box to the right. Alternatively, click the "Upload datasets" button and select your datasets to upload.</p>
            <p><em>Permitted dataset formats are Genbank or EMBL.</em></p>
            <p>Once your datasets have been uploaded, select them by clicking in the box to the right. Hold Ctrl to select multiple. Hold Shift to select a range.</p>
            <p>Now that you have selected your data to compare, make any necessary changes to the analysis parameters, and click submit.</p>
            <p>The pending job will appear below. Once complete a "Visualize" button will appear along with the option to download the analysis.</p>
        </section>
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
                   @workflow-invoked="$refs.jobs.update()"
        >
            <template v-slot:workflow_params="params">
                <label>Minimum island size<input type="number" min="0" required v-model.number.lazy="params.minimum_island_size" /></label>
                <label>Minimum homologous region<input type="number" min="0" required v-model.number.lazy="params.minimum_homologous_region" /></label>
                <label>Minimum cluster size<input type="number" min="1" required v-model.number.lazy="params.minimum_cluster_size" /></label>
            </template>
        </JobRunner>
        <Jobs v-if="workflow" v-bind:workflows="[workflow]" ref="jobs">
            <template v-slot:functions="slot" >
                <a v-if="slot.model.history.state === 'ok' && slot.model.outputs['IslandCompare Result']"
                   v-bind:href="'/plugins/visualizations/islandcompare/show?key=admin&dataset_id='+slot.model.outputs['IslandCompare Result'].id">
                    Visualize
                </a>
            </template>
        </Jobs>
        <div ref="toast" class="toast"></div>
    </div>
</template>

<script>
    import * as galaxy from '@/galaxy'
    import JobRunner from './components/JobRunner.vue'
    import Jobs from './components/Jobs'
    export default {
        name: 'app',
        components: {
            JobRunner,
            Jobs,
        },
        data() { return {
            workflow_id: "1cd8e2f6b131e891",
            fetchedHistories: galaxy.histories.History.$fetch(),
        }},
        methods: {
            toast(e) {
                //TODO
                console.log(e); // eslint-disable-line no-console
            },
            async fetchHistories() {

            }
        },
        asyncComputed: {
            async workflow(){
                let response = await galaxy.workflows.StoredWorkflow.$get({params: {id: this.workflow_id}}); //eslint-disable-line
                response = await galaxy.workflows.WorkflowInvocation.$fetch({params: {url: response.url}, /*query: {view: "element"}*/});
                await this.fetchedHistories;

                //await Promise.all(missing.map(a=>galaxy.histories.History.$get({params: {id: a}})));
                return galaxy.workflows.StoredWorkflow.query().with('invocations', invocations=>{
                    invocations.whereHas('history', history=>{
                        history.where('deleted', false)
                    }).with('workflow').with('history')
                }).find(this.workflow_id);
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
                    history = galaxy.histories.History.query().with('datasets.history').find(history.id);
                }
                return history;
            }
        },
        mounted() {
            //TODO wait for galaxy server
        },
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
        margin: 10em;
    }

    #app .JobRunner {
        justify-self: stretch;
    }

    #app .Jobs {
        grid-area: jobs;
    }

    #app .help {
        max-width: 30em;
    }
</style>
