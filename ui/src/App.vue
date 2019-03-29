<template>
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
        <Jobs v-if="workflow" v-bind:workflows="[workflow]" ref="jobs"/>
        <div ref="toast" class="toast"></div>
    </div>
</template>

<script>
    const workflow_id = "1cd8e2f6b131e891";
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
            history: null,
            workflow: null,
        }},
        methods: {
            toast(e) {
                //TODO
                console.log(e); // eslint-disable-line no-console
            },
        },
        mounted() {
            //TODO wait for galaxy server
            const self = this;
            galaxy.histories.History.$fetch().then(() => {
                self.history = galaxy.histories.History.query().where('tags', tags=>tags.includes('user_data')).withAll().first();
                if (!self.history) {
                    galaxy.histories.History.$create({
                        data: {
                            name: "Uploaded data", //TODO replace with app name
                        }
                    }).then(response => {
                        self.history = galaxy.histories.History.find(response.id);
                        self.history.tags.push('user_data');
                        self.history.upload();
                    });
                }

                //Load workflow
                galaxy.workflows.StoredWorkflow.$get({params: {id: workflow_id}}).then(response => {
                    galaxy.workflows.WorkflowInvocation.$fetch({params: {url: response.url}}).then(response=> {
                        self.workflow = galaxy.workflows.StoredWorkflow.find(workflow_id);
                        if (!self.workflow) {
                            throw "Cant find workflow id: " + workflow_id;
                        }
                        //TODO bandaid to fix incorrect workflow_id
                        self.workflow.invocations = galaxy.workflows.WorkflowInvocation.query().whereIdIn(response.map(a=>a.id)).with('history').get();
                        for (let i of self.workflow.invocations) i.workflow = self.workflow;
                    }); //TODO catch
                }); //TODO catch
            }); //TODO catch
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
