<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <div class="JobHistory">
        <h1>Job History</h1>
        <Jobs v-bind:invocations="invocations">
            <template v-slot:functions="slot">
                <template v-if="slot.done && slot.model.outputs['Results']">
                    <b-link v-bind:to="`/visualize/${slot.model.outputs['Results'].id}` | auth">Visualize</b-link>
                    <WorkflowInvocationOutputDownload :outputs="slot.outputs" :url_xform="url_xform" />
                </template>
            </template>
        </Jobs>
    </div>
</template>

<script>
    import Jobs from "../components/Jobs";
    import {getConfiguredWorkflow, getInvocations} from "../app";
    import {updateRoute} from "../auth";
    import {fetchState} from "../app";
    import WorkflowInvocationOutputDownload from "galaxy-client/src/workflows/WorkflowInvocationOutputDownload";

    export default {
        name: "JobHistory",
        components: {WorkflowInvocationOutputDownload, Jobs },
        data() {return{
            auth_fail: false,
        }},
        methods: {
            init(force) {
                if (this.auth_fail || force) {
                    this.auth_fail = false;
                    fetchState().then(()=>{ // gidPromise can fail and then be reassigned by fetch. This happens late and is a race condition with below
                        updateRoute(this.$router, this.$route);
                    }).catch(() => { // TODO emit an event on $root and listen for it instead
                        this.auth_fail = true;
                    });
                }
            },
            url_xform(x) {
                return this.$options.filters.auth(this.$options.filters.galaxybase(x))
            }
        },
        computed: {
            invocations() {
                const workflow = getConfiguredWorkflow();
                if (this.auth_fail) return [];
                if (!workflow || !workflow.invocationsFetched) return null;
                return getInvocations(workflow);
            }
        },
        activated() {
            this.init();
        },
        created() {
            this.init(true);
        }
    }
</script>

<style scoped>
    .Jobs >>> thead {
        border-bottom: 1px solid;
    }

    .Jobs >>> .galaxy-workflow-invocation-state {
        display: none;
    }

    .Jobs >>> .galaxy-workflow-output-download > * {
      padding: 0;
      border: none;
    }
</style>
