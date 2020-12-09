<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <div class="JobHistory">
        <h1>Job History</h1>
        <Jobs v-bind:invocations="invocations">
            <template v-slot:functions="slot">
                <template v-if="slot.done && slot.model.outputs['Results']">
                    <b-link v-bind:to="`/visualize/${slot.model.outputs['Results'].id}`">Visualize</b-link>
                    <b-link v-bind:href="$filters.galaxybase($filters.auth(`/api/histories/${slot.model.history_id}/contents/${slot.model.outputs['Results'].id}/display?to_ext=gff3`))"><i class="icon-download"></i>Download</b-link>
                    <b-link v-if="slot.model.outputs['Genomic Islands']" v-bind:href="$filters.galaxybase($filters.auth(`/api/histories/${slot.model.history_id}/contents/${slot.model.outputs['Genomic Islands'].id}/display?to_ext=gff3&filename=${encodeURIComponent(slot.model.history.name)}.gff3`))"><i class="icon-download"></i>Genomic Islands</b-link>
                    <b-link v-if="slot.model.outputs['Newick']" v-bind:href="$filters.galaxybase($filters.auth(`/api/histories/${slot.model.history_id}/contents/${slot.model.outputs['Newick'].id}/display?to_ext=newick`))"><i class="icon-download"></i>Phylo Tree</b-link>
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

    export default {
        name: "JobHistory",
        components: { Jobs },
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
    .Jobs :deep( thead ) {
        border-bottom: 1px solid;
    }

    .Jobs :deep( .galaxy-workflow-invocation-state ) {
        display: none;
    }
</style>
