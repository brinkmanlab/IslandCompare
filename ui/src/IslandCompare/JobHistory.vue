<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <div class="JobHistory">
        <h1>Job History</h1>
        <Jobs v-bind:invocationsPromise="invocationsPromise">
            <template v-slot:functions="slot">
                <template v-if="slot.done && slot.outputs['Results']">
                    <b-link v-bind:to="`/visualize/${slot.outputs['Results'].id}`">Visualize</b-link>
                    <b-link v-bind:href="`/api/histories/${slot.model.history_id}/contents/${slot.outputs['Results'].id}/display?to_ext=gff3&filename=${encodeURIComponent(slot.model.history.name)}.gff3` | auth | galaxybase"><i class="icon-download"></i>Download</b-link>
                    <b-link v-if="slot.outputs['Genomic Islands']" v-bind:href="`/api/histories/${slot.model.history_id}/contents/${slot.outputs['Genomic Islands'].id}/display?to_ext=gff3&filename=${encodeURIComponent(slot.model.history.name)}.gff3` | auth | galaxybase"><i class="icon-download"></i>Genomic Islands</b-link>
                </template>
            </template>
        </Jobs>
    </div>
</template>

<script>
    import Jobs from "@/components/Jobs";
    import {getConfiguredWorkflow, getInvocations} from "@/app";
    import { getUUID } from "@/auth";

    export default {
        name: "JobHistory",
        components: { Jobs },
        data() {return{
        }},
        computed: {
            invocationsPromise() {
                if (getUUID()) return getInvocations(getConfiguredWorkflow());
                else return Promise.resolve([]); // eslint-disable-line vue/no-async-in-computed-properties
            },
        },
        activated() {
            // Force uuid into url when navigating to this page
            if (!('uuid' in this.$route.query)) {
                this.$router.replace({query: {uuid: getUUID()}});
            }
        },
    }
</script>

<style scoped>
    .JobHistory {
        padding-left: 1vw;
        padding-right: 1vw;
        padding-bottom: 1vh;
    }

    .Jobs >>> thead {
        border-bottom: 1px solid;
    }

    .Jobs >>> .galaxy-workflow-invocation-state {
        display: none;
    }
</style>
