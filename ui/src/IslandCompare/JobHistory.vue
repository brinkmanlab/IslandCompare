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
    import {updateRoute} from "@/auth";

    export default {
        name: "JobHistory",
        components: { Jobs },
        data() {return{
        }},
        asyncComputed: {
            invocationsPromise: {
                async get() {
                    if (await getUUID()) return getInvocations(getConfiguredWorkflow());
                    return Promise.resolve([]);
                },
                default: Promise.resolve([]),
            }
        },
        activated() {
            // Force uuid into url when navigating to this page
            updateRoute(this.$router, this.$route);
        },
    }
</script>

<style scoped>
    .Jobs >>> thead {
        border-bottom: 1px solid;
    }

    .Jobs >>> .galaxy-workflow-invocation-state {
        display: none;
    }
</style>
