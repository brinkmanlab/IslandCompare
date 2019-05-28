<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <div class="JobHistory">
        <h1>Job History</h1>
        <Jobs v-bind:workflow="workflow">
            <template v-slot:functions="slot">
                <template v-if="slot.done">
                    <b-link v-bind:href="`/plugins/visualizations/islandcompare/show?dataset_id=${slot.outputs['IslandCompare Result'].id}` | auth">Visualize</b-link>
                    <b-link v-bind:href="`/datasets/${slot.outputs['IslandCompare Result'].id}/display?to_ext=gff3` | auth">Download</b-link>
                </template>
            </template>
        </Jobs>
    </div>
</template>

<script>
    import Jobs from "@/components/Jobs";
    import { getConfiguredWorkflow } from "@/app";

    export default {
        name: "JobHistory",
        components: {Jobs},
        data() {return{
        }},
        asyncComputed: {
            workflow: getConfiguredWorkflow,
        }
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
</style>