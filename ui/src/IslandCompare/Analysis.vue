<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <JobManager>
        <template v-slot:workflow_params="params">
            <label>Minimum island size<input type="number" min="0" required v-model.number.lazy="params.minimum_island_size" /></label>
            <label>Minimum homologous region<input type="number" min="0" required v-model.number.lazy="params.minimum_homologous_region" /></label>
            <label>Minimum cluster size<input type="number" min="1" required v-model.number.lazy="params.minimum_cluster_size" /></label>
        </template>
        <template v-slot:invocation_functions="slot">
            <template v-if="slot.done && slot.outputs['IslandCompare Result']">
                <a v-bind:href="`/plugins/visualizations/islandcompare/show?dataset_id=${slot.outputs['IslandCompare Result'].id}` | auth">Visualize</a>
                <a v-bind:href="`/datasets/${slot.outputs['IslandCompare Result'].id}/display?to_ext=gff3` | auth">Download</a>
            </template>
        </template>
    </JobManager>
</template>

<script>
    import JobManager from "@/components/JobManager";
    export default {
        name: "Analysis",
        components: {
            JobManager
        },
        methods: {
        }
    }
</script>

<style scoped>

</style>