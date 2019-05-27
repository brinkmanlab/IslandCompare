<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <JobManager v-if="user_id"
                v-bind:permitted_file_extensions="permitted_file_extensions"
                v-bind:selection_validator="selection=>selection.length<2?'You must select more than one dataset for comparison':null"
    >
        <template v-slot:workflow_params="params">
            <label>Minimum island size<b-form-input type="number" min="0" required v-model.number.lazy="params.minimum_island_size" /></label>
            <label>Minimum homologous region<b-form-input type="number" min="0" required v-model.number.lazy="params.minimum_homologous_region" /></label>
            <label>Minimum cluster size<b-form-input type="number" min="1" required v-model.number.lazy="params.minimum_cluster_size" /></label>
        </template>
        <template v-slot:invocation_functions="slot">
            <template v-if="slot.done">
                <b-link v-bind:href="`/plugins/visualizations/islandcompare/show?dataset_id=${slot.outputs['IslandCompare Result'].id}` | auth">Visualize</b-link>
                <b-link v-bind:href="`/datasets/${slot.outputs['IslandCompare Result'].id}/display?to_ext=gff3` | auth">Download</b-link>
            </template>
        </template>
    </JobManager>
</template>

<script>
    import JobManager from "@/components/JobManager";
    import { permitted_file_extensions } from "@/app.config";
    import { getOrCreateUUID } from "@/auth";

    export default {
        name: "Analysis",
        components: {
            JobManager
        },
        data() {return{
            permitted_file_extensions: permitted_file_extensions,
        }},
        asyncComputed: {
            user_id: getOrCreateUUID,
        },
    }
</script>

<style scoped>

</style>