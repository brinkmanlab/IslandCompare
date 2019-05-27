<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <JobManager v-if="user_id"
                v-bind:permitted_file_extensions="permitted_file_extensions"
                v-bind:selection_validator="selection=>selection.length<2?'You must select more than one dataset for comparison':null"
    >
        <template v-slot:workflow_params="params">
            <label v-b-popover.hover="'Filter detected islands by size'">Minimum island size<b-form-input type="number" min="0" required v-model.number.lazy="params.minimum_island_size" /></label>
            <label v-b-popover.hover="'Filter alignments by size. Small alignments are usually meaningless and only bloat the resulting dataset.'">Minimum homologous region<b-form-input type="number" min="0" required v-model.number.lazy="params.minimum_homologous_region" /></label>
            <label v-b-popover.hover="'Filter clusters by cardinality.'">Minimum cluster size<b-form-input type="number" min="1" required v-model.number.lazy="params.minimum_cluster_size" /></label>
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