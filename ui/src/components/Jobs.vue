<template>
    <table class="Jobs">
        <thead>
        <slot name="header"/>
        <tr>
            <th>Workflow output history name</th>
            <th>Status</th>
            <th>Updated</th>
            <th></th>
        </tr>
        </thead>
        <WorkflowInvocations v-bind:workflows="workflows"/>
    </table>
</template>

<script>
    import WorkflowInvocations from "@/components/workflows/WorkflowInvocations";
    export default {
        name: "Jobs",
        components: {
            WorkflowInvocations,
        },
        props: {
            workflows: {
                type: Array,
                required: true,
            }
        },
        data: ()=>{return {
        }},
        computed: {
        },
        methods: {
        },
        mounted() {
        },
    }
</script>

<style scoped>
    .Jobs >>> .WorkflowInvocations {
        display: table-row-group;
    }

    .Jobs >>> .WorkflowInvocation {
        display: table-row;
        width: auto;
        clear: both;
    }

    .Jobs >>> .WorkflowInvocation > * {
        display: table-cell;
        text-align: right;
    }
</style>