<template>
    <div class="WorkflowInvocations">
        <WorkflowInvocation
            v-for="invocation of invocations"
            v-bind:key="invocation.id"
            v-bind:model="invocation"
        />
    </div>
</template>

<script>
    import * as galaxy from '@/galaxy'
    import WorkflowInvocation from "./WorkflowInvocation";
    export default {
        name: "WorkflowInvocations",
        components: { WorkflowInvocation },
        props: {
            workflows: {
                type: Array,
                required: true,
                validator: value=>value.every(a=>a instanceof galaxy.workflows.StoredWorkflow),
            },
            orderBy: {
                type: Function,
                default: (a,b)=>(Date.parse(a)-Date.parse(b)),
            },
        },
        computed: {
            invocations: function() { return this.workflows[0].invocations.filter(a=>a.history).concat(...(this.workflows.slice(1).map(a=>a.invocations.filter(a=>a.history)))).sort(this.orderBy) },
        },
        mounted() {
        },
    }
</script>

<style scoped>
    table {
        width: 100% ;
    }
</style>