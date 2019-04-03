<template>
    <table class="Jobs">
        <thead>
        <slot name="header"/>
        <tr>
            <th>Label</th>
            <th>Status</th>
            <th>Updated</th>
            <th></th>
        </tr>
        </thead>
        <tbody>
            <WorkflowInvocation
                    v-for="invocation of invocations"
                    v-bind:key="invocation.id"
                    v-bind:model="invocation"
            >
                <template v-slot:functions="slot">
                    <slot name="functions" v-bind="slot" />
                </template>
            </WorkflowInvocation>
        </tbody>
    </table>
</template>

<script>
    //import * as galaxy from '@/galaxy'
    import WorkflowInvocation from "./workflows/WorkflowInvocation";
    export default {
        name: "Jobs",
        components: {
            WorkflowInvocation,
        },
        props: {
            workflows: {
                type: Array,
                required: true,
            },
            orderBy: {
                type: Function,
                default: (a,b)=>(Date.parse(a)-Date.parse(b)),
            },
        },
        data: ()=>{return {
        }},
        computed: {
            invocations: function() { return this.workflows[0].invocations.concat(...(this.workflows.slice(1).map(a=>a.invocations))).sort(this.orderBy) },
        },
        methods: {
            update() {
                //TODO
                //for (let workflow in this.workflows) {
                //    galaxy.workflows.WorkflowInvocation.$fetch({params: {url: workflow.url}});
                //}
            },
        },
        mounted() {
        },
    }
</script>

<style scoped>
    .Jobs tbody > * {
        display: table-row;
        width: auto;
        clear: both;
    }

    .Jobs >>> .WorkflowInvocation > * {
        display: table-cell;
        text-align: center;
    }
</style>