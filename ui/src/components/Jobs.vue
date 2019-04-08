<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
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
            <tr v-if="!invocations.length"><td colspan="4" style="text-align: center">No jobs found</td></tr>
            <WorkflowInvocation v-for="invocation of invocations"
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
    import * as galaxy from '@/galaxy'
    import WorkflowInvocation from "./workflows/WorkflowInvocation";
    export default {
        name: "Jobs",
        components: {
            WorkflowInvocation,
        },
        props: {
            workflow: {
                type: galaxy.workflows.StoredWorkflow,
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
            invocations() {
                return galaxy.workflows.WorkflowInvocation.query().has('history').with('history', q=>q.where('deleted', false)).with('workflow').where('workflow_id', this.workflow.id).get();
                //galaxy.histories.History.query().where('deleted', false).where('tags', tags=>tags.includes(this.workflow.id)).get();
            },
        },
        methods: {
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

    .Jobs >>> .WorkflowInvocation > *, .Jobs >>> .WorkflowInvocation .History > * {
        display: table-cell;
        text-align: center;
    }
</style>