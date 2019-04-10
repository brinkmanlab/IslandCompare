<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <div class="Jobs">
        <div class="header">
            <slot name="header"/>
            <div>
                <span>Label</span>
                <span>Status</span>
                <span>Updated</span>
                <span></span>
            </div>
        </div>
        <div>
            <div v-if="!invocations.length"><span style="grid-column: 1 / end;">No jobs found</span></div>
            <WorkflowInvocation v-for="invocation of invocations"
                                v-bind:key="invocation.id"
                                v-bind:model="invocation"
            >
                <template v-slot:functions="slot">
                    <slot name="functions" v-bind="slot" />
                </template>
            </WorkflowInvocation>
        </div>
    </div>
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
    .Jobs {
        display: grid;
        grid-template-columns: [start label] auto [status] minmax(10em, min-content) [updated] minmax(10em, auto) [functions] minmax(10em, auto) [end];
        grid-auto-rows: min-content;
        grid-row-gap: 0.5em;
    }

    .Jobs > *, .Jobs > * > *, .WorkflowInvocation, .WorkflowInvocation >>> .History {
        display: contents;
    }

    .Jobs .header > :first-child {
        display: block;
        grid-column: start / end;
    }

    .Jobs > * > * > *, .WorkflowInvocation >>> .History > * {
        display: block;
        text-align: center;
    }

    .Jobs >>> .functions {
        text-align: right;
    }

    .Jobs >>> .functions > :not(:last-child) {
        margin-right: 2em;
    }
</style>