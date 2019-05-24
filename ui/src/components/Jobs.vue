<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <table class="Jobs">
        <thead>
            <tr><th colspan="4"><slot name="header"/></th></tr>
            <tr class="header">
                <th>Label</th>
                <th>Updated</th>
                <th>Status</th>
                <th></th>
            </tr>
        </thead>
        <tbody>
        <tr v-if="workflow === null"><td colspan="4">Loading jobs</td></tr>
        <tr v-else-if="!invocations.length"><td colspan="4">No jobs found</td></tr>
        <WorkflowInvocation v-for="invocation of invocations"
                            v-bind:key="invocation.id"
                            v-bind:model="invocation"
                            v-bind:class="row_class(invocation.aggregate_state())"

                >
                <template v-slot:functions="slot">
                    <slot name="functions" v-bind="slot"/>
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
            workflow: [galaxy.workflows.StoredWorkflow, null],
            orderBy: {
                type: Function,
                default: (a,b)=>(Date.parse(a)-Date.parse(b)),
            },
        },
        data: ()=>{return {
            col_names: ['galaxy-history-label', 'galaxy-workflow-invocation-state', 'galaxy-history-updated', 'galaxy-history-functions']
        }},
        computed: {
            invocations() {
                if (this.workflow === null) return [];
                return galaxy.workflows.WorkflowInvocation.query().has('history').with('history', q=>q.where('deleted', false)).with('workflow').where('workflow_id', this.workflow.id).with('steps').get();
                //galaxy.histories.History.query().where('deleted', false).where('tags', tags=>tags.includes(this.workflow.id)).get();
            },
        },
        methods: {
            row_class(state) {
                if (state === "new") return "table-primary";
                if (state === "done") return "table-success";
                if (state === "error") return "table-danger";
                if (state === "running") return "table-info";
                return "table-secondary";
            },
            element_formatter(elem) {
                let i = this.col_names.indexOf(elem);
                if (i !== -1) return elem + " col order-" + (i+1).toString();
                return elem + " hidden";
            },
        },
        mounted() {
        },
    }
</script>

<style scoped>
    .Jobs {
        width: 100%;
        /*display: grid;
        grid-template-columns: [start label] auto [status] minmax(10em, min-content) [updated] minmax(10em, auto) [functions] minmax(10em, auto) [end];
        grid-auto-rows: min-content;
        grid-row-gap: 0.5em;*/
    }

    .header {
        text-align: center;
    }

    .galaxy-workflow-invocation {
        display: table-row;
    }

    .galaxy-workflow-invocation >>> .History {
        display: contents;
    }

    .galaxy-workflow-invocation >>> .History, .galaxy-workflow-invocation >>> .History > * {
        background-color: inherit;
    }

    >>> .hidden {
        display: none;
    }
/*
    .Jobs .header > :first-child {
        display: block;
        grid-column: start / end;
    }
    */
    .galaxy-workflow-invocation >>> .History > * {
        display: table-cell;
        text-align: center;
    }
/*
    .Jobs > * > * > * {
        grid-row: 2;
    }
*/
    .galaxy-workflow-invocation >>> .galaxy-history-state {
        display: none;
    }
/*
    .galaxy-workflow-invocation >>> .galaxy-history-label {
        grid-column: label;
        order: 1;
    }

    .galaxy-workflow-invocation >>> .galaxy-workflow-invocation-state {
        grid-column: status;
        order: 2;
    }

    .galaxy-workflow-invocation >>> .galaxy-history-updated {
        grid-column: updated;
        order: 3;
    }

    .galaxy-workflow-invocation >>> .galaxy-history-functions {
        grid-column: functions;
        order: 4;
    }*/

    .Jobs >>> .galaxy-history-functions {
        text-align: right;
    }

    .Jobs >>> .galaxy-history-functions > :not(:last-child) {
        margin-right: 2em;
    }
</style>