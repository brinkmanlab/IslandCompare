<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <!-- TODO change to bootstrap table -->
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
            <tr v-if="invocations === null"><td colspan="4">Loading jobs</td></tr>
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
    import WorkflowInvocation from "@/galaxy/src/workflows/WorkflowInvocation";
    import { WorkflowInvocation as WorkflowInvocationModel } from "@/galaxy/src/api/workflows";

    export default {
        name: "Jobs",
        components: {
            WorkflowInvocation,
        },
        props: {
            invocationsPromise: {
                type: Promise,
                required: true,
            },
            orderBy: {
                type: Function,
                default: (a,b)=>(Date.parse(a)-Date.parse(b)),
            },
        },
        data: ()=>{return {
            col_names: ['galaxy-history-label', 'galaxy-workflow-invocation-state', 'galaxy-history-updated', 'galaxy-history-functions'],
        }},
        asyncComputed: {
            async initial_invocations() {
                return this.invocationsPromise; // TODO find way to make returned value reactive
                //return WorkflowInvocationModel.query().whereHas('history', q => q.where('deleted', false)).with('history|workflow|steps.jobs').get();
            },
        },
        computed: {
            invocations() {
                // This is in place of initial_invocations as the full query needs to be within the computed property :(
                if (this.initial_invocations === null) return null;
                return WorkflowInvocationModel.query().whereHas('history', q => q.where('deleted', false)).with('history|workflow|steps.jobs').get();
            }
        },
        methods: {
            row_class(state) {
                if (state === "new") return "table-primary";
                if (state === "done") return "table-success";
                if (state === "error") return "table-danger";
                if (state === "running") return "table-info";
                return "table-secondary";
            },
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

    .galaxy-workflow-invocation >>> .galaxy-history {
        display: contents;
    }

    .galaxy-workflow-invocation >>> .galaxy-history, .galaxy-workflow-invocation >>> .galaxy-history > * {
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
    .galaxy-workflow-invocation >>> .galaxy-history > * {
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