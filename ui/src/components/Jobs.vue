<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <!-- TODO change to bootstrap table -->
    <table class="Jobs">
        <thead>
            <tr><th colspan="4"><slot name="header"/></th></tr>
            <tr class="header">
                <th>Label</th>
                <th>Status</th>
                <th>Updated</th>
                <th></th>
            </tr>
        </thead>
        <tbody>
            <tr v-if="invocations === null"><td colspan="4"><div class="text-center loading-spinner">
                <b-spinner class="align-middle"></b-spinner>
                <strong>Loading...</strong>
            </div></td></tr>
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
    import {workflows} from "galaxy-client";

    export default {
        name: "Jobs",
        components: {
            WorkflowInvocation: ()=>workflows.WorkflowInvocation,
        },
        props: {
            invocations: {
                validator(prop) {return prop instanceof Array || prop === null},
                required: true,
            },
        },
        data: ()=>{return {
            col_names: ['galaxy-history-label', 'galaxy-workflow-invocation-state', 'galaxy-history-updated', 'galaxy-history-functions'],
        }},
        methods: {
            row_class(state) {
                if (state === "new") return "table-primary";
                if (state === "done") return "table-success";
                if (state === "error" || state === "failed") return "table-danger";
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

    .galaxy-workflow-invocation >>> .galaxy-history > * {
        display: table-cell;
        text-align: center;
    }

    .galaxy-workflow-invocation >>> .galaxy-workflow-invocation-state:not(.new) {
        display: none;
    }

    .galaxy-workflow-invocation >>> .galaxy-workflow-invocation-progress.new {
        display: none;
    }

    .Jobs >>> .galaxy-history-functions {
        text-align: right;
      white-space: nowrap;
    }

    .Jobs >>> .galaxy-history-functions > :not(:last-child) {
        margin-right: 2em;
    }

    .loading-spinner {
        margin: 1em;
    }

    .loading-spinner strong {
        padding-left: 1em;
    }
</style>
