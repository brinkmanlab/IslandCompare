<template>
    <b-container class="JobRunner">
        <b-row>
            <b-input-group prepend="Project">
                <b-form-select v-model="project" :options="projects" :disabled="this.histories === null" disabled-field="Loading.." />
            </b-input-group>
        </b-row>
        <b-row>
            <b-input-group prepend="Analysis Label">
                <b-form-input ref="invocation_name" class="invocation-name" v-model="invocation_name" placeholder="Label this job to identify it from others" required />
                <b-input-group-append>
                    <b-button text="Submit" variant="success" @click="submit" class="invocation-submit">Submit</b-button>
                </b-input-group-append>
            </b-input-group>
        </b-row>
        <WorkflowParameters ref="workflow_parameters"
                            v-bind:history="history"
                            v-bind:workflow="workflow"
                            @input="onInput"
        />
    </b-container>
</template>

<script>
    import {workflows, api} from "galaxy-client";

    export default {
        name: "JobRunner",
        components: {
            WorkflowParameters: ()=>workflows.WorkflowParameters,
        },
        props: {
            histories: {
                validator(prop) {return prop === null || (prop.every && prop.every(p=>p instanceof api.histories.History))},
                required: true,
            },
            workflow: {
                validator(prop) {return prop === null || prop instanceof api.workflows.StoredWorkflow},
                required: true,
            },
            selection_validator: {
                type: Function,
                default: selection => selection.length === 0 ? "Invalid dataset selection" : null,
            },
        },
        data() {
            return {
                invocation_name: '',
                params: {},
                project: null,
            }
        },
        computed: {
            history() {
                if (!this.histories) return null;
                return this.histories.find(h=>h.id === this.project) || null;
            },
            projects() {
                if (!this.histories) return [];
                return this.histories.map(h=>({text: h.name, value: h.id}));
            }
        },
        methods: {
            submit() {
                if (!this.workflow) throw("Workflow not ready");
                if (!this.$refs.invocation_name.reportValidity()) return;  // Trigger label field validator
                if (!this.$refs.workflow_parameters.reportValidity()) return;  // Trigger workflow params validator

                // TODO https://github.com/brinkmanlab/IslandCompare/issues/87

                // Vue observers cause a race condition where the values are reset before invoke can access them.
                // The following causes a deep copy to avoid that
                const params = JSON.parse(JSON.stringify(this.params));

                const invocation = this.workflow.invoke(params, undefined, this.invocation_name);
                this.$refs.workflow_parameters.reset();
                this.invocation_name = '';
                this.$emit('galaxy-workflow-invocation', invocation);
            },
            onInput(value) {
                this.params = value;
            },
        },
        watch: {
            histories(new_val) {
                // Set current project to most recent updated on load
                if (new_val && new_val.length > 0 && !new_val.some(h=>h.id === this.project)) {
                    this.project = new_val.reduce((acc, cur)=>new Date(cur.update_time) > new Date(acc.update_time) ? cur : acc, new_val[0]).id;
                }
            },
        }
    }
</script>

<style scoped>
    .row, >>> .galaxy-workflow-parameters > .row {
        margin: 0;
        margin-bottom: 1vh;
    }

    .JobRunner {
        padding-top: 1vh;
        padding-bottom: 1vh;
    }

    >>> .galaxy-workflow-parameters > .row {
        //margin-bottom: 3vh;
    }

    >>> .galaxy-workflow-parameters .card-header {
        font-weight: bold;
        height: 2em;
    }

    >>> .galaxy-workflow-parameters {
        display: contents;
    }

    >>> .galaxy-workflow-parameters .card-body {
        padding: 0;
    }

    >>> .galaxy-workflow-parameters .card-header {
        padding: 0 1rem;
    }

    >>> .galaxy-workflow-parameters .card-footer {
        padding: 0 1rem;
        font-size: 0.8em;
    }

    >>> .galaxy-workflow-parameters-steps {
        display: none;
    }

    >>> .galaxy-workflow-parameters > .row > * {
        width: 100%;
    }

    >>> .galaxy-workflow-parameters .table-sm td {
        padding: 0;
    }
</style>
