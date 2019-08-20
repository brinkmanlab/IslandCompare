<template>
    <b-container class="JobRunner">
        <b-row>
            <b-input-group prepend="Analysis Label">
                <b-form-input ref="invocation_name" v-model="invocation_name" placeholder="Label this job to identify it from others" />
                <b-input-group-append>
                    <b-button text="Submit" variant="success" @click="submit">Submit</b-button>
                </b-input-group-append>
            </b-input-group>
        </b-row>
        <WorkflowParameters ref="workflow_parameters"
                            v-bind:historyPromise="historyPromise"
                            v-bind:workflowPromise="workflowPromise"
        />
    </b-container>
</template>

<script>
    import { invokeConfiguredWorkflow } from "@/app";
    import WorkflowParameters from "@/galaxy/workflows/WorkflowParameters";

    export default {
        name: "JobRunner",
        components: {
            WorkflowParameters,
        },
        props: {
            historyPromise: {
                type: Promise,
                required: true,
            },
            workflowPromise: {
                type: Promise,
                required: true,
            },
            workflow_params: {
                type: Object,
                default() {
                    return {}
                },
            },
            selection_validator: {
                type: Function,
                default: selection => selection.length === 0 ? "Invalid dataset selection" : null,
            },
        },
        data() {
            return {
                invocation_name: "",
                params: this.workflow_params,
            }
        },
        asyncComputed: {
            history() { return this.historyPromise },
        },
        methods: {
            submit() {
                if (!this.$refs.invocation_name.reportValidity()) return;  // Trigger label field validator
                if (!this.$refs.workflow_parameters.reportValidity()) return;  // Trigger workflow params validator

                this.$refs.workflow_parameters.reset();

                //TODO move to StoredWorkflow.invoke()
                let invocation = invokeConfiguredWorkflow(selected, this.invocation_name, this.params);
                this.$emit('galaxy-workflow-invocation', invocation);
            },
            onInput(val, prop) {
                //v-model cant bind via slot, this implements explicitly
                this.params[prop] = val;
            },
        },
    }
</script>

<style scoped>

</style>