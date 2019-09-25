<template>
    <b-container class="JobRunner">
        <b-row>
            <b-input-group prepend="Analysis Label">
                <b-form-input ref="invocation_name" v-model="invocation_name" placeholder="Label this job to identify it from others" required />
                <b-input-group-append>
                    <b-button text="Submit" variant="success" @click="submit">Submit</b-button>
                </b-input-group-append>
            </b-input-group>
        </b-row>
        <WorkflowParameters ref="workflow_parameters"
                            v-bind:historyPromise="historyPromise"
                            v-bind:workflowPromise="workflowPromise"
                            @input="onInput"
        />
    </b-container>
</template>

<script>
    import WorkflowParameters from "@/galaxy/src/workflows/WorkflowParameters";

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
            selection_validator: {
                type: Function,
                default: selection => selection.length === 0 ? "Invalid dataset selection" : null,
            },
        },
        data() {
            return {
                invocation_name: '',
                params: {},
            }
        },
        asyncComputed: {
            //history() { return this.historyPromise },
            workflow() { return this.workflowPromise },
        },
        methods: {
            submit() {
                if (!this.$refs.invocation_name.reportValidity()) return;  // Trigger label field validator
                if (!this.$refs.workflow_parameters.reportValidity()) return;  // Trigger workflow params validator

                // Vue observers cause a race condition where the values are reset before invoke can access them.
                // The following causes a deep copy to avoid that
                const params = JSON.parse(JSON.stringify(this.params));

                this.workflow.invoke(params, undefined, this.invocation_name);
                this.$refs.workflow_parameters.reset();
                this.invocation_name = '';
                this.$emit('galaxy-workflow-invocation');
            },
            onInput(value) {
                this.params = value;
            },
        },
    }
</script>

<style scoped>

</style>