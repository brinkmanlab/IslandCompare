<template>
    <b-container class="galaxy-workflow-parameters">
        <b-row no-gutters>
            <slot name="before" v-bind="this"></slot>
        </b-row>
        <b-row v-for="(input, index) of workflow.inputs"
               v-bind:key="index"
               ref="params"
               no-gutters
        >
            <WorkflowParameter v-bind:label="input.label"
                               v-bind:type="type(index)"
                               v-bind:uuid="input.uuid"
                               v-bind:annotation="annotation(index)"
                               v-bind:historyPromise="historyPromise"
                               @input="onInput"
                               @upload="upload=>$emit('upload', upload)"
                               no-body
            />
        </b-row>
        <b-row no-gutters class="galaxy-workflow-parameters-steps">
            <b-card no-body>
                <b-card-header header-tag="header" role="tab">
                    <b-button block href="#" v-b-toggle="`${_uid}-step-parameters`" variant="info">Step Parameters</b-button>
                </b-card-header>
                <b-collapse v-bind:id="`${_uid}-step-parameters`">
                    <!-- Step parameters -->
                    <b-card v-for="(step, index) of workflow.steps" v-bind:key="index">
                        <!-- TODO need to load tool info to get name and load stored workflow -->
                        <b-card-header header-tag="header">{{step.tool_id}}</b-card-header>
                        <!-- TODO move this out to a tool view -->
                        <!--WorkflowParameter v-for="[label, input] of step.tool_inputs" v-bind:key="label" /-->
                    </b-card>
                </b-collapse>
            </b-card>
        </b-row>
        <b-row>
            <slot name="after" v-bind="this"></slot>
        </b-row>
    </b-container>
</template>

<script>
    import WorkflowParameter from "./WorkflowParameter";

    export default {
        name: "WorkflowParameters",
        components: {WorkflowParameter},
        props: {
            workflowPromise: {
                type: Promise,
                required: true,
            },
            historyPromise: {
                type: Promise,
                required: true,
            },
        },
        data() {return {
            inputs: {},
        }},
        asyncComputed: {
            workflow: {
                get() { return this.workflowPromise },
                default: {inputs: [], steps: []},
            },
        },
        methods: {
            type(index) {
                // TODO Temporary until https://github.com/galaxyproject/galaxy/issues/7496
                if (this.workflow.steps[index].type === 'parameter_input') {
                    if ('__dbkey__' in this.workflow.steps[index].annotation) return 'dbkey_input';
                    return 'parameter_input';
                } else {
                    return this.workflow.steps[index].type;
                }
            },
            annotation(index) {
                if (this.workflow.steps)
                    return this.workflow.steps[index].annotation || '';
                return '';
            },
            onInput(value) {
                Object.assign(this.inputs, value);
            },
            checkValidity() {
                /* Returns a Boolean that is false if the element is a candidate for constraint validation,
                and it does not satisfy its constraints. In this case, it also fires an invalid event at the element.
                It returns true if the element is not a candidate for constraint validation, or if it satisfies
                its constraints.
                 */
                return this.$refs.params.every(param=>param.checkValidity());
            },
            reportValidity() {
                /* Runs the checkValidity() method, and if it returns false (for an invalid input or no pattern
                attribute provided), then it reports to the user that the input is invalid in the same manner
                as if you submitted a form.
                 */
                return this.$refs.params.every(param=>param.reportValidity());
            },
            reset() {
                /* Reset all parameters to their default state */
                this.$refs.params.each(param=>param.reset());
            }
        },
    }
</script>

<style scoped>

</style>