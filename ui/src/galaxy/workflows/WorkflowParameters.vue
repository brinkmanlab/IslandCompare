<template>
    <b-container class="galaxy-workflow-parameters">
        <b-row no-gutters>
            <slot name="before" v-bind="this"></slot>
        </b-row>

        <!-- Workflow inputs -->
        <b-row v-for="(input, index) of workflow_inputs"
               v-bind:key="index"
               no-gutters
        >
            <WorkflowParameter v-bind="input"
                               v-bind:historyPromise="historyPromise"
                               @input="onInput"
                               no-body
                               ref="params"
                               v-bind:class="input.label.replace(' ', '-')"
            />
        </b-row>

        <!-- Workflow Steps -->
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

        <!-- Loading Spinner -->
        <div v-if="workflow_inputs.length === 0" class="text-center">
            <b-spinner class="align-middle"></b-spinner>
            <strong>Loading...</strong>
        </div>

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
        computed: {
            workflow_inputs() {
                let inputs = [];
                for (const [index, input] of Object.entries(this.workflow.inputs)) {
                    const param = {index: index, uuid: input.uuid, label: input.label, type: this.workflow.steps[index].type, annotation: '', optional: false, order: 0};
                    if (this.workflow.steps[index].annotation) {
                        // TODO Temporary until https://github.com/galaxyproject/galaxy/issues/7496
                        // TODO Temporary until galaxy has native support for optional workflow inputs
                        try {
                            Object.assign(param, JSON.parse(this.workflow.steps[index].annotation));
                            if (param.subinputs && param.subinputs.length) {
                                let i = 0;
                                for (; i < param.subinputs.length - 1; ++i) {
                                    inputs.push({...param, ...param.subinputs[i]})
                                }
                                Object.assign(param, param.subinputs[i]);
                            }
                        } catch (e) {
                            if (!(e instanceof SyntaxError)) throw e;
                        }
                    }
                    inputs.push(param);
                }
                inputs = inputs.sort((a,b)=>(a.order || 0) - (b.order || 0));
                return inputs;
            },
        },
        methods: {
            annotation(index) {
                if (this.workflow.steps)
                    return this.workflow.steps[index].annotation || '';
                return '';
            },
            onInput(value) {
                for (const [index, input] of Object.entries(value)) {
                    if (input.subinput) {
                        // TODO remove when Galaxy optional workflow inputs implemented
                        // nest subinput under input for later processing
                        if (!this.inputs[index]) {
                            this.inputs[index] = {
                                src: 'new_collection',
                                collection_type: 'list',
                                name: this.workflow.inputs[index].label,
                                element_identifiers: []
                            };
                        }

                        const subindex = this.inputs[index].element_identifiers.findIndex(i=>i.name === input.subinput);
                        if (Object.entries(input).length === 1) {
                            // Null selection from DatasetParameter
                            if (subindex !== -1) delete this.inputs[index].element_identifiers[subindex];
                            if (this.inputs[index].element_identifiers.length === 0 || this.inputs[index].element_identifiers.every(i=>i === undefined)) delete this.inputs[index];
                        } else {
                            // Add/Update subinput
                            delete input.subinput;
                            if (subindex !== -1) this.inputs[index].element_identifiers[subindex] = input;
                            else this.inputs[index].element_identifiers.push(input);
                        }
                    } else if (input === null) {
                        delete this.inputs[index];
                    } else {
                        this.inputs[index] = input;
                    }
                }

                //Object.assign(this.inputs, value);
                this.$emit('input', this.inputs);
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
                this.$refs.params.forEach(param=>param.reset());
            }
        },
    }
</script>

<style scoped>

</style>