<template>
    <div class="JobRunner">
        <HistoryContents ref="history_contents" v-bind:model="history" v-bind:upload_callback="upload_callback"/>
        <label class="btn btn-primary UploadButton" v-if="history">Upload datasets
            <input type="file" hidden multiple
                   v-if="history"
                   v-bind:accept="permitted_file_extensions.map(ext=>'.'+ext)"
                   @input.prevent="evt=>$refs.history_contents.uploadHandler({dataTransfer:{files: evt.target.files}})"
            />
        </label>
        <div class="WorkflowParams">
            <label v-b-popover.hover="'Label your this job to identify it from others.'">Analysis label<b-form-input type="text" name="invocation_name" ref="invocation_name" v-model="invocation_name" required/></label>
            <slot name="workflow_params" v-bind="params"/>
        </div>
        <b-button @click.prevent="submit()" class="submit" variant="primary">Submit</b-button>
    </div>
</template>

<script>
    import * as galaxy from '@/galaxy'
    import HistoryContents from './histories/HistoryContents';
    import { permitted_file_extensions } from "@/app.config";
    import { invokeConfiguredWorkflow } from "@/app";

    export default {
        name: "JobRunner",
        components: {
            HistoryContents
        },
        props: {
            workflow: [galaxy.workflows.StoredWorkflow, null],
            history: [galaxy.histories.History, null],
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
            upload_callback: {
                type: Function,
                default: u=>u,
            }
        },
        data() {
            return {
                invocation_name: "",
                params: this.workflow_params,
                permitted_file_extensions: permitted_file_extensions,
            }
        },
        methods: {
            submit() {
                if (!this.$refs['invocation_name'].reportValidity()) return;
                let selected = this.$refs.history_contents.getSelectedItems();
                this.$refs.history_contents.clearSelectedItems();

                let error = this.selection_validator(selected);
                if (error) {
                    throw error;
                }

                invokeConfiguredWorkflow(selected, this.invocation_name, this.params);
            }
        },
    }
</script>

<style scoped>
    .JobRunner {
        height: 100%;
        display: grid;
        grid-template-areas: "history_contents history_contents history_contents" "upload params submit";
        grid-template-rows: minmax(20em, auto) auto;
        grid-template-columns: min-content auto min-content;
    }

    .JobRunner .HistoryContents {
        grid-area: history_contents;
        border: solid 1px lightgray;
    }

    .JobRunner .WorkflowParams {
        grid-area: params;
        display: flex;
        flex-direction: row;
        justify-content: space-between;
        flex-wrap: wrap;
        margin-left: 1em;
        margin-right: 1em;
        align-items: flex-start;
        height: min-content;
    }

    .JobRunner .WorkflowParams label input {
        margin-left: 1em;
    }

    .JobRunner .WorkflowParams label input[type="number"] {
        width: 5em;
    }

    .JobRunner .UploadButton {
        grid-area: upload;
        align-self: start;
        white-space: nowrap;
    }

    .JobRunner .submit {
        grid-area: submit;
        align-self: start;
    }
</style>