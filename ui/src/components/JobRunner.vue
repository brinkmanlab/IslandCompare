<template>
    <div class="JobRunner">
        <HistoryContents ref="history_contents" v-bind:model="history" v-bind:upload_callback="upload_callback"/>
        <label class="btn btn-primary UploadButton" v-if="history">Upload datasets
            <input type="file" v-if="history" hidden @input.prevent="evt=>$refs.history_contents.uploadHandler({dataTransfer:{files: evt.target.files}})"/>
        </label>
        <div class="WorkflowParams">
            <label>Analysis label<b-form-input type="text" name="invocation_name" v-model="invocation_name"/></label>
            <slot name="workflow_params" v-bind="params"/>
        </div>
        <b-button @click.prevent="submit()" class="submit">Submit</b-button>
    </div>
</template>

<script>
    import * as galaxy from '@/galaxy'
    import HistoryContents from './histories/HistoryContents';

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
            }
        },
        watch: {
            workflow(val, oldVal) {
                if (oldVal === null && val && this.invocation_name === "") this.invocation_name = val.name;
            },
        },
        computed: {},
        methods: {
            async submit() {
                //Invoke workflow, creating a new history and adding selected datasets to a collection first

                let selected = this.$refs.history_contents.getSelectedItems();
                this.$refs.history_contents.clearSelectedItems();

                let error = this.selection_validator(selected);
                if (error) {
                    throw error;
                }

                //Create history to store run
                let response = await galaxy.histories.History.$create({
                    data: {
                        name: this.invocation_name,
                    }
                });

                let run_history = galaxy.histories.History.find(response.id);
                if (!run_history) {
                    error = "Failed to create a invocation history.";
                    throw error;
                }
                run_history.tags.push(this.workflow.id);
                run_history.upload();

                //Create collection of inputs in new history
                response = await galaxy.history_contents.HistoryDatasetCollectionAssociation.$create({
                    params: {
                        url: run_history.contents_url,
                    },
                    data: {
                        name: "Selected datasets",
                        type: 'dataset_collection',
                        collection_type: 'list',
                        element_identifiers: selected.map(model => ({
                            src: (model instanceof galaxy.history_contents.HistoryDatasetAssociation) ? model.hda_ldda : 'hdca', //TODO else 'hdca' is fragile
                            name: model.hid + model.name,
                            id: model.id,
                        })),
                    }
                });
                //let params = Object.entries(this.params).reduce((a,[k,v])=>{a[k]=(v.toString ? v.toString() : v); return a}, {}); //https://github.com/galaxyproject/galaxy/issues/7654
                //Invoke workflow
                try {
                    response = await galaxy.workflows.WorkflowInvocation.$create({
                        params: {
                            url: this.workflow.url,
                        },
                        data: {
                            //workflow_id: this.workflow.id,
                            history_id: run_history.id,
                            parameters: {
                                //TODO dynamically generate simple inputs, currently hardcoded to islandcompare
                                0:{"input":{"values":[{"src":"hdca","tags":[],"hid":8,"id":response.id}],"batch":false}},
                                13: {
                                    cond: `c5-c4>${this.params.minimum_island_size}`,
                                    "header_lines": "0"
                                },
                                27: {
                                    "envs_0|name": "minimum_homologous_region",
                                    "envs_0|val": this.params.minimum_homologous_region.toString(),
                                    "envs_1|name": "min_cluster_size",
                                    "envs_1|val": this.params.minimum_cluster_size.toString(),
                                }
                            },
                            parameters_normalized: true,
                            batch: true,
                            no_add_to_history: true,
                        }
                    });
                    //galaxy.workflows.WorkflowInvocation.find(response.id).start_polling();
                } catch (e) {
                    galaxy.histories.History.$delete({params: {id: run_history.id}});
                    throw e;
                }

                this.$emit('workflow-invoked', galaxy.workflows.WorkflowInvocation.find(response.id));
            },
        },
        mounted() {
        },
    }
</script>

<style scoped>
    .JobRunner {
        height: 100%;
        display: grid;
        grid-template-areas: "history_contents history_contents history_contents" "upload params submit";
        grid-template-rows: minmax(30em, auto) auto;
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