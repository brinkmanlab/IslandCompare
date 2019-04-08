<template>
    <div class="JobRunner">
        <HistoryContents v-if="history" ref="history_contents" v-bind:model="history" filter="dataset" @operation-fail="e=>this.$emit('toast', e)"/>
        <label class="UploadButton" v-if="history">Upload datasets
            <input type="file" v-if="history" v-show="false" @input.prevent="evt=>$refs.history_contents.uploadHandler({dataTransfer:{files: evt.target.files}})"/>
        </label>
        <div class="WorkflowParams">
            <label>Invocation label<input type="text" name="invocation_name" v-model="invocation_name"/></label>
            <slot name="workflow_params" v-bind="params"/>
        </div>
        <input type="submit" @click.prevent="submit()"/>
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
            workflow: {
                type: galaxy.workflows.StoredWorkflow,
                required: true,
            },
            history: {
                type: galaxy.histories.History,
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
            }
        },
        data() {
            return {
                invocation_name: this.workflow.name,
                params: this.workflow_params,
            }
        },
        computed: {},
        methods: {
            async submit() {
                let selected = this.$refs.history_contents.getSelectedItems();
                let error = this.selection_validator(selected);
                if (error) {
                    this.$emit('toast', error);
                    throw error;
                }

                //Create history to store run
                let response = await galaxy.histories.History.$create({
                    data: {
                        name: this.invocation_name,
                    }
                });

                let run_history = galaxy.histories.History.find(response.id);
                if (!run_history) throw "Failed to create a invocation history.";
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
    }

    .JobRunner .WorkflowParams {
        grid-area: params;
        display: flex;
        flex-direction: row;
        justify-content: space-around;
        flex-wrap: wrap;
    }

    .JobRunner .UploadButton {
        grid-area: upload;
        display: inline-block;
        color: #444;
        border: 1px solid #CCC;
        background: #DDD;
        box-shadow: 0 0 5px -1px rgba(0, 0, 0, 0.2);
        cursor: pointer;
        vertical-align: middle;
        max-width: 100px;
        padding: 5px;
        text-align: center;
    }

    .JobRunner .UploadButton:active {
        box-shadow: 0 0 5px -1px rgba(0, 0, 0, 0.6);
    }

    .JobRunner input[type=submit] {
        grid-area: submit;
    }
</style>