<template>
    <b-container class="galaxy-history-contents">
        <b-row>
            <b-input-group size="sm" prepend="Filter">
                <b-form-input
                    v-model="user_filter"
                    type="search"
                    placeholder="Type to Search"
                />
                <b-input-group-append>
                    <b-button :disabled="!user_filter" @click="user_filter = ''">Clear</b-button>
                </b-input-group-append>
            </b-input-group>

            <FunctionIcon label="Upload" description="Select a dataset to upload" icon="icon-file-upload" v-bind:action="show_upload" />
            <input type="file" hidden multiple
                   ref="upload"
                   v-bind:accept="permitted_file_extensions.map(ext=>'.'+ext)"
                   @change.prevent="uploadHandler"
            />
        </b-row>

        <b-row>
            <b-table
                     show-empty
                     small
                     borderless
                     stacked="md"
                     per-page="0"
                     :items="items"
                     :busy="isLoading"
                     selectable
                     select-mode="range"
                     :fields="['item']"
                     :filter="user_filter"
                     :filter-function="filter"
                     empty-text="Drag and drop files here to upload"
                     thead-class="hidden_header"
                     @row-selected="items=>$emit('input', items.map(item=>item.item))"
                     @dragstart.native.stop.prevent @dragover.native.prevent="upload_dragging=true" @dragleave.native="upload_dragging=false" @dragexit.native="upload_dragging=false" @drop.native.prevent="uploadHandler"
                     ref="table"
            >
                <template slot="[item]" slot-scope="row">
                    <DatasetItem
                            v-if="row.value.history_content_type==='dataset'"
                            v-bind:model="row.value"
                            v-bind:class="{'table-danger': row.value.hid === -1}"
                    />
                    <CollectionItem
                            v-else-if="row.value.history_content_type==='dataset_collection'"
                            v-bind:model="row.value"
                    />
                </template>

                <!--template slot="row-details" slot-scope="row">
                    <b-card>
                        TODO add peek and other options
                    </b-card>
                </template-->
                <template slot="table-busy" class="text-center my-2">
                    <b-spinner class="align-middle"></b-spinner>
                    <strong>Loading...</strong>
                </template>
                <template slot="[]">Unexpected column</template>
            </b-table>
        </b-row>
    </b-container>
</template>

<script>
    import * as galaxy from '@/galaxy';
    import DatasetItem from './HistoryItems/Dataset';
    import CollectionItem from "./HistoryItems/Collection";
    import FunctionIcon from "@/galaxy/misc/FunctionIcon";

    export default {
        name: "HistoryContents",
        components: {
            FunctionIcon,
            CollectionItem,
            DatasetItem,
        },
        props: {
            historyPromise: {
                type: Promise,
                required: true,
            },
            permitted_file_extensions: {
                type: Array,
                default(){ return [] },
            },
            filter: {
                type: Function,
                default: (row, filter)=>{
                    if (typeof filter === "string")
                        return row.item.hid.toString().includes(filter) || row.item.name.includes(filter);
                    else if (filter.hasOwnProperty('test'))
                        return filter.test(row.item.hid.toString()) || filter.test(row.item.name);
                    return true;
                },
            },
            value: {
                type: Array,
                default() { return [] },
            },
        },
        data() {return{
            upload_dragging: false,
            user_filter: '',
        }},
        asyncComputed: {
            model() { return this.historyPromise },
        },
        computed: {
            items() {
                if (this.model === null || this.upload_dragging) return [];
                let history = galaxy.histories.History.query().with('datasets.history').find(this.model.id);
                //TODO when features available replace concat with v-for..of or model.morphTo(element property)
                return history.datasets
                    .concat(history.collections)
                    .filter(item=>item.deleted) //TODO make this optional, needs ui control
                    .sort((a,b)=>(a.hid === 0)?-1:b.hid-a.hid)
                    .reduce((acc, cur)=>{acc.push({item: cur}); return acc;}, []);
            },
            isLoading() {
                return this.model === null;
            }
        },
        methods: {
            show_upload() {
                this.$refs.upload.click();
            },
            uploadHandler(evt) {
                this.upload_dragging=false;
                if (this.model === null) return; //If user tries to upload before finished loading, do nothing.
                this.$emit('upload', {
                    history: this.model,
                    files: evt.dataTransfer ? evt.dataTransfer.files : evt.target.files,
                    permitted_file_extensions: this.permitted_file_extensions,
                    default() {this.model.fileUpload(evt.files);},
                });
            },
            clearSelected() { this.$refs.table.clearSelected(); this.$emit('input', this.items.map(item=>item.item))},
        },
        mounted() {
            if (this.value.length) {
                // Select any items passed in value prop
                for (const [index, item] of this.items.elements()) {
                    if (this.value.includes(item)) this.$refs.table.selectRow(index);
                }
            }
        },
    }
</script>

<style>
    .galaxy-history-contents .hidden_header {
        display: none;
    }

    .galaxy-history-contents .row:first-child {
        flex-wrap: nowrap;
        align-items: center;
    }

    .galaxy-history-contents .row:first-child .galaxy-function i {
        font-size: 1.5rem;
        padding-left: 0.5vw;
    }
</style>