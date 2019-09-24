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
                   v-bind:accept="accepted_upload_types.map(ext=>'.'+ext)"
                   @input.prevent="uploadHandler"
            />
        </b-row>

        <b-row>
            <b-table
                     show-empty
                     small
                     borderless
                     sticky-header="10rem"
                     per-page="0"
                     primary-key="key"
                     :items="items"
                     :busy="isLoading"
                     selectable
                     select-mode="range"
                     :fields="['item']"
                     :filter="user_filter"
                     :filter-function="filterFunc"
                     sort-by="hid"
                     empty-text="Drag and drop files here to upload"
                     thead-class="hidden_header"
                     @row-selected="onInput"
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
                <template slot="table-busy">
                    <div class="text-center">
                        <b-spinner class="align-middle"></b-spinner>
                        <strong>Loading...</strong>
                    </div>
                </template>
                <template slot="[]">Unexpected column</template>
            </b-table>
        </b-row>
    </b-container>
</template>

<script>
    import { History } from "@/galaxy/api/histories";
    import { HistoryDatasetAssociation } from "@/galaxy/api/history_contents";
    import DatasetItem from './HistoryItems/Dataset';
    import CollectionItem from "./HistoryItems/Collection";
    import FunctionIcon from "@/galaxy/misc/FunctionIcon";

    const temporary_extension_to_datatype_map = {
        "genbank": "genbank", "gbk": "genbank", "embl": "embl", "gbff": "genbank"
    };

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
            filter: {
                type: Function,
                default: ()=>true,
            },
            value: {
                type: Array,
                default() { return [] },
            },
            accepted_upload_types: {
                type: Array,
                default() { return [] },
            }
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
                //TODO when features available replace concat with v-for..of or model.morphTo(element property)
                const history = History.query().with('datasets.history').find(this.model.id); // TODO the reactivity system fails to update if .with(datasets) in historyPromise
                return history.datasets
                    .concat(history.collections)
                    .filter(item=>!item.deleted && this.filter(item)) //TODO make deleted optional, needs ui control
                    .sort((a,b)=>(a.hid === 0)?-1:b.hid-a.hid)
                    .reduce((acc, cur)=>{acc.push({item: cur, hid: cur.hid, key: cur.id}); return acc;}, []);
            },
            isLoading() {
                return this.model === null;
            }
        },
        methods: {
            filterFunc(row, filter) {
                if (typeof filter === "string")
                    return row.item.hid.toString().includes(filter) || row.item.name.includes(filter);
                else if (filter.hasOwnProperty('test'))
                    return filter.test(row.item.hid.toString()) || filter.test(row.item.name);
                return true;
            },
            show_upload() {
                this.$refs.upload.click();
            },
            onInput(items) {
                this.$emit('input', items.map(item=>item.item));
            },
            uploadHandler(evt) {
                this.upload_dragging=false;
                if (this.model === null) return; //If user tries to upload before finished loading, do nothing.
                const files = evt.dataTransfer ? evt.dataTransfer.files : evt.target.files;

                // Enforce accepted file types
                const accepted_files = [];
                for (const file of files) {
                    let ext = file.name.match(/[^.]+$/);
                    console.log(ext);
                    if (this.accepted_upload_types.length && !(ext && this.accepted_upload_types.includes(ext[0]))) {
                        let tmp_id = file.name + Math.floor(Math.random() * 10 ** 16).toString();
                        HistoryDatasetAssociation.insert({
                            data: {
                                id: tmp_id,
                                file: file,
                                name: "Incorrect file format: " + file.name,
                                hid: -1,
                                history_id: this.model.id,
                                extension: this.accepted_upload_types[0],
                            }
                        });
                    } else {
                        accepted_files.push(file);
                    }
                }

                // Hand off upload event to parent for further processing
                const self = this;
                this.$emit('upload', {
                    history: this.model,
                    files: accepted_files,
                    default() { for (const file of accepted_files) {
                        // TODO handle specifying extension elsewhere
                        // TODO hardcoded extension mapping to avoid sniff for now
                        let ext = file.name.match(/[^.]+$/);
                        /*if (ext && temporary_extension_to_datatype_map.hasOwnProperty(ext[0])) self.model.fileUpload(file, temporary_extension_to_datatype_map[ext[0]]);
                        else*/ self.model.fileUpload(file);
                    }},
                });
            },
            clearSelected() { this.$refs.table.clearSelected(); /*this.$emit('input', this.items.map(item=>item.item))*/},
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

    .galaxy-history-contents .b-table-sticky-header {
        width: 100%;
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