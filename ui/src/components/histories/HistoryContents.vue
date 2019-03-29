<template>
    <div class="HistoryContents" @mousedown.stop.prevent @mousemove.stop.prevent>
        <div class="DropZone" @dragover.prevent="upload_dragging=true" @dragleave="upload_dragging=false" @dragexit="upload_dragging=false" @drop.prevent="uploadHandler"></div>
        <p v-if="items.length === 0 || upload_dragging">Drag and drop files here to upload</p>
        <ul v-else>
            <template v-for="(item, index) of items">
                <DatasetItem
                    v-if="item.history_content_type==='dataset'"
                    v-bind:key="index"
                    v-bind:model="item"
                    v-bind:class="{selected: selection[index]}"
                    @click.native.left.prevent="selectHandler($event, index)"
                    @history-item-deleted="update"
                />
                <CollectionItem
                    v-else-if="item.history_content_type==='dataset_collection'"
                    v-bind:key="index"
                    v-bind:model="item"
                    v-bind:class="{selected: selection[index]}"
                    @click.native.left.prevent="selectHandler($event, index)"
                    @history-item-deleted="update"
                />
            </template>
        </ul>
    </div>
</template>

<script>
    import * as galaxy from '@/galaxy';
    import DatasetItem from './HistoryItems/Dataset';
    import CollectionItem from "./HistoryItems/Collection";
    import {HistoryDatasetAssociation, HistoryDatasetCollectionAssociation} from "@/api/history_contents";

    export default {
        name: "HistoryContents",
        components: {
            CollectionItem,
            DatasetItem,
        },
        props: {
            model: {
                type: galaxy.histories.History,
                required: true,
            },
            filter: {
                type: String,
                default: '',
            }
        },
        data() {return{
            upload_dragging: false,
            upload_pending: [],
            selection: Array(this.$props.model.datasets.length + this.$props.model.collections.length).fill(false), //TODO move logic to item holding selection state
            last_selected: null,
        }},
        computed: {
            items() {
                //Changing sort order will require reworking this.selection grow/shrink logic
                return this.$props.model.datasets.concat(this.$props.model.collections, this.upload_pending).sort((a,b)=>(a.hid === 0)?-1:b.hid-a.hid); //TODO when features available replace with v-for..of or model.morphTo(element property)
            },
        },
        methods: {
            update() {
                this.upload_pending = this.upload_pending.filter(model=>model.file);
                this.model = History.query().with('datasets').with('collections').find(this.model.id);
            },
            uploadHandler(evt) {
                this.upload_dragging=false;
                this.upload_pending = this.upload_pending.filter(model=>model.file);
                let i, file;
                //let items = [];
                let self = this;
                if (evt.dataTransfer.items) {
                    // Use DataTransferItemList interface to access the file(s)
                    for (i = 0; i < evt.dataTransfer.items.length; i++) {
                        // If dropped items aren't files, reject them
                        if (evt.dataTransfer.items[i].kind === 'file') {
                            file = evt.dataTransfer.items[i].getAsFile();
                            let model = new galaxy.history_contents.HistoryDatasetAssociation();
                            model.name = file.name;
                            model.history_id = self.model.id;
                            model.file_name = file.name;
                            model.file = file;
                            self.upload_pending.push(model);
                        }
                    }
                } else {
                    // Use DataTransfer interface to access the file(s)
                    for (i = 0; i < evt.dataTransfer.files.length; i++) {
                        file = evt.dataTransfer.files[i];
                        let model = new galaxy.history_contents.HistoryDatasetAssociation();
                        model.name = file.name;
                        model.history_id = self.model.id;
                        model.file_name = file.name;
                        model.file = file;
                        self.upload_pending.push(model);
                    }
                }
                //Grow selection, this depends on the sort order of items being descending
                this.selection = Array(this.upload_pending.length).fill(false).concat(this.selection);
            },
            selectHandler(evt, index) {
                let start = Math.min(this.last_selected, index);
                let stop = Math.max(this.last_selected, index);
                let dist = stop - start;
                if (evt.ctrlKey && evt.shiftKey) {
                    this.selection.splice(start, dist, ...Array(dist).fill(true));
                    this.last_selected = index;
                } else if (evt.ctrlKey) {
                    this.last_selected = index;
                    this.$set(this.selection, this.last_selected, !this.selection[this.last_selected]);
                } else if (evt.shiftKey) {
                    this.selection.fill(false);
                    this.selection.splice(start, dist, ...Array(dist).fill(true));
                    this.last_selected = index;
                } else {
                    this.last_selected = index;
                    let selected = this.selection[this.last_selected];
                    this.selection.fill(false);
                    this.$set(this.selection, this.last_selected, !selected);
                }
            },
            getSelectedItems() {
                return this.items.filter((e,i)=>this.selection[i]);
            }
        },
        mounted() {
            let self = this;
            if (this.model.contents_url) {
                //Manually populate datasets and collections, in case this.model is not in store.
                galaxy.history_contents.HistoryDatasetAssociation.$fetch({
                    params: {
                        url: this.model.contents_url,
                    }
                }).then(()=>{
                    self.model.datasets = HistoryDatasetAssociation.query().where('history_id', self.model.id).get();
                });
                galaxy.history_contents.HistoryDatasetCollectionAssociation.$fetch({
                    params: {
                        url: this.model.contents_url,
                    }
                }).then(()=>{
                    self.model.collections = HistoryDatasetCollectionAssociation.query().where('history_id', self.model.id).get();
                });
            }
        },
    }
</script>

<style scoped>
    .HistoryContents {
        position: relative;
        min-height: 10em;
        height: 100% border-box;
        overflow-y: scroll;
    }

    .HistoryContents * {
        position: absolute;
        width: 100%;
        height: 100%;
        box-sizing: border-box;
    }

    .HistoryContents .DropZone {

    }

    .HistoryContents p {
        left: 50%;
        top: 50%;
        transform: translate(-50%, -50%);
        z-index: -1;
    }

    .HistoryContents ul * {
        position: relative;
        height: auto;
        border-bottom: grey 1px;
    }

    .HistoryContents .selected {
        background-color: azure;
    }

</style>