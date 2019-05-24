<template>
    <div class="HistoryContents" @mousedown.stop.prevent @mousemove.stop.prevent @dragover.prevent="upload_dragging=true" @dragleave="upload_dragging=false" @dragexit="upload_dragging=false" @drop.prevent="uploadHandler">
        <div v-if="model === null">Loading data</div>
        <div v-else-if="all_deleted() || upload_dragging">Drag and drop files here to upload</div>
        <ul v-else>
            <template v-for="(item, index) of items">
                <DatasetItem
                    v-if="item.history_content_type==='dataset' && !item.deleted"
                    v-bind:key="index"
                    v-bind:model="item"
                    v-bind:class="{selected: selection[index]}"
                    @click.native.left.prevent="selectHandler($event, index)"
                />
                <CollectionItem
                    v-else-if="item.history_content_type==='dataset_collection' && !item.deleted"
                    v-bind:key="index"
                    v-bind:model="item"
                    v-bind:class="{selected: selection[index]}"
                    @click.native.left.prevent="selectHandler($event, index)"
                />
            </template>
        </ul>
    </div>
</template>

<script>
    import * as galaxy from '@/galaxy';
    import DatasetItem from './HistoryItems/Dataset';
    import CollectionItem from "./HistoryItems/Collection";
    export default {
        name: "HistoryContents",
        components: {
            CollectionItem,
            DatasetItem,
        },
        props: {
            model: [ galaxy.histories.History, null],
            filter: {
                type: String,
                default: '',
            }
        },
        data() {return{
            upload_dragging: false,
            last_selected: null,
            selected: new Map(),
        }},
        computed: {
            items() {
                if (this.model === null) return;
                let history = galaxy.histories.History.query().with('datasets.history').find(this.model.id);
                //TODO when features available replace with v-for..of or model.morphTo(element property)
                return history.datasets
                    .concat(history.collections)
                    .sort((a,b)=>(a.hid === 0)?-1:b.hid-a.hid);
            },
        },
        methods: {
            uploadHandler(evt) {
                this.upload_dragging=false;
                let i, file;
                //let items = [];
                let self = this;
                if (evt.dataTransfer.items) {
                    // Use DataTransferItemList interface to access the file(s)
                    for (i = 0; i < evt.dataTransfer.items.length; i++) {
                        // If dropped items aren't files, reject them
                        if (evt.dataTransfer.items[i].kind === 'file') {
                            file = evt.dataTransfer.items[i].getAsFile();
                            galaxy.history_contents.HistoryDatasetAssociation.$upload(file, self.model.id);
                        }
                    }
                } else {
                    // Use DataTransfer interface to access the file(s)
                    for (i = 0; i < evt.dataTransfer.files.length; i++) {
                        file = evt.dataTransfer.files[i];
                        galaxy.history_contents.HistoryDatasetAssociation.$upload(file, self.model.id);
                    }
                }
            },
            selectHandler(evt, index) {
                let start = Math.min(this.last_selected, index);
                let stop = Math.max(this.last_selected, index);
                let dist = stop - start + 1;
                let ctrlKey = evt.ctrlKey;
                if (navigator.appVersion.indexOf("Mac") !== -1) ctrlKey = evt.metaKey; // Handle Macs screwey key layout
                if (ctrlKey && evt.shiftKey) {
                    this.selection.splice(start, dist, ...Array(dist).fill(true));
                    this.last_selected = index;
                } else if (ctrlKey) {
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
        },
    }
</script>

<style scoped>
    .HistoryContents {
        position: relative;
        min-height: 10em;
        height: 100% border-box;
        overflow-y: scroll;
        user-select: none;
    }

    .HistoryContents > * {
        position: absolute;
        box-sizing: border-box;
    }

    .HistoryContents > div:first-child {
        left: 50%;
        top: 50%;
        transform: translate(-50%, -50%);
        z-index: -1;
    }

    .HistoryContents > ul {
        width: 100%;
        height: 100%;
        margin: 0;
        padding: 0.5em;
    }

    .HistoryContents > ul > * {
        position: relative;
        height: auto;
    }

    .HistoryContents .selected {
        background-color: azure;
    }

</style>