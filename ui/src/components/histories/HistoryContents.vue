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
                    v-bind:class="{'table-primary': selected.get(item.id) === true, 'table-danger': item.hid === -1}"
                    @click.native.left.prevent="selectHandler($event, index)"
                />
                <CollectionItem
                    v-else-if="item.history_content_type==='dataset_collection' && !item.deleted"
                    v-bind:key="index"
                    v-bind:model="item"
                    v-bind:class="{'table-primary': selected.get(item.id) === true}"
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
            },
            upload_callback: {
                type: Function,
                default: u=>u,
            },
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
                if (this.model === null) return; //If user tries to upload before finished loading, do nothing.
                let i, file;
                //let items = [];
                let self = this;
                if (evt.dataTransfer.items) {
                    // Use DataTransferItemList interface to access the file(s)
                    for (i = 0; i < evt.dataTransfer.items.length; i++) {
                        // If dropped items aren't files, reject them
                        if (evt.dataTransfer.items[i].kind === 'file') {
                            file = evt.dataTransfer.items[i].getAsFile();
                            file = this.upload_callback(file);
                            if (file)
                                galaxy.history_contents.HistoryDatasetAssociation.$upload(file, self.model.id);
                        }
                    }
                } else {
                    // Use DataTransfer interface to access the file(s)
                    for (i = 0; i < evt.dataTransfer.files.length; i++) {
                        file = evt.dataTransfer.files[i];
                        file = this.upload_callback(file);
                        if (file)
                            galaxy.history_contents.HistoryDatasetAssociation.$upload(file, self.model.id);
                    }
                }
            },
            selectHandler(evt, index) {
                let start = Math.min(this.last_selected, index);
                let stop = Math.max(this.last_selected, index);
                this.last_selected = index;
                let ctrlKey = evt.ctrlKey;
                if (navigator.appVersion.indexOf("Mac") !== -1) ctrlKey = evt.metaKey; // Handle Macs screwey key layout
                for (let [ind, item] of this.items.entries()) {
                    if (item.hid <= 0) continue; //Prevent selecting ghost items
                    if (ctrlKey && evt.shiftKey) {
                        this.selected.set(item.id, (this.selected.get(item.id) === true) || (start <= ind && ind <= stop));
                    } else if (ctrlKey) {
                        if (index === ind) this.selected.set(item.id, !(this.selected.get(item.id) === true));
                    } else if (evt.shiftKey) {
                        this.selected.set(item.id, (start <= ind && ind <= stop));
                    } else {
                        this.selected.set(item.id, (index === ind));
                    }
                }
                this.$forceUpdate();
            },
            getSelectedItems() {
                return this.items.filter(i=>this.selected.get(i.id)===true && i.hid !== 0);
            },
            clearSelectedItems() {
                this.selected = new Map();
            },
            all_deleted() {
                // Iteratables dont have every() function. This is a stand-in.
                for (let item of this.items.keys()) if (!item.deleted) return false;
                return true;
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

    .HistoryContents > div:first-child:nth-child(-n+1) {
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

</style>