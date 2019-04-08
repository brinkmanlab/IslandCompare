<template>
    <div class="HistoryContents" @mousedown.stop.prevent @mousemove.stop.prevent @dragover.prevent="upload_dragging=true" @dragleave="upload_dragging=false" @dragexit="upload_dragging=false" @drop.prevent="uploadHandler">
        <p v-if="items.length === 0 || upload_dragging">Drag and drop files here to upload</p>
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
            selection: Array(this.model.datasets.length + this.model.collections.length).fill(false), //TODO move logic to item holding selection state
            last_selected: null,
        }},
        computed: {
            items() {
                let history = galaxy.histories.History.query().with('datasets.history').find(this.model.id);
                //Changing sort order will require reworking this.selection grow/shrink logic
                return history.datasets.concat(history.collections).sort((a,b)=>(a.hid === 0)?-1:b.hid-a.hid); //TODO when features available replace with v-for..of or model.morphTo(element property)
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
        width: 100%;
        height: 100%;
        box-sizing: border-box;
    }

    .HistoryContents > p {
        left: 50%;
        top: 50%;
        transform: translate(-50%, -50%);
        z-index: -1;
    }

    .HistoryContents > ul * {
        position: relative;
        height: auto;
        border-bottom: grey 1px;
    }

    .HistoryContents .selected {
        background-color: azure;
    }

</style>