<template>
    <div class="HistoryContents">
        <div class="DropZone" @dragover.prevent="upload_dragging=true" @dragleave="upload_dragging=false" @dragexit="upload_dragging=false" @drop.prevent="dropHandler"></div>
        <p v-if="items.length === 0 || upload_dragging">Drag and drop files here to upload</p>
        <ul v-else>
            <template v-for="(item, index) of items">
                <DatasetItem
                    v-if="item.history_content_type==='dataset'"
                    v-bind:key="item.id + '_' + item.history_content_type"
                    v-bind:id="item.id"
                    v-bind:file="item.file || null"
                    v-bind:initial_data="item"
                    @history-item-deleted="items.splice(index)"
                />
                <CollectionItem
                    v-else-if="item.history_content_type==='dataset_collection'"
                    v-bind:key="item.id + '_' + item.history_content_type"
                    v-bind:id="item.id"
                    v-bind:initial_data="item"
                    @history-item-deleted="items.splice(index)"
                />
            </template>
        </ul>
    </div>
</template>

<script>
    import * as Model from '@/api/history_contents';
    import { Model as HistoryModel } from '@/api/histories';
    import DatasetItem from './HistoryItems/Dataset';
    import CollectionItem from "./HistoryItems/Collection";

    export default {
        name: "HistoryContents",
        components: {
            CollectionItem,
            DatasetItem
        },
        props: {
            model: {
                type: HistoryModel,
                required: true,
            },
            filter: {
                type: String,
                default: '',
            }
        },
        data: ()=>{return{
            upload_dragging: false,
        }},
        computed: {
            items() {
                return this.$props.model.datasets.concat(this.$props.model.collections).sort((a,b)=>a.hid-b.hid); //TODO when features available replace with v-for..of or model.morphTo(element property)
            }
        },
        methods: {
            dropHandler(evt) {
                this.upload_dragging=false;
                var i, file;
                var items = [];
                if (evt.dataTransfer.items) {
                    // Use DataTransferItemList interface to access the file(s)
                    for (i = 0; i < evt.dataTransfer.items.length; i++) {
                        // If dropped items aren't files, reject them
                        if (evt.dataTransfer.items[i].kind === 'file') {
                            file = evt.dataTransfer.items[i].getAsFile();
                            items.push({type: 'dataset', id: '', file: file});
                        }
                    }
                } else {
                    // Use DataTransfer interface to access the file(s)
                    for (i = 0; i < evt.dataTransfer.files.length; i++) {
                        file = evt.dataTransfer.files[i];
                        items.push({type: 'dataset', id: '', file: file});
                    }
                }
                this.items = [...items, ...this.items];
            },
        },
        mounted() {
            if (this.model.contents_url) {
                Model.HistoryDatasetAssociation.$fetch({
                    params: {
                        contents_url: this.model.contents_url,
                    }
                });
                Model.HistoryDatasetCollectionAssociation.$fetch({
                    params: {
                        contents_url: this.model.contents_url,
                    }
                });
            }
        },
    }
</script>

<style scoped>
    .HistoryContents {
        position: relative;
        height: 500px;
        overflow-y: scroll;
    }

    .HistoryContents * {
        position: absolute;
    }

    .HistoryContents .DropZone {
        width: 100%;
        Height: 100%;
    }

    .HistoryContents p {
        left: 50%;
        top: 50%;
        transform: translate(-50%, -50%);
        z-index: -1;
    }

    .HistoryContents input[type=file] {
        /*opacity: 0;*/
        z-index: 1;
    }

    .HistoryContents ul * {
        position: relative;
    }
</style>