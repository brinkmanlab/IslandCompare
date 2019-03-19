<template>
    <div class="HistoryContents">
        <div class="DropZone" @dragover.prevent="upload_dragging=true" @dragleave="upload_dragging=false" @dragexit="upload_dragging=false" @drop.prevent="dropHandler"></div>
        <p v-if="items.length === 0 || upload_dragging">Drag and drop files here to upload</p>
        <ul v-else>
            <template v-for="(item, index) in items">
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
    import axios from 'axios';
    import DatasetItem from './HistoryItems/Dataset';
    import CollectionItem from "./HistoryItems/Collection";
    export default {
        name: "HistoryContents",
        components: {
            CollectionItem,
            DatasetItem
        },
        props: {
            history_id: {
                type: String,
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
            items: {
                get: ()=>this.$store.state.history_contents[this.$props.history_id],
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
            refresh() {
                axios.get(`/api/histories/${this.history_id}/contents/?key=${this.$store.state.galaxy.api_key}&type=${this.filter}`)
                    .then(response => {
                        this.$data.items = response.data;
                    })
                    .catch(error => {
                        //TODO
                        console.log(error.status); // eslint-disable-line no-console
                    });
            }
        },
        mounted() {
            this.refresh();
        }
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