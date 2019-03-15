<template>
    <div class="HistoryList">
        <div class="DropZone" @dragover.prevent="upload_dragging=true" @dragleave="upload_dragging=false" @dragexit="upload_dragging=false" @drop.prevent="dropHandler"></div>
        <p v-if="items.length === 0 || upload_dragging">Drag and drop files here to upload</p>
        <ul v-else>
            <template v-for="(item, i) in items">
                <DatasetItem
                    v-if="item.type==='dataset'"
                    v-bind:key="i"
                    v-bind:id="item.id"
                    v-bind:file="item.file || null"
                />
                <!--CollectionItem v-else-if="item.type==='collection'"/-->
            </template>
        </ul>
    </div>
</template>

<script>
    import DatasetItem from './HistoryItems/Dataset';
    export default {
        name: "HistoryList",
        components: {
            DatasetItem
        },
        props: {
            history_id: {
                type: String,
                default: '',
                required: true,
            }
        },
        data: ()=>{return{
            upload_dragging: false,
            items: [],
        }},
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
                            console.log('... file[' + i + '].name = ' + file.name); // eslint-disable-line no-console
                        }
                    }
                } else {
                    // Use DataTransfer interface to access the file(s)
                    for (i = 0; i < evt.dataTransfer.files.length; i++) {
                        file = evt.dataTransfer.files[i];
                        items.push({type: 'dataset', id: '', file: file});
                        console.log('... file[' + i + '].name = ' + file.name); // eslint-disable-line no-console
                    }
                }
                this.items = [...items, ...this.items];
            },
        }
    }
</script>

<style scoped>
    .HistoryList {
        position: relative;
        height: 500px;
        width: 50%;
    }

    .HistoryList * {
        position: absolute;
    }

    .HistoryList .DropZone {
        width: 100%;
        Height: 100%;
    }

    .HistoryList p {
        left: 50%;
        top: 50%;
        transform: translate(-50%, -50%);
        z-index: -1;
    }

    .HistoryList input[type=file] {
        /*opacity: 0;*/
        z-index: 1;
    }

    .HistoryList ul * {
        position: relative;
    }
</style>