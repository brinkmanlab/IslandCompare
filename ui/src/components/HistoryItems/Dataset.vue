<template>
    <HistoryItem>
        <span>{{state}}</span>
        <HistoryItemFunctions>
            <slot name="history_item_functions"></slot>
            <RemoveHistoryItem />
        </HistoryItemFunctions>
        <!--div class="HistoryItemAdv">
            <HistoryItemFunctions />
            <DatasetPreview />
        </div-->
    </HistoryItem>
</template>

<script>
    import axios from 'axios';
    import HistoryItem from "../HistoryItem";
    import HistoryItemFunctions from "../HistoryItemFunctions";
    import RemoveHistoryItem from "../HistoryItemFunctions/Remove";

    export default {
        extends: HistoryItem,
        name: "Dataset",
        components: {
            HistoryItem,
            HistoryItemFunctions,
            RemoveHistoryItem,
            //DatasetPreview,
        },
        data: ()=>{return{
            state: '',
        }},
        props: {
            file: {
                type: File,
                default: null,
            },
        },
        computed: {
        },
        methods: {
        },
        mounted() {
            if (this.file) {
                this.$nextTick(() => {
                    var formData = new FormData();
                    var key = this.$vnode.key;
                    formData.append('key', this.$store.state.galaxy.api_key);
                    formData.append('history_id', this.$store.state.galaxy.current_history);
                    formData.append('inputs', '{"dbkey":"?","file_type":"auto","' + key + '|type":"upload_dataset","' + key + '|space_to_tab":null,"' + key + '|to_posix_lines":"Yes","' + key + '|NAME":"' + this.$props.file.name + '"}');
                    formData.append('tool_id', 'upload1');
                    formData.append(key + '|file_data', this.$props.file, this.$props.file.name);
                    axios.post('http://galaxy.brinkman.mbb.sfu.ca:8000/api/uploads/', formData, {
                        headers: {'Content-Type': 'multipart/form-data'},
                        onUploadProgress: progressEvent => this.$data.progress = progressEvent.loaded
                    }).then(response => {
                        if (response.status == 200) {
                            console.log(response.statusText); // eslint-disable-line no-console
                        }
                    }).catch((error)=>{
                        this.$data.state = error.status;
                    })
                })
            }
        },
    }
</script>

<style scoped>

</style>