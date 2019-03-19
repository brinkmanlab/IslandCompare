<template>
    <HistoryItem ref="history_item" v-bind="$props">
        <span>{{state}}</span>
        <!--div class="HistoryItemAdv">
            <HistoryItemFunctions />
            <DatasetPreview />
        </div-->
    </HistoryItem>
</template>

<script>
    import axios from 'axios';
    import HistoryItem from "../HistoryItem";

    export default {
        extends: HistoryItem,
        name: "Dataset",
        components: {
            HistoryItem,
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
                    formData.append('key', this.$store.state.galaxy.api_key);
                    formData.append('history_id', this.$store.state.galaxy.current_history);
                    var inputs = {
                        'files_0|NAME'  : this.$props.file.name,
                        'files_0|type'  : 'upload_dataset',
                        'dbkey'         : '?',
                        'file_type'     : 'auto',
                        'ajax_upload'   : 'true',
                    };
                    formData.append('inputs', JSON.stringify(inputs));
                    formData.append('tool_id', 'upload1');
                    formData.append('files_0|file_data', this.$props.file, this.$props.file.name);
                    axios.post('/api/tools/', formData, {
                        headers: {'Content-Type': 'multipart/form-data'},
                        onUploadProgress: progressEvent => {
                            //console.log(progressEvent.total); // eslint-disable-line no-console
                            this.$refs.history_item.progress = progressEvent.loaded * 100 / progressEvent.total;
                        }
                    }).then(response => {
                        if (response.status == 200) {
                            Object.assign(this.$refs.history_item.$data, response.data.outputs[0]);
                        }
                        //TODO else
                    }).catch(error=>{
                        //TODO
                        this.$refs.history_item.$data.state = error.status;
                    })
                })
            }
        },
    }
</script>

<style scoped>

</style>