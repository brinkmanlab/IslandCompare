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
    //import * as galaxy from '@/galaxy';
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
        },
        computed: {
        },
        methods: {
        },
        mounted() {
            if (this.model.file) {
                this.$nextTick(() => {
                    this.model.upload(this.model.file, {
                        onUploadProgress: progressEvent => {
                            //console.log(progressEvent.total); // eslint-disable-line no-console
                            this.$refs.history_item.progress = progressEvent.loaded * 100 / progressEvent.total;
                        }
                    }).then(model=>{
                        this.model = model;
                        this.$refs.history_item.progress = 0;
                    });
                    delete this.model.file;
                })
            }
        },
    }
</script>

<style scoped>

</style>