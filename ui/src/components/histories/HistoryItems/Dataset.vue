<template>
    <HistoryItem ref="history_item" v-bind="$props">
        <!--span>{{model.state}}</span>
        <div class="HistoryItemAdv">
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
        }},
        props: {
        },
        computed: {
        },
        methods: {
        },
        mounted() {
        },
    }
</script>

<style scoped>

</style>