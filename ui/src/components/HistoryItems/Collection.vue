<template>
    <HistoryItem ref="history_item" v-bind="$props">
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
    import HistoryItem from "../HistoryItem";
    import HistoryItemFunctions from "../HistoryItemFunctions";
    import RemoveHistoryItem from "../HistoryItemFunctions/Remove";

    export default {
        extends: HistoryItem,
        name: "Collection",
        components: {
            HistoryItem,
            HistoryItemFunctions,
            RemoveHistoryItem,
        },
    }
</script>

<style scoped>

</style>