<template>
    <b-button-group class="galaxy-history-item-functions">
        <slot v-bind:item="item"></slot>
        <RenameHistoryItem v-bind:item="item" v-on="$listeners"/>
        <RemoveHistoryItem v-bind:item="item" v-on="$listeners"/>
    </b-button-group>
</template>

<script>
    import RemoveHistoryItem from "./HistoryItemFunctions/Remove";
    import RenameHistoryItem from "./HistoryItemFunctions/Rename";
    import HistoryItem from "@/galaxy/histories/HistoryItem";

    export default {
        name: "HistoryItemFunctions",
        components: {RemoveHistoryItem, RenameHistoryItem},
        props: {
            item: {
                type: HistoryItem,
                required: true,
            },
        }
    }
</script>

<style scoped>
    .galaxy-history-item-functions {
        display: inline-flex;
        flex-direction: row;
        align-items: stretch;
        margin-left: auto;
    }
</style>