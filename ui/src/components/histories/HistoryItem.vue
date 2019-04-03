<template>
    <li class="HistoryItem">
        <slot name="before"/>
        <span class="name">{{ model.hid }}: {{ model.name }}</span>
        <progress max="100" v-bind:value="model.progress" v-if="model.progress<100"></progress>
        <slot></slot>
        <HistoryItemFunctions v-bind:item="this">
            <template v-slot:default="slot">
                <slot name="history_item_functions" v-bind:item="slot.item"></slot>
                <RemoveHistoryItem v-bind:item="slot.item" />
            </template>
        </HistoryItemFunctions>
        <slot name="after"/>
    </li>
</template>

<script>
    import * as galaxy from '@/galaxy'
    import HistoryItemFunctions from "./HistoryItemFunctions";
    import RemoveHistoryItem from "./HistoryItemFunctions/Remove";

    export default {
        name: "HistoryItem",
        components: {
            HistoryItemFunctions,
            RemoveHistoryItem,
        },
        data: ()=>{return{
        }},
        props: {
            model: {
                type: [galaxy.history_contents.HistoryDatasetAssociation, galaxy.history_contents.HistoryDatasetCollectionAssociation],
                required: true,
            },
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
.HistoryItem {
    display: flex;
    flex-direction: row;
    align-items: stretch;
    width: 100%;
}
</style>