<template>
    <b-button-group class="galaxy-history-item-functions">
        <slot v-bind:item="item"></slot>
    </b-button-group>
</template>

<script>
    export default {
        name: "HistoryItemFunctions",
        props: {
            item: {
                type: Object,
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