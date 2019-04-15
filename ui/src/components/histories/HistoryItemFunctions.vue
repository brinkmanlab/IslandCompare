<template>
    <div class="HistoryItemFunctions">
        <slot v-bind:item="item"></slot>
    </div>
</template>

<script>
    export default {
        name: "HistoryItemFunctions",
        props: {
            item: {
                type: Object,
                required: true,
            }
        }
    }
</script>

<style scoped>
    .HistoryItemFunctions {
        display: inline-flex;
        flex-direction: row;
        align-items: stretch;
        margin-left: auto;
    }
</style>