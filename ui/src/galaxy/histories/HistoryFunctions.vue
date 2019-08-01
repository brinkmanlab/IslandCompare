<template>
    <b-button-group class="galaxy-history-functions">
        <slot v-bind:item="item"></slot>
    </b-button-group>
</template>

<script>
    export default {
        name: "HistoryFunctions",
        props: {
            item: {
                type: Object,
                required: true,
            },

        }
    }
</script>

<style scoped>
    .galaxy-history-functions {
        display: inline-flex;
        flex-direction: row;
        align-items: stretch;
        margin-left: auto;
    }
</style>