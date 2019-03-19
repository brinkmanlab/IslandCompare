<template>
    <div>
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

</style>