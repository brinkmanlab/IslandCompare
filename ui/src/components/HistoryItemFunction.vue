<template>
    <button class="HistoryItemFunction" v-bind:title="description">
        <i v-bind:class="icon"></i>
        <label>{{label}}</label>
    </button>
</template>

<script>
    export default {
        name: "HistoryItemFunction",
        props: {
            icon: {
                type: String,
                default: '',
            },
            label: {
                type: String,
                default: '',
            },
            description: {
                type: String,
                default: '',
            },
        }
    }
</script>

<style scoped>

</style>