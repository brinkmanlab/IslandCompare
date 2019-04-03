<template>
    <button class="HistoryItemFunction" v-bind:title="description" @click.stop="onClick">
        <i v-bind:class="icon"></i>
        <label>{{label}}</label>
    </button>
</template>

<script>
    export default {
        name: "HistoryItemFunction",
        props: {
            item: {
                type: Object,
                required: true,
            },
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
            success: {
                type: Function,
                default: response=>this.$emit('operation-success', response),
            },
            fail: {
                type: Function,
                default: error=>this.$emit('operation-fail', error),
            },
        },
        methods: {
        }
    }
</script>

<style scoped>

</style>