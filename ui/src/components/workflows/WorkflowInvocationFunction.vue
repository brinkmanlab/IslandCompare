<template>
    <b-link class="galaxy-workflow-invocation-function" v-bind:title="description" @click.stop.prevent="onClick" href="">
        <i class="galaxy-workflow-invocation-function-icon" v-bind:class="icon"></i>
        {{label}}
    </b-link>
</template>

<script>
    export default {
        name: "WorkflowInvocationFunction",
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
                default(){ return response=>this.$emit('operation-success', response) },
            },
            fail: {
                type: Function,
                default(){ return error=>this.$emit('operation-fail', error) },
            },

        },
        methods: {
        }
    }
</script>

<style scoped>

</style>