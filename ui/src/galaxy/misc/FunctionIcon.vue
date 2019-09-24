<template>
    <b-link class="galaxy-function" v-bind:title="description" @click.stop.prevent="$props.action" href="">
        <i class="galaxy-function-icon" v-bind:class="icon"></i>
        {{label}}
    </b-link>
</template>

<script>
    export default {
        name: "FunctionIcon",
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
            action: {
                type: Function,
                required: true,
            }
        },
    }
</script>

<style scoped>

</style>