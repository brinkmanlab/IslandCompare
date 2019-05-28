<template>
    <div v-bind:class="name">
        <iframe v-bind:src="src"></iframe>
    </div>
</template>

<script>
    export default {
        name: "IFrameContent",
        props: {
            src: {
                type: String,
                required: true,
            },
            name: {
                type: String,
                required: true,
            }
        }
    }
</script>

<style scoped>

</style>