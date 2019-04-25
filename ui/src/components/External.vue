<template>
    <iframe class="External"></iframe>
</template>

<script>
    export default {
        name: "External",
        props: {
            path: {
                type: String,
                required: true,
            }
        }
    }
</script>

<style scoped>
    .External {
        border: 0;
    }
</style>