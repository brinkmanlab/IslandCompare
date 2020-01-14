<template>
    <div v-html="content" class="static_content" :class="slug"></div>
</template>

<script>
    export default
    {
        name: 'HTMLFragment',
        components: {
        },
        props: ["content", "slug"],
        // TODO https://github.com/vuejs/vue-router/issues/48
    }
</script>