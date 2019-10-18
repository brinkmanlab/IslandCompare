<template>
    <div v-html="content" class="static_content"></div>
</template>

<script>
    import 'markdown-it-vue/dist/markdown-it-vue.css'
    export default
    {
        name: 'HTMLFragment',
        components: {
        },
        props: ["content"],
        // TODO https://github.com/vuejs/vue-router/issues/48
    }
</script>