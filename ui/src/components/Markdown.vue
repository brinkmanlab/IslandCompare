<template>
    <markdown-it-vue class="md-body" :content="content"/>
</template>

<script>
    import axios from 'axios/index'
    import MarkdownItVue from 'markdown-it-vue'
    import 'markdown-it-vue/dist/markdown-it-vue.css'
    export default
    {
        name: 'Markdown',
        components: {
            MarkdownItVue,
        },
        props: ["src"],
        asyncComputed: {
            content: {
                get() {
                    return axios.get(this.src).then(response => response.data)
                },
                default: "# Loading.."
            }
        },
    }
</script>