<template>
    <div v-html="content"></div>
</template>

<script>
    import axios from 'axios/index'
    import 'markdown-it-vue/dist/markdown-it-vue.css'
    export default
    {
        name: 'HTMLFragment',
        components: {
        },
        props: ["src"],
        asyncComputed: {
            content: {
                get() {
                    return axios.get(this.src).then(response => response.data)
                },
                default: ""
            }
        },
    }
</script>