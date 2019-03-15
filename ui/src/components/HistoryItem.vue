<template>
    <li>
        <span class="name">{{index}}: {{name}}</span>
        <progress max="100" v-bind:value="progress"></progress>
        <slot></slot>
    </li>
</template>

<script>
    export default {
        name: "HistoryItem",
        data: ()=>{return{
            progress: 0,
        }},
        computed: {
            name: ()=>{ return "History Item" },
            index: ()=>{ return 0 },
        }
    }
</script>

<style scoped>

</style>