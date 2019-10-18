<template>
    <b-card-group deck class="News">
        <b-card v-for="(item, index) of news" v-bind:key="index" v-bind:header="(new Date(item.date)).toDateString()" v-bind:title="item.title"><b-card-text v-html="item.body"></b-card-text></b-card>
    </b-card-group>
</template>

<script>
    export default {
        name: "News",
        props: {
            news: {
                type: Array,
                required: true,
            }
        }
    }
</script>

<style scoped>

</style>