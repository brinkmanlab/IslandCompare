<template>
    <b-container class="Home">
        <b-row>
            <b-col class="floating welcome">
                <HTMLFragment :content="welcome" />
            </b-col>
        </b-row>
       <b-row>
            <b-col class="floating examples">
                <h1>Getting Started</h1>
                <HTMLFragment :content="examples" />
            </b-col>
            <b-col class="floating">
                <h1>News and Updates</h1>
                <News :news="news"/>
                <b-link to="/news">View more...</b-link>
            </b-col>
        </b-row>
    </b-container>
</template>

<script>
    import HTMLFragment from "@/components/HTMLFragment";
    import News from "@/components/News";

    import news from "@/assets/news";
    import welcome from "html-loader!@/assets/welcome.htm"
    import examples from "html-loader!@/assets/examples.htm"

    export default {
        name: "Home",
        components: {News, HTMLFragment},
        data() {return{
            origin,
            news,
            welcome,
            examples,
        }},
    }
</script>

<style scoped>
    .col, .col-8 {
        margin-bottom: 15px;
    }

    .col.floating, .col-8.floating {
        margin-left: 15px;
        margin-right: 15px;
        padding-top: 0.5em;
        padding-bottom: 0.5em;
    }

    .col > h1, .col-6 > h1 {
        border-bottom: 1px solid black;
    }

    /*.welcome {
        background: #d4dcef;
        color: black;
    }*/

    >>> .buttons {
        text-align: center;
    }

    .welcome >>> .static_content > img {
        /*max-width: 15em;*/
        max-height: 12em;
        float: right;
        margin-left: 1rem;
    }

    .welcome >>> .static_content > :last-child {
        margin-bottom: unset;
    }

    .News {
        flex-flow: column;
        overflow-y: hidden;
        max-height: 70vh;
    }

    >>> .row {
        align-items: flex-start;
    }

    .News >>> .card {
        margin-bottom: 15px;
        min-width: 15em;
        border: 0;
    }

    .News >>> .card-header {
        padding: 0.5em;
    }

    .News >>> .card-body {
        padding: 0.5em;
    }

    .News + a {
        display: block;
        border-radius: 0.25rem;
        /*background: #d4dcef;*/
        text-align: center;
        margin-left: -2px;
        margin-right: -2px;
        color: grey;
        text-decoration: none;
        background: rgba(0, 0, 0, 0.03);
    }
</style>