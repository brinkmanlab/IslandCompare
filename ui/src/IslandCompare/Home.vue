<template>
    <b-container class="Home">
        <b-row>
            <b-col class="floating welcome">
                <HTMLFragment :content="welcome" />
            </b-col>
        </b-row>
        <b-row>
            <b-col class="floating examples" cols="8">
                <HTMLFragment :content="examples" />
            </b-col>
            <b-col>
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
    //align-self: center;
    }

    .col.floating, .col-8.floating {
        margin-left: 15px;
        margin-right: 15px;
        padding-top: 0.5em;
        padding-bottom: 0.5em;
    }

    .welcome {
    //min-width: 40vw;
        background: #d4dcef;
        color: black;
    }

    >>> .buttons {
        text-align: center;
    }

    >>> .buttons > a {
        background: #7bbdf7;
        color: white;
        height: 5em;
        width: 10em;
        border-radius: 10px;
        margin: 10px;
        text-align: center;
        display: inline-block;
        padding: 4px;
        text-decoration: none;
    }

    >>> .buttons .icon {
        font-size: 3em;
        display: block;
        margin-top: 2px;
    }

    .welcome >>> .static_content > img:first-child {
        //max-width: 15em;
        max-height: 11em;
        float: left;
        margin-right: 1em;
    }

    .welcome >>> .static_content > :last-child {
        margin-bottom: unset;
    }

    .News {
        //flex-direction: column;
        overflow-y: hidden;
        max-height: 50vh;
        margin-top: -10px;
        padding-top: 10px;
    }

    >>> .row {
        align-items: flex-start;
    }

    .News >>> .card {
        margin-bottom: 15px;
        min-width: 15em;
    }

    .News + a {
        display: block;
        border-radius: 0.25rem;
        background: #d4dcef;
        text-align: center;
        margin-left: -2px;
        margin-right: -2px;
        color: grey;
        text-decoration: none;
    }
</style>