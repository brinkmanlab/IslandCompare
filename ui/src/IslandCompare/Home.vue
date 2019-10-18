<template>
    <b-container class="Home">
        <b-row>
            <b-col>
                <Slides :slides="slides" :loop="false"/>
                <HTMLFragment :content="welcome" />
            </b-col>
        </b-row>
        <b-row>
            <b-col>
                <News :news="news"/>
            </b-col>
        </b-row>
    </b-container>
</template>

<script>
    import HTMLFragment from "@/components/HTMLFragment";
    import News from "@/components/News";
    import Slides from "@/components/Slides";

    import news from "@/assets/news";
    import slides from "@/assets/slides";
    import welcome from "html-loader!@/assets/home.htm"

    slides.forEach(item=>{if (item.img) item.img = require('@/../static/slides/images/' + item.img)}); // Webpack will slurp up the images

    export default {
        name: "Home",
        components: {Slides, News, HTMLFragment},
        data() {return{
            origin,
            news,
            slides,
            welcome,
        }},
    }
</script>

<style scoped>
    .Slides >>> .carousel-item {
        height: 40vh;
        font-size: 0.8em;
    }
</style>