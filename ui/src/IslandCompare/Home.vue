<template>
    <b-container class="Home">
        <b-row>
            <b-col class="floating cite">
                For any publications resulting in the use of our services see <b-link to="/publications">how to cite this service</b-link>.
            </b-col>
        </b-row>
        <b-row>
            <b-col class="floating welcome">
                <HTMLFragment :content="welcome" />
            </b-col>
            <b-col>
                <News :news="news"/>
            </b-col>
        </b-row>
    </b-container>
</template>

<script>
    import HTMLFragment from "@/components/HTMLFragment";
    import News from "@/components/News";

    import news from "@/assets/news";
    import welcome from "html-loader!@/assets/home.htm"

    export default {
        name: "Home",
        components: {News, HTMLFragment},
        data() {return{
            origin,
            news,
            welcome,
        }},
    }
</script>

<style scoped>

</style>