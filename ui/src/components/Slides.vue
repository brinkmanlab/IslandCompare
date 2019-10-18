<template>
    <b-carousel class="Slides" controls indicators :background="background" :interval="interval" :no-wrap="!loop">
        <b-carousel-slide v-for="(item, index) of slides" :key="index" :class="item.slug">
            <template v-slot:img>
                <img v-if="item.img" class="slide-img" :src="item.img" :alt="item.alt"/>
                <h1 v-if="item.title">{{ item.title }}</h1>
                <div style="display: contents" v-html="item.body"></div>
            </template>
        </b-carousel-slide>
    </b-carousel>
</template>

<script>
    export default {
        name: "Slides",
        props: {
            slides: {
                type: Array,
                required: true,
            },
            interval: {
                type: Number,
                default: 3000
            },
            loop: {
                type: Boolean,
                default: true,
            },
            background: {
                type: String,
                default: '#FFFFFF',
            },
        }
    }
</script>

<style scoped>
    .carousel {
        padding: 0.2em;
        overflow: hidden;
    }

    img {
        float: right;
        margin: 1em;
        width: 50%;
        background: white;
    }

    .carousel-item:nth-child(even) img {
        float: left;
    }

    .carousel >>> .carousel-control-prev-icon {
        background-image: url("data:image/svg+xml,%3csvg xmlns='http://www.w3.org/2000/svg' fill='%237bbdf7' viewBox='0 0 8 8'%3e%3cpath d='M5.25 0l-4 4 4 4 1.5-1.5-2.5-2.5 2.5-2.5-1.5-1.5z'/%3e%3c/svg%3e");
    }

    .carousel >>> .carousel-control-next-icon {
        background-image: url("data:image/svg+xml,%3csvg xmlns='http://www.w3.org/2000/svg' fill='%237bbdf7' viewBox='0 0 8 8'%3e%3cpath d='M2.75 0l-1.5 1.5 2.5 2.5-2.5 2.5 1.5 1.5 4-4-4-4z'/%3e%3c/svg%3e");
    }

    .carousel >>> .carousel-control-prev, .carousel >>> .carousel-control-next {
        opacity: 0.2;
    }

    .carousel:hover >>> .carousel-control-prev, .carousel:hover >>> .carousel-control-next {
        opacity: 0.8;
    }

    .carousel >>> .carousel-indicators li {
        background-color: #7bbdf7;
    }
</style>