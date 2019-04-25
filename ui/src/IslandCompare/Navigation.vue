<template>
    <div class="Navigation pure-menu pure-menu-horizontal">
        <ul class="pure-menu-list">
            <li class="pure-menu-item pure-menu-has-children pure-menu-allow-hover" v-bind:class="[$route.path === home.path ? 'pure-menu-selected' : '' ]">
                <router-link to="/" class="pure-menu-link">{{ home.name }}</router-link>
                <ul class="pure-menu-children">
                    <li v-for="service in services"
                        v-bind:key="service.name"
                        class="pure-menu-item">
                        <a v-bind:href="service.path" class="pure-menu-link">{{ service.name }}</a>
                    </li>
                </ul>
            </li>
            <li v-for="page in pages"
                v-bind:key="page.name"
                class="pure-menu-item"
                v-bind:class="[$route.path === page.path ? 'pure-menu-selected' : '' ]"
            >
                <router-link v-bind:to="page.path" class="pure-menu-link" v-bind:key="page.name">{{ page.name }}</router-link>
            </li>
        </ul>
    </div>
</template>

<script>
    export default {
        name: "Navigation",
        data() {return{
            home: this.$router.options.routes.find(o=>o.path === "/"),
            services: [
                {name: "IslandViewer", path: "#"},
                {name: "PSORTdb", path: "#"},
                {name: "Pseudomonas", path: "#"},
            ],
            pages: this.$router.options.routes.filter(o=>o.path !== "/"),
        }},
    }
</script>

<style scoped>

</style>