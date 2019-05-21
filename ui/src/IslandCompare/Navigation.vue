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
                {name: "IslandViewer", path: "http://www.pathogenomics.sfu.ca/islandviewer/"},
		{name: "PSORTb v.3.0", path: "https://www.psort.org/psortb/index.html"},
                {name: "PSORTdb", path: "https://db.psort.org"},
                {name: "Pseudomonas Genome Database", path: "http://www.pseudomonas.com"},
		{name: "Burkholderia Genome Database", path: "http://www.burkholderia.com"},
            ],
            pages: this.$router.options.routes.filter(route=>route.path !== "/" && (!route.hasOwnProperty('meta') || !route.meta.hasOwnProperty('navbar') || route.meta.navbar === true)),
        }},
    }
</script>

<style scoped>

</style>