<template>
    <b-navbar :sticky="true" variant="primary" type="dark" toggleable="sm">
        <b-navbar-toggle target="nav-collapse"/>
        <b-navbar-nav>
            <b-collapse is-nav id="nav-collapse">
                <b-nav-item-dropdown v-bind:text="$root.appName">
                    <!-- TODO v-bind:class="[$route.path === home.path ? 'pure-menu-selected' : '' ]" -->
                    <b-dropdown-item to="/">{{ home.name }}</b-dropdown-item>
                    <b-dropdown-item v-for="service in services" v-bind:key="service.name" v-bind:href="service.path">{{ service.name }}</b-dropdown-item>
                </b-nav-item-dropdown>
                <b-nav-item v-for="page in pages" v-bind:key="page.name" v-bind:to="page.path">
                    <!-- TODO v-bind:class="[$route.path === page.path ? 'pure-menu-selected' : '' ]" -->
                    {{ page.name }}
                </b-nav-item>
            </b-collapse>
        </b-navbar-nav>
    </b-navbar>
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