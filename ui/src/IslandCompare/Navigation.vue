<template>
    <b-navbar :sticky="true" variant="primary" type="dark" toggleable="sm">
        <b-navbar-toggle target="nav-collapse"/>
        <slot name="before" />
        <b-navbar-nav>
            <b-collapse is-nav id="nav-collapse">
                <!--b-nav-item-dropdown v-bind:text="$root.appName">
                    <b-dropdown-item to="/">{{ home.name }}</b-dropdown-item>
                    <b-dropdown-item v-for="service in services" v-bind:key="service.name" v-bind:href="service.path">{{ service.name }}</b-dropdown-item>
                </b-nav-item-dropdown-->
                <b-nav-item v-for="page in pages" v-bind:key="page.name" v-bind:to="page.path" v-bind:active="$route.path === page.path">
                    {{ page.name }}
                </b-nav-item>
                <slot />
            </b-collapse>
        </b-navbar-nav>
        <slot name="after" />
    </b-navbar>
</template>

<script>
    import services from "@/brinkman_services";
    export default {
        name: "Navigation",
        data() {return{
            home: this.$router.options.routes.find(o=>o.path === "/"),
            services: services,
        }},
        computed: {
            pages() {
                return this.$router.options.routes.filter(route=>
                    /*route.path !== "/" // Don't include home button
                    && */(!route.hasOwnProperty('meta') || !route.meta.hasOwnProperty('navbar') || route.meta.navbar === true) // Allow to be hidden using meta tag
                );
            }
        },
    }
</script>

<style scoped>

</style>