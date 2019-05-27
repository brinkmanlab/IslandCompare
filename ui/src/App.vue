<template>
    <div id="app">
        <Navigation></Navigation>
        <router-view class="content"/>
        <b-navbar fixed="bottom" variant="primary" type="dark" size="sm" class="tosfooter">
            <b-nav-text class="mx-auto tosfooter_content">
                By using this service you agree to the <b-link to="/terms">terms of service</b-link>.
                For any publications resulting in the use of our services see <b-link to="/publications">how to cite this service</b-link>.
            </b-nav-text>
        </b-navbar>
    </div>
</template>

<script>
    import Navigation from '@/IslandCompare/Navigation'
    export default {
        name: 'app',
        components: {
            Navigation,
        },
        data() { return {
            message: '',
        }},
        methods: {
        },
        mounted() {
        },
        errorCaptured(err, vm, info) { //eslint-disable-line
            this.$bvModal.msgBoxOk(err, {
                title: 'Error',
                size: 'sm',
                buttonSize: 'sm',
                okVariant: 'danger',
                headerClass: 'p-2 border-bottom-0',
                footerClass: 'p-2 border-top-0',
                centered: true
            });
            return true;
        }
    }
</script>

<style>

</style>
