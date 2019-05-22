<template>
    <div id="app">
        <Navigation></Navigation>
        <router-view class="content"/>
        <b-navbar fixed="bottom" variant="dark" type="dark" size="sm" class="tosfooter">
            <b-nav-text class="mx-auto tosfooter_content">
                By using this service you agree to the <b-link to="/terms">terms of service</b-link>.
                For any publications resulting in the use of our services see <b-link to="/publications">how to cite this service</b-link>.
            </b-nav-text>
        </b-navbar>
        <Toast ref="toast" />
    </div>
</template>

<script>
    import Toast from '@/components/Toast'
    import Navigation from '@/IslandCompare/Navigation'
    export default {
        name: 'app',
        components: {
            Toast,
            Navigation,
        },
        data() { return {
        }},
        methods: {
            toast(e) {
                this.$refs['toast'].show(e);
                console.log(e); // eslint-disable-line no-console
            },
        },
        mounted() {
            //TODO wait for galaxy server
        },
        //errorCaptured(err, vm, info) { //eslint-disable-line
        //    this.toast(err);
        //    return true;
        //}
    }
</script>

<style>

</style>
