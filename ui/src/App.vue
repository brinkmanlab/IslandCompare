<template>
    <div id="app">
        <Navigation></Navigation>
        <router-view class="content"/>
        <footer>
            <p>If you use IslandCompare please cite (a publication is underway, in the meantime please cite IslandViewer):</p>
            <p>Bertelli, C. et al. "IslandViewer 4: Expanded prediction of genomic islands for larger-scale datasets" Nucleic Acids Research. 2017 May 2.
                doi: <b-link href="https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkx343">10.1093/nar/gkx343</b-link></p>
        </footer>
        <b-navbar fixed="bottom" variant="primary" type="dark" size="sm" class="tosfooter">
            <b-nav-text class="mx-auto tosfooter_content">
                For any publications resulting in the use of our services see <b-link to="/publications">how to cite this service</b-link>.
                By using this service you agree to the <b-link to="/terms">terms of service</b-link>.
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
            window.onerror = this.errorCaptured
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
