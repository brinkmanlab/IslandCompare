<template>
    <div id="app">
        <Navigation />
        <keep-alive>
            <router-view class="content"/>
        </keep-alive>
        <footer class="container">
            <b-row>
                <b-col>
                    <p>If you use IslandCompare please cite (a publication is underway, in the meantime please cite IslandViewer):</p>
                    <p>Bertelli, C., Gray, K.L. et al. “Enabling genomic island prediction and comparison in multiple genomes to investigate bacterial evolution and outbreaks” Microbial Genomics. 2022.
                        doi: <b-link href="https://doi.org/10.1099/mgen.0.000818">10.1099/mgen.0.000818</b-link></p>
                </b-col>
                <b-col cols="2">
                    <p>Managed by the <b-link href="https://www.brinkman.mbb.sfu.ca/">Brinkman Lab</b-link> and the Bertelli Lab</p>
                </b-col>
            </b-row>
            <b-row>
                <b-col>
                    <div class="logos">
                        <img src="/logos/cihr.png" />
                        <img src="/logos/genome_atlantic.jpg" />
                        <img src="/logos/genome-bc.png" />
                        <img src="/logos/genome_canada.jpg" />
                        <img src="/logos/computecanada.png" />
                        <img src="/logos/sfu.png" />
                    </div>
                </b-col>
            </b-row>
        </footer>
        <!--b-navbar fixed="bottom" variant="primary" type="dark" size="sm" class="tosfooter">
            <b-nav-text class="mx-auto tosfooter_content">
                For any publications resulting in the use of our services see <b-link to="/publications">how to cite this service</b-link>.
                By using this service you agree to the <b-link to="/terms">terms of service</b-link>.
            </b-nav-text>
        </b-navbar-->
        <v-tour name="tour" :steps="tour.steps($tours['tour'])" :callbacks="tour.callbacks" />
    </div>
</template>

<script>
    import Navigation from '@/IslandCompare/Navigation';
    import Tour from './tour'
    export default {
        name: 'app',
        components: {
            Navigation,
        },
        data() { return {
            message: '',
            tour: Tour,
        }},
        methods: {
            start_tour() {
                this.$router.push({path: '/analysis', query: {tour: 'tour'}});
            },
            show_error(err) {
                if (err && err.hasOwnProperty('name') && err.name === 'NavigationDuplicated') return;
                if (!err) return;
                this.$bvModal.msgBoxOk(err.message || err, {
                    title: 'Error',
                    size: 'sm',
                    buttonSize: 'sm',
                    okVariant: 'danger',
                    headerClass: 'p-2 border-bottom-0',
                    footerClass: 'p-2 border-top-0',
                    centered: true
                });
            }
        },
        computed: {
            showTour() {
                return (!this.$route.hasOwnProperty('meta') || !this.$route.meta.hasOwnProperty('showTour') || this.$route.meta.showTour === true); // Allow to be hidden using meta tag
            }
        },
        created() {
            window.onerror = this.show_error.bind(this);
            window.onunhandledrejection = (error) => { window.onerror(error.reason || error)};
        },
        errorCaptured(err, vm, info) { //eslint-disable-line
            this.show_error(err);
            return true;
        }
    }
</script>

<style>
    .v-step {
        /* bootstrap sets some elements to z=2 */
        z-index: 10 !important;
    }

    .logos {
        display: flex;
        flex-direction: row;
        justify-content: space-evenly;
        background: white;
        padding: 2px;
        border-radius: 2px;
    }

    .logos img {
        max-height: 7vh;
        width: auto;
    }
</style>
