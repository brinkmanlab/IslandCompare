<template>
    <div id="app">
        <Navigation></Navigation>
        <router-view />
        <footer>
            <span>Please read our <a href="#">license</a> for any publications resulting in the use of our services.</span>
        </footer>
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
    :root {
        --color-primary-0: #9282DF;	/* Main Primary color */
        --color-primary-1: #FFFFFF;
        --color-primary-2: #C7BEF6;
        --color-primary-3: #6452BC;
        --color-primary-4: #3C2B92;

        --color-secondary-1-0: #FFC686;	/* Main Secondary color (1) */
        --color-secondary-1-1: #FFFFFF;
        --color-secondary-1-2: #FFE0BD;
        --color-secondary-1-3: #FFB25B;
        --color-secondary-1-4: #D48329;

        --color-secondary-2-0: #76BBD9;	/* Main Secondary color (2) */
        --color-secondary-2-1: #FFFFFF;
        --color-secondary-2-2: #B7E1F4;
        --color-secondary-2-3: #438EAF;
        --color-secondary-2-4: #1F6685;

        --color-complement-0: #FFE786;	/* Main Complement color */
        --color-complement-1: #FFFFFF;
        --color-complement-2: #FFF2BD;
        --color-complement-3: #FFDF5B;
        --color-complement-4: #D4B229;
    }

    #app {
        font-family: 'Avenir', Helvetica, Arial, sans-serif;
        -webkit-font-smoothing: antialiased;
        -moz-osx-font-smoothing: grayscale;
        margin: 10%;
    }

    #app .Toast {
        position: fixed;
        top: 20%;
        left: 50%;
        transform: translate(-50%, -50%);
    }

    #app .Navigation {
        position: fixed;
        background-color: var(--color-secondary-2-2);
        font-size: 1.2em;
        top: 0;
        left: 0;
        width: 100%;
    }

    #app footer {
        grid-area: footer;
        text-align: center;
    }
</style>
