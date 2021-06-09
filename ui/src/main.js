import Vue from 'vue'
import VueRouter from 'vue-router'
import AsyncComputed from 'vue-async-computed'
import BootstrapVue from 'bootstrap-vue'
import VueAnalytics from 'vue-analytics'
import VueTour from 'vue-tour';

import Galaxy from 'galaxy-client'

import App from './App.vue'
import store from './store'

import './assets/main.scss';
import 'vue-tour/dist/vue-tour.css';

// Create a filter for links to galaxy
import {base_path, galaxy_path} from "./app.config";
Vue.filter('galaxybase', value=>galaxy_path+value);

const AsyncAnalysis = () => import("@/IslandCompare/Analysis");
const AsyncHistory = () => import("@/IslandCompare/JobHistory");
import Home from '@/IslandCompare/Home';
import HTMLFragment from '@/components/HTMLFragment'
import IFrameContent from "@/components/IFrameContent";
import News from "@/components/News";

//TODO import i18n from './i18n'

const isProd = process.env.NODE_ENV === 'production';

Vue.config.productionTip = false;

Vue.use(Galaxy, {store, baseURL: galaxy_path});
Vue.use(VueRouter);
Vue.use(AsyncComputed, {useRawError: true,});
Vue.use(BootstrapVue);
Vue.use(VueTour);

const router = new VueRouter({
    mode: 'history',
    base: base_path,
    routes: [ // The order of these determines the order of they appear in the menu bar
        { path: '/', component: Home, name: "IslandCompare" },
        { path: '/analysis', component: AsyncAnalysis, name: "Analyse", props: route=>({tour: route.query.tour || ''}) },
        { path: '/history', component: AsyncHistory, name: "Job History", meta: {showTour: false}},
        { path: '/about', component: HTMLFragment, name: "About", props: {content: require('html-loader!@/assets/about.htm'), slug:'about'} },
        { path: '/faq', component: HTMLFragment, name: "FAQ", props: {content: require('html-loader!@/assets/faq.htm'), slug: 'faq'} },
        { path: '/download', component: HTMLFragment, name: "Download", props: {content: require('html-loader!@/assets/download.htm')} },
        { path: '/publications', component: HTMLFragment, name: "Publications", props: {content: require('html-loader!@/assets/publications.htm'), slug: 'publications'} },
        { path: '/contact', component: HTMLFragment, name: "Contact", props: {content: require('html-loader!@/assets/contact.htm'), slug: 'contact'} },
        //{ path: '/terms', component: HTMLFragment, name: "Terms of Use", props: {content: require('html-loader!@/assets/terms.htm')}, meta: {navbar: false} },
        { path: '/news', component: News, name: "News", props: {news: require('@/assets/news')}, meta: {navbar: false} },
        {
            path: '/visualize',
            component: IFrameContent,
            meta: {navbar: false, showTour: false},
            props: route => ({
                src: `${galaxy_path}/plugins/visualizations/islandcompare/static/index.html?src=${route.query.src}`,
                name: 'visualize'
            }),
        },
        {
            path: '/visualize/:id',
            component: IFrameContent,
            meta: {navbar: false, showTour: false},
            props: route=>({
                src: `${galaxy_path}/plugins/visualizations/islandcompare/static/index.html?src=${galaxy_path}/api/histories/any/contents/${route.params.id}/display?key=${route.query.key}`,
                name: 'visualize'
            }),
        },
    ]
});

const analytics_id = {
    'islandcompare.pathogenomics.ca': 'UA-46024702-13',
    'islandcompare.pathogenomics.sfu.ca': 'UA-46024702-14',
}[window.location.hostname];

if (analytics_id) {
    Vue.use(VueAnalytics, {
        id: analytics_id,
        router,
        autoTracking: {
            exception: true,
        },
        debug: {
            enabled: !isProd,
            sendHitTask: isProd,
        },
    });
} else console.log('Unknown analytics hostname'); //eslint-disable-line

new Vue({
    store,
    router,
    //  i18n,
    render: h => h(App),
    data() {return{
        appName: "IslandCompare",
    }},
}).$mount('#app');
