import Vue from 'vue'
import VueRouter from 'vue-router'
import AsyncComputed from 'vue-async-computed'
import BootstrapVue from 'bootstrap-vue'
//import VueAnalytics from 'vue-analytics'
import VueTour from 'vue-tour';

import App from './App.vue'
import { store } from './store'

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

//TODO import i18n from './i18n'

//const isProd = process.env.NODE_ENV === 'production';

Vue.config.productionTip = false;

Vue.use(VueRouter);
Vue.use(AsyncComputed);
Vue.use(BootstrapVue);
Vue.use(VueTour);

const router = new VueRouter({
    mode: 'history',
    base: base_path,
    routes: [
        { path: '/', component: Home, name: "IslandCompare" },
        { path: '/about', component: HTMLFragment, name: "About", props: {content: require('html-loader!@/assets/about.htm')} },
        { path: '/faq', component: HTMLFragment, name: "FAQ", props: {content: require('html-loader!@/assets/faq.htm')} },
        { path: '/publications', component: HTMLFragment, name: "Publications", props: {content: require('html-loader!@/assets/publications.htm')} },
        //{ path: '/download', component: HTMLFragment, name: "Download", props: {content: require('html-loader!@/assets/download.htm')} },
        { path: '/contact', component: HTMLFragment, name: "Contact", props: {content: require('html-loader!@/assets/contact.htm')} },
        { path: '/analysis', component: AsyncAnalysis, name: "Run Analysis", props: route=>({tour: route.query.tour || ''}) },
        { path: '/history', component: AsyncHistory, name: "Job History", /*meta: {navbar: getUUID() == true}*/},
        { path: '/terms', component: HTMLFragment, name: "Terms of Use", props: {content: require('html-loader!@/assets/terms.htm')}, meta: {navbar: false} },
        {
            path: '/visualize',
            component: IFrameContent,
            props: route => ({
                src: `${galaxy_path}/plugins/visualizations/islandcompare/static/islandcompare.html?src=${route.query.src}`,
                name: 'visualize'
            }),
        },
        {
            path: '/visualize/:id',
            component: IFrameContent,
            props: route=>({
                src: `${galaxy_path}/plugins/visualizations/islandcompare/static/islandcompare.html?src=${galaxy_path}/datasets/${route.params.id}/display`,
                name: 'visualize'
            }),
        },
    ]
});

/*Vue.use(VueAnalytics, {
    id: 'UA-46024702-13',
    router,
    autoTracking: {
        exception: true,
    },
    debug: {
        enabled: isProd,
        sendHitTask: isProd,
    },
});*/

new Vue({
    store,
    router,
    //  i18n,
    render: h => h(App),
    data() {return{
        appName: "IslandCompare",
    }},
}).$mount('#app');