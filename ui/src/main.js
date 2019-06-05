import Vue from 'vue'
import VueRouter from 'vue-router'
import AsyncComputed from 'vue-async-computed'
import BootstrapVue from 'bootstrap-vue'
import VueAnalytics from 'vue-analytics'

import App from './App.vue'
import { store } from './store'

import './assets/main.scss'

// Create a filter for links to galaxy
import {base_path, galaxy_path} from "@/app.config";
Vue.filter('galaxybase', value=>galaxy_path+value);

const AsyncAnalysis = () => import("@/IslandCompare/Analysis");
const AsyncHistory = () => import("@/IslandCompare/JobHistory");
//import Markdown from "./components/Markdown"
import HTMLFragment from './components/HTMLFragment'
import IFrameContent from "@/components/IFrameContent";

//TODO import i18n from './i18n'

const isProd = process.env.NODE_ENV === 'production';

Vue.config.productionTip = false;

Vue.use(VueRouter);
Vue.use(AsyncComputed);
Vue.use(BootstrapVue);

const router = new VueRouter({
    mode: 'history',
    base: base_path,
    routes: [
        { path: '/', component: HTMLFragment, name: "IslandCompare", props: {src: 'home.htm'} },
        { path: '/about', component: HTMLFragment, name: "About", props: {src: 'about.htm'} },
        { path: '/analysis', component: AsyncAnalysis, name: "Run Analysis" },
        { path: '/history', component: AsyncHistory, name: "Job History", /*meta: {navbar: getUUID() == true}*/},
        { path: '/faq', component: HTMLFragment, name: "FAQ", props: {src: 'faq.htm'} },
        //{ path: '/download', component: HTMLFragment, name: "Download", props: {src: 'download.htm'} },
        { path: '/publications', component: HTMLFragment, name: "Publications", props: {src: 'publications.htm'} },
        { path: '/contact', component: HTMLFragment, name: "Contact", props: {src: 'contact.htm'} },
        { path: '/terms', component: HTMLFragment, name: "Terms of Use", props: {src: 'terms.htm'}, meta: {navbar: false} },
        { path: '/visualize/:id', component: IFrameContent, props: route=>({src: `${galaxy_path}/plugins/visualizations/islandcompare/show?dataset_id=${route.params.id}`, name: 'visualize'})},
    ]
});

Vue.use(VueAnalytics, {
    id: 'UA-46024702-13',
    router,
    autoTracking: {
        exception: true,
    },
    debug: {
        enabled: isProd,
        sendHitTask: isProd,
    },
});

new Vue({
    store,
    router,
    //  i18n,
    render: h => h(App),
    data() {return{
        appName: "IslandCompare",
    }},
}).$mount('#app');