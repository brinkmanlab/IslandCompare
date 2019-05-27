import Vue from 'vue'
import VueRouter from 'vue-router'
import AsyncComputed from 'vue-async-computed'
import BootstrapVue from 'bootstrap-vue'

import App from './App.vue'
import { store } from './store'

import './assets/main.scss'

const AsyncAnalysis = () => import("@/IslandCompare/Analysis");
const AsyncHistory = () => import("@/IslandCompare/JobHistory");
//import Markdown from "./components/Markdown"
import HTMLFragment from './components/HTMLFragment'

//TODO import i18n from './i18n'

Vue.config.productionTip = false;

Vue.use(VueRouter);
Vue.use(AsyncComputed);
Vue.use(BootstrapVue);

const router = new VueRouter({
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
    ]
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