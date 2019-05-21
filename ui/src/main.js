import Vue from 'vue'
import VueRouter from 'vue-router'
import AsyncComputed from 'vue-async-computed'

import App from './App.vue'
import { getStore } from './store'

import Analysis from "@/IslandCompare/Analysis";
import Markdown from "./components/Markdown"

//TODO import i18n from './i18n'

Vue.config.productionTip = false;

Vue.use(VueRouter);
Vue.use(AsyncComputed);

const router = new VueRouter({
    routes: [
        { path: '/', component: Markdown, name: "IslandCompare", props: {src: 'home.md'} },
        { path: '/about', component: Markdown, name: "About", props: {src: 'about.md'} },
        { path: '/analysis', component: Analysis, name: "Run Analysis" },
        { path: '/faq', component: Markdown, name: "FAQ", props: {src: 'faq.md'} },
        { path: '/download', component: Markdown, name: "Download", props: {src: 'download.md'} },
        { path: '/publications', component: Markdown, name: "Publications", props: {src: 'publications.md'} },
        { path: '/contact', component: Markdown, name: "Contact", props: {src: 'contact.md'} },
        { path: '/terms', component: Markdown, name: "Terms of Use", props: {src: 'terms.md'}, meta: {navbar: false} },
    ]
});

getStore().then(store=>{
    new Vue({
        store,
        router,
        //  i18n,
        render: h => h(App),
    }).$mount('#app');
});