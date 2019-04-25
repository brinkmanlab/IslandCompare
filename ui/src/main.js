import Vue from 'vue'
import VueRouter from 'vue-router'
import AsyncComputed from 'vue-async-computed'
import App from './App.vue'
import { getStore } from './store'

import About from "@/IslandCompare/About";
import Analysis from "@/IslandCompare/Analysis";
import FAQ from "@/IslandCompare/FAQ";
import Download from "@/IslandCompare/Download";
import Publications from "@/IslandCompare/Publications";
import Contact from "@/IslandCompare/Contact";
import Home from "@/IslandCompare/Home";

//TODO import i18n from './i18n'

Vue.config.productionTip = false;

Vue.use(VueRouter)
Vue.use(AsyncComputed);

//TODO Replace static .vue files with single vue that imports flat html for about, faq, download, pub, contact
const router = new VueRouter({
    routes: [
        { path: '/', component: Home, name: "IslandCompare" },
        { path: '/about', component: About, name: "About" },
        { path: '/analysis', component: Analysis, name: "Run Analysis" },
        { path: '/faq', component: FAQ, name: "FAQ" },
        { path: '/download', component: Download, name: "Download" },
        { path: '/publications', component: Publications, name: "Publications" },
        { path: '/contact', component: Contact, name: "Contact" },
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