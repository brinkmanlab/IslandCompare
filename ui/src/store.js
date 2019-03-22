import Vue from 'vue'
import Vuex, { Store } from 'vuex'
import VuexORM from '@vuex-orm/core'
import database from './galaxy'
import VueCookies from 'vue-cookies'
//TODO vuex-persistedstate for current session
Vue.use(VueCookies);

// set default config
VueCookies.config('7d');

Vue.use(Vuex);

const store = new Store({
    namespaced: true,
    state: {
    },
    mutations: {
    },
    actions: {
    },
    modules: {
    },
    plugins: [
        VuexORM.install(database, { namespace: 'galaxy' }),
    ]
});

export {
    store as default,
}

