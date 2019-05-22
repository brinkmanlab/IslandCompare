import Vue from 'vue'
import Vuex, { Store } from 'vuex'
import VuexORM from '@vuex-orm/core'
import { getDatabase } from './galaxy'

//TODO vuex-persistedstate for current session


Vue.use(Vuex);

async function getStore() {
    let database = await getDatabase();
    return new Store({
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
}

export {
    getStore,
}

