import Vue from 'vue'
import Vuex, { Store } from 'vuex'

//TODO vuex-persistedstate for current session


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
    ]
});

// Async load galaxy ORM as it is BIG
const galaxy_load = import('@/galaxy').then(module=>{
    //Lazy load galaxy ORM as it is BIG
    module.register(store);
    return module;
});

export {
    store,
    galaxy_load,
}

