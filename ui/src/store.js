import Vue from 'vue'
import Vuex, { Store } from 'vuex'
import VuexPersistence from 'vuex-persist'
import { galaxy_path } from "@/app.config";
import Galaxy, {register} from "@/galaxy/src";

Vue.use(Vuex);

const persistence = new VuexPersistence({
    reducer(state) {
        return {
            galaxy: {
                $name: state.galaxy.$name,
                Genome: state.galaxy.Genome,
                StoredWorkflow: state.galaxy.StoredWorkflow,
            }
        };
    }
});

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
        //persistence.plugin,
    ]
});

register(store, {baseURL: galaxy_path});

// Async load galaxy ORM as it is BIG
const galaxy_load = Promise.resolve(Galaxy);/*import('@/galaxy/src/').then(module=>{
    //Lazy load galaxy ORM as it is BIG
    module.register(store, {baseURL: galaxy_path});
    return module;
});*/

export {
    store,
    galaxy_load,
}

