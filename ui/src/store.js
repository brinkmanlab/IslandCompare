import Vue from 'vue'
import Vuex, { Store } from 'vuex'
//import VuexPersistence from 'vuex-persist'

Vue.use(Vuex);

/*const persistence = new VuexPersistence({
    reducer(state) {
        return {
            galaxy: {
                $name: state.galaxy.$name,
                Genome: state.galaxy.Genome,
                StoredWorkflow: state.galaxy.StoredWorkflow,
            }
        };
    }
});*/

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

export {
    store as default,
}

