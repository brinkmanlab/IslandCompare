import { createStore } from 'vuex';

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

const store = createStore({
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

