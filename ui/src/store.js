import Vue from 'vue'
import Vuex from 'vuex'

Vue.use(Vuex);

export default new Vuex.Store({
    state: {
        galaxy: {
            api_key: '',
            current_history: '',
        }
    },
    mutations: {
    },
    actions: {
        getHistoryItems(history_id) { if (history_id) return []; return [] },
    }
})
