import {BaseModule, BaseState, BaseMutations, BaseActions, BaseGetters, query} from './_common' //eslint-disable-line no-unused-vars

let VuexModule = {
    ...BaseModule,
    state: BaseState({}),
    mutations: {
        ...BaseMutations,
        delete(state, payload) {
            var history_content = state[payload.history_id];
            var index = history_content.findIndex(i=>i.id === payload.id && i.history_content_type === payload.history_content_type);
            if (index >= 0) history_content.splice(index, 1);
            else throw "State desync detected while deleting " + payload.id;
        }
    },
    actions: {
        ...BaseActions,
        async delete(context, history_id, id, history_content_type) {
            await query({
                method: 'delete',
                url: `/api/histories/${history_id}/contents/${history_content_type}/${id}`,
            });
            context.commit({
                type: 'delete',
                history_id: history_id,
                id: id,
                history_content_type: history_content_type,
            });
        },
        async get(context, history_id, filter='') {
            if (!context.state[history_id]) {
                await query({
                    method: 'get',
                    url: `/api/histories/${this.history_id}/contents/`,
                    params: {
                        type: filter,
                    },
                });
            }
        },
    },
    getters: {
        ...BaseGetters,
    },
};

export default VuexModule;