import { BaseModule, BaseState, BaseMutations, BaseActions, BaseGetters } from './_common' //eslint-disable-line no-unused-vars

let VuexModule = {
    ...BaseModule,
    state: BaseState({}),
    mutations: {
        ...BaseMutations,
    },
    actions: {
        ...BaseActions,
    },
    getters: {
        ...BaseGetters,
    },
};

export default VuexModule;