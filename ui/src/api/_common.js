import axios from 'axios';

export function query(options) {
    return axios(options);
}

export let BaseModule = {
    namespaced: true,
};

export function BaseState(state) {
    if (state instanceof Function) state = state();
    return ()=>{return {
        ...state,
    }};
}

export let BaseMutations = {  };
export let BaseActions = {  };
export let BaseGetters = {  };