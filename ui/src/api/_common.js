import { Model as VuexModel } from '@vuex-orm/core';

class Model extends VuexModel {

    static fields() {
        return {

        }
    }

    /*async download() {
        let params = {
            id: this[this.prototype.primaryKey],
        };
        if (parent.hasOwnProperty('url')) params.url = parent.url;
        await this.prototype.$get({
            params: params,
        });
        return this;
    }*/

    async upload() {
        return await this.prototype.$update({
            params: {
                id: this[this.prototype.primaryKey],
            },
            data: this.$toJson(),
        });
    }
}

const Module = {
    namespaced: true,
};

function State(state) {
    if (state instanceof Function) state = state();
    return ()=>{return {
        ...state,
    }};
}

const Mutations = {  };
const Actions = {  };
const Getters = {  };

export {
    Model,
    Module,
    State,
    Mutations,
    Actions,
    Getters,
}