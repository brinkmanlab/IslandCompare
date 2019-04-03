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

    async upload(params = {}) {
        return await this.constructor.$update({
            params: {
                id: this[this.constructor.primaryKey],
                ...params,
            },
            data: this.$toJson(),
        });
    }

    async delete(params = {}) {
        return await this.constructor.$delete({
            params: {
                id: this[this.constructor.primaryKey],
                ...params,
            },
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