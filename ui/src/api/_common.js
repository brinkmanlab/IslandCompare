import { Model as VuexModel } from '@vuex-orm/core';

class Model extends VuexModel {

    static fields() {
        return {
            _pollHandle: this.attr(null),
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

    get_base_url() {
        return '';
    }

    start_polling(stop_criteria=null, interval=10000) {
        if (this._pollHandle === null) {
            this._pollHandle = setInterval(() => {
                this.constructor.$get({
                    params: {
                        url: this.get_base_url(),
                        id: this.id,
                    },
                }).then(() => {
                    if (typeof stop_criteria === "function" && stop_criteria()) {
                        this.stop_polling();
                    }
                });
            }, interval);
        }
    }

    stop_polling() {
        if (this._pollHandle) {
            clearInterval(this._pollHandle);
            this._pollHandle = null;
        }
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