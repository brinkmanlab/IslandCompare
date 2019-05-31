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
        this.stop_polling();
        if (this.hid <= 0) {
            //Delete locally if ghost item
            this.constructor.delete(this.id);
            return;
        }
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
        if (!window.hasOwnProperty('pollHandles')) window.pollHandles = new Map();
        if (!window.pollHandles.has(this.id)) {
            let pollHandle = setInterval(() => {
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
            window.pollHandles.set(this.id, pollHandle);
        }
    }

    stop_polling() {
        if (window.hasOwnProperty('pollHandles')) {
            let pollHandle = window.pollHandles.get(this.id);
            if (pollHandle !== undefined) {
                clearInterval(pollHandle);
                window.pollHandles.delete(this.id);
            }
        }
    }

    static beforeDelete(model) {
        model.stop_polling();
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