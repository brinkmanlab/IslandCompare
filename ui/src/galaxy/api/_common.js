import { Model as VuexModel } from '@vuex-orm/core';

const HasState = {
    /**
     * This mixin expects the class to have a static list of end_states and to extend Model
     */

    /**
     * Poll model state until it is an end_state
     * @param callback Called when state changes to end_state, return true to stop polling, false to continue
     * @param extra Additional parameters passed to Model.start_polling()
     */
    poll_state(callback = ()=>true, ...extra) {
        let self = this;
        if (!this.constructor.end_states.includes(self.state)) {
            this.start_polling(()=>{
                self = self.constructor.find(self.id); // TODO recover from the reactivity system failing
                if (self.constructor.end_states.includes(self.state)) {
                    return callback();
                }
                return false;
            }, ...extra);
        }
    }
};

class Model extends VuexModel {

    static fields() {
        return {
        }
    }

    async post(fields, params = {}) {
        let data = this;
        if (fields !== undefined) {
            //Only keep requested fields
            data = Object.fromEntries(
                Object.entries(data).filter(([k,_])=>fields.includes(k)) //eslint-disable-line
            );
        }
        if (data.hasOwnProperty('$toJson')) data = data.$toJson();
        else data = JSON.stringify(data);
        return await this.constructor.$update({
            params: {
                id: this[this.constructor.primaryKey],
                url: this.get_base_url(),
                ...params,
            },
            data: data,
        });
    }

    // TODO find all places this can be used and replace code
    async reload(options={}) {
        return await this.constructor.$get({
            ...options,
            params: {
                url: this.get_base_url(),
                id: this.id,
                ...options.params,
            },
        });
    }

    async delete(options = {}) {
        this.stop_polling();
        if (this.hid === -1) {
            //Delete locally if ghost item
            this.constructor.delete(this.id);
            return;
        }
        this.constructor.delete(this.id); // TODO this shouldn't be needed, the next command calls it but fails
        await this.constructor.$delete({
            ...options,
            params: {
                id: this[this.constructor.primaryKey],
                ...options.params,
            },
        });
    }

    get_base_url() {
        return '';
    }

    start_polling(stop_criteria=null, options={}, interval=10000) {
        const self = this;
        if (!window.hasOwnProperty('pollHandles')) window.pollHandles = new Map();
        if (!window.pollHandles.has(this.id)) {
            const f = ()=>{ //TODO this can be done better :/
                self.reload(options).then(() => {
                    if (typeof stop_criteria === "function" && stop_criteria()) {
                        this.stop_polling();
                    } else {
                        // Reschedule after reload
                        window.pollHandles.set(this.id, setTimeout(f, interval));
                    }
                });
            };
            f();
        }
    }

    stop_polling() {
        if (window.hasOwnProperty('pollHandles')) {
            const pollHandle = window.pollHandles.get(this.id);
            if (pollHandle !== undefined) {
                clearTimeout(pollHandle);
                window.pollHandles.delete(this.id);
            }
        }
    }

    static beforeDelete(model) {
        model.stop_polling();
        super.beforeDelete(model);
    }

    // TODO find all places this can be used and replace code
    static async findOrLoad(id, url='', options={}) {
        const result = this.find(id);
        if (result) return result;
        const option_params = 'params' in options ? options.params : {};
        await this.$get({
            ...options,
            params: {
                id: id,
                url: url,
                ...option_params,
            },
        });
        return this.find(id);
    }

    is(obj) {
        if (obj.constructor.hasOwnProperty('primaryKey')) {
            //TODO are ids unique across Galaxy?
            return this[this.constructor.primaryKey] === obj[obj.constructor.primaryKey];
        } else {
            return super.is(obj);
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
    HasState,
}