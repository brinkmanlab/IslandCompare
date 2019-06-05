import * as Common from "./_common";
import {History} from "@/api/histories";



class Job extends Common.Model {
    static entity = 'Jobs';
    static primaryKey = 'id';
    static end_states = ['ok', 'error', 'deleted', 'deleted_new'];

    static fields() {
        return {
            ...super.fields(),
            id: this.string(null),
            tool_id: this.string().nullable(),
            update_time: this.string().nullable(),
            history_id: this.string().nullable(),
            exit_code: this.number().nullable(),
            state:	this.string().nullable(), // ‘new’, ‘upload’, ‘waiting’, ‘queued’, ‘running’, ‘ok’, ‘error’, ‘paused’, ‘deleted’, ‘deleted_new’
            create_time: this.string().nullable(),
            model_class: this.string("Job"),

            //ORM Only
            history: this.belongsTo(History, 'history_id'),
            workflow_invocation_step_id: this.string().nullable(),
        }
    }

    //Vuex ORM Axios Config
    static methodConf = {
        http: {
            url: '/api/jobs'
        },
        methods: {
            $fetch: {
                name: 'fetch',
                http: {
                    url: '',
                    method: 'get',
                },
            },
            $get: {
                name: 'get',
                http: {
                    url: '/:id',
                    method: 'get',
                },
            },
            //$create: {
            //    name: 'create',
            //    http: {
            //        url: '', //TODO
            //        method: 'post',
            //    },
            //},
            //$update: {
            //    name: 'update',
            //    http: {
            //        url: '/:id', //TODO
            //        method: 'put',
            //    },
            //},
            //$delete: {
            //    name: 'delete',
            //    http: {
            //        url: '/:id', //TODO
            //        method: 'delete',
            //    },
            //},
        }
    }
}

const Module = {
    ...Common.Module,
    state: Common.State({

    }),
    mutations: {
        ...Common.Mutations,
    },
    actions: {
        ...Common.Actions,
    },
    getters: {
        ...Common.Getters,
    },
};

function register(database) {
    database.register(Job, Module);
}

export {
    Job,
    Module,
    register,
};