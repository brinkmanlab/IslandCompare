import * as Common from "./_common";



class Model extends Common.Model {
    static entity = 'users';
    static primaryKey = 'id';

    static fields() {
        return {
            ...super.fields(),
            id: this.string(null).nullable(),
            username: this.string(null).nullable(),
            quota_percent: this.number(null).nullable(),
            preferences: this.attr({}),
            total_disk_usage: this.number(0),
            deleted: this.boolean(false),
            purged: this.boolean(false),
            nice_total_disk_usage: this.string("None"),
            quota: this.string("unlimited"),
            email: this.string(null).nullable(),
            is_admin: this.boolean(false),
            tags_used: this.attr([]),
        }
    }

    static async getCurrent() {
        let response = await this.$get({
            params: {
                id: 'current',
            }
        });
        return this.find(response.id);
    }

    //TODO add methods to operate on rest of users api

    //Vuex ORM Axios Config
    static methodConf = {
        http: {
            url: '/api/users'
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
            $create: {
                name: 'create',
                http: {
                    url: '',
                    method: 'post',
                },
            },
            $update: {
                name: 'update',
                http: {
                    url: '/:id',
                    method: 'put',
                },
            },
            $delete: {
                name: 'delete',
                http: {
                    url: '/:id',
                    method: 'delete',
                },
            },
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
    database.register(Model, Module);
}

export {
    Model,
    Module,
    register,
};